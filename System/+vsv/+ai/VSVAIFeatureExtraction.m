classdef (Sealed) VSVAIFeatureExtraction < handle
% VSVAIFeatureExtraction Extracts features from image and reconstructs
% image from map data
%
% @type vsv.ai.VSVAIFeatureExtraction
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.


    properties(SetAccess=private)

        % The kernel providing the size of image cropped out of image for
        % feature extraction @type numeric
        kernel;

        % Number of pixels to move kernelSize through image @type numeric
        stepSize;

    end

    methods

        function this = VSVAIFeatureExtraction(kernel, stepSize)
        % VSVAIFeatureExtraction Contructs @type vsv.ai.VSVAIFeatureExtraction
        % class
        %
        % Several errors may be thrown in case kernel or stepSize are not
        % valid:
        %       - both must be numeric and not empty
        %       - step size must be a scalar and larger than 0
        %
        % @param kernel   - @type numeric kernel must be a 2D matrix
        %                   providing ones where to extract image data for
        %                   features
        % @param stepSize - @type stepSize, a scalar providing how many
        %                   steps the kernel is moving

            if nargin < 2
                stepSize = 1;
            end

            if isnumeric(kernel) && isnumeric(kernel)

                if isscalar(stepSize) && ~isempty(kernel)

                    if stepSize > 0
                        this.kernel = kernel;
                        this.stepSize   = round( stepSize );
                    else
                        error('VSVAIFeatureExtraction:invalidArgument', ....
                            'kernel size and step size must be larger than 0' );
                    end

                else
                    error('VSVAIFeatureExtraction:invalidArgument', ....
                            'kernel size must not be empty and stepSize must be a scalar' );
                end
            else
                error('VSVAIFeatureExtraction:invalidArgument', ...
                      'Kernel size and stepSize must be numeric' );
            end

        end

        function features = imagePatches(this, image)
        % imagePatches Extracts image patches from an image based on a
        % given kernel, step size and image.
        %
        % This extracts kernel sized images. The number of
        % image-patches that are derived from image depends on the kernel
        % size and the step size.
        %
        % The images are then reshaped in to a vector that can be used as
        % feature vector for classification.
        %
        % A single image patch is considered to be a set of NF features
        % where NF = prod( size(image) ). Feature 1 is image-patch(1,1) and
        % feature NF is image-patch( end, end). The features are
        % linearlized based on matlab rows first order.
        %
        % The final features vector consists of #image-patch feature
        % vectors, in other words features = [ #image-patch, NF ].
        %
        %
        % @param image @type numeric
        %
        % @return features - @type numeric, the feature image of size,
        %                    [ rows*cols, NF ] where NF is the number of
        %                    features, default is kernelsize(1)*kernelsize(2)
        %

            [ks, si]     = this.kernelImageSize(image);
            [rows, cols] = this.rowsCols(ks, si);
            ks           = uint32(ks);
            si           = uint32(si);
            NF           = ks(1)*ks(2);
            features     = zeros(rows*cols, NF, class(image));
            indexCols    = 1;
            indexRows    = 1;
            rowStep         = (indexRows:rows-1)'.*this.stepSize;
            rowFeatureIndex = uint32(1:length(rowStep))';

            r            = indexRows:indexRows + ks(1) - 1;

            for ci = 1:cols
                c                                       = indexCols:indexCols + ks(2) - 1;
                [C, R]                                  = meshgrid(c, r);
                imgIndex                                = R + (C-1).*si(1);
                featureImgIndex                         = imgIndex(:)' + rowStep;
                features(rowFeatureIndex + (ci-1).*rows+1, :)   = image(featureImgIndex);
                indexCols                               = indexCols + this.stepSize;
            end

        end

        function map = patchRecon(this, mapData, image)
        % patchRecon Reconstructs classification map from classification
        % score
        %
        % @param mapData - @type numeric, map classification scores
        % @param image - @type numeric, the image
        %
        % @return map - @type numeric, the map that was reconstructed, will
        %                have the same size as input image

            [ks, si]     = this.kernelImageSize(image);
            [rows, cols] = this.rowsCols(ks, si);
             map         = zeros(si, class(image));
             NF          = ks(1)*ks(2);

             indexRows   = 1;
             indexCols   = 1;
             index       = 1;

             % For image reconstruction directly from imagePatches
             if size( mapData, 2 ) == NF
                 for ci = 1:cols

                    ic = indexCols:indexCols + ks(2) - 1;
                    for ri = 1:rows

                        ir = indexRows:indexRows + ks(1) - 1;
                        map( ir, ic) = map(ir, ic) + reshape(mapData(index, :), ks);
                        indexRows    = indexRows + this.stepSize;
                        index        = index + 1;
                    end
                    indexCols = indexCols + this.stepSize;
                    indexRows = 1;
                 end

             else % For reconstruction of map from classLabel
                for ci = 1:cols

                    ic = indexCols:indexCols + ks(2) - 1;
                    for ri = 1:rows

                        ir = indexRows:indexRows + ks(1) - 1;
                        map(ir, ic) = map(ir, ic) + mapData(index);
                        indexRows   = indexRows + this.stepSize;
                        index       = index + 1;
                    end

                    indexCols = indexCols + this.stepSize;
                    indexRows = 1;
                end
             end

        end


        function T = transformSparse(this, sparseFilteringModel, X)%#ok
        % transformSparse Transforms data with sparse filter without
        % normalizing based on the size of the input data
        %
        % @param sparseFilteringModel @type SparseFiltering,
        %               parseFilteringModel created in training the classifiers
        % @param X @type numeric, input data to be transformed by sparse
        %                filter

            % 1. Compute a tentative N-by-Q matrix of new features.
            delta = 1e-8;
            W     = sparseFilteringModel.TransformWeights;
            T     = sqrt(delta + (X*W).^2);

        end

        function [rows,cols] = rowsCols(this, ks, si)
        % returns the number of rows and columns used for feature
        % extraction and classification
        %
        % @param ks - @type numeric the kernel size
        % @param si - @type numeric the image size
        %
        % @return rows - @type numeric, number of rows
        % @return cols - @type numeric, number of columns

            rows = uint32(ceil((si(1) - ks(1) + 1)/this.stepSize));
            cols = uint32(ceil((si(2) - ks(2) + 1)/this.stepSize));
        end

        function [ks,si] = kernelImageSize(this, image)
        % returns the kernel size and size of the image
        %
        % @param image - @type numeric the image size
        %
        % @return ks - @type numeric, the kernel size
        % @return si - @type numeric, the image size

            ks = size(this.kernel);
            si = size(image);
        end

    end


end

