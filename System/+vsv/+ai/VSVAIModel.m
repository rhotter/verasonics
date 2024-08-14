classdef (Sealed) VSVAIModel < handle
% VSVAIModel creates a class that stores the model and all parameters
% needed for classification in one location
%
% @type vsv.ai.VSVAIModel
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.


    properties(SetAccess=private)

       % Sparse filtering model used in classification
       % @type SparseFiltering a class from the machine learning toolbox
       sparseFilteringModel;

       % Model to be used for classification
       % @type classreg.learning.classif.FullClassificationModel
       classificationModel;

       % Parameters for classification
       % @type vsv.ai.VSVAIParameters
       parameters;

    end

    methods

        function this = VSVAIModel(varargin)
        % VSVAIModel constructs @type vsv.ai.VSVAIModel class
        %
        % @param varargin @type vsv.ai.AIParameters, SparseFiltering,
        % classreg.learning.classif.FullClassificationModel

            for k = 1:length(varargin)
                input = varargin{k};
                if isa(input, 'vsv.ai.VSVAIParameters')
                    this.parameters = input;
                elseif isa(input, 'SparseFiltering')
                    this.sparseFilteringModel = input;
                elseif isa(input, 'classreg.learning.classif.FullClassificationModel')
                    this.classificationModel  = input;
                else
                    error( 'VSVAIModel:invalidArgument', ...
                           'Input does not contain correct models');
                end
            end

        end

        function makeModelCompact(this)
        % this will transform the classification model to a compact version
        %
        % The standard models will save all the training data which can
        % take a lot of disk space. If you only interested in using the
        % model for classification and prediction this will help to
        % compress the model data to only the information needed for
        % performing the classification
        %
            this.classificationModel = compact(this.classificationModel);
        end


        function [results] = createClassificationResult(this, ...
                                                        image,...
                                                        mapThresh)
        % scoreImage Returns the classification score and reformats them
        % into structs with boundary information of a region classified
        % to have desired object in it
        %
        % @param image - @type numeric, the image where to find the objects
        %                in
        % @param mapThresh - @type numeric, the map threshold, a threshold to
        %                  exclude all scores < mapThreshold from display

            [desiredRows, desiredCols, bn1] = this.parameters.desiredRowsCols(image);

            imgResize           = this.parameters.properResize(image(desiredRows, desiredCols));
            [dRR, dCR]          = this.desiredRowsColsResize(desiredRows, desiredCols);
            score               = this.classifyResizedImage(imgResize);
            map                 = this.parameters.featureExtraction.patchRecon(score(:, 2), imgResize);
            mapFull             = zeros(size(imgResize));
            mapFull( dRR, dCR ) = map;
            map1                = mapFull >= mapThresh;
            bn                  = bwboundaries(this.parameters.undoResize(map1), 'noholes');

            results = struct( 'scoreMapThresholded',  map1, ...
                              'scoreMapFull',         mapFull, ...
                              'classBoundary',        {bn}, ...
                              'classificationRegion', {bn1});
        end


    end


    methods(Access=private)

        function [desiredRowsResize, desiredColsResize] = desiredRowsColsResize(this, desiredRows, desiredCols)
        % desiredRowsColsResize adjusts the rows and columns for
        % classification for resized images
        %
        % @param desiredRows @type numeric
        % @param desiredCols @type numeric

            imrf = this.parameters.imageResizeFactor;

            desiredRowsResize = ceil(desiredRows(1)*imrf):ceil(desiredRows(end)*imrf);
            desiredColsResize = ceil(desiredCols(1)*imrf):ceil(desiredCols(end)*imrf);
        end

        function score = classifyResizedImage(this, image)
        % classifyImage uses a trained classifier to detect and identify
        % objects in the current ultrasound image.
        %
        % @param image @type numeric

            filtImg             = this.parameters.imageFilter.filterImage(image);
            features            = this.parameters.featureExtraction.imagePatches(single(filtImg));
            featureData         = this.parameters.featureExtraction.transformSparse(this.sparseFilteringModel, features);
            [~, score] = this.classificationModel.predict(featureData);

        end

    end

end
