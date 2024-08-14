classdef (Sealed) VSVAIParameters < handle
% VSVAIParameters A container set of parameters that are used for
% configuring the AI process
%
% @type vsv.ai.VSVAIParameters
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    properties(SetAccess=private)

        % A valid file path to a mat file containing a model for
        % classifying images, @type char
        modelPath = '';

        % Size of image cropped out of image for feature extraction @type
        % numeric, this determines the image-patch size
        kernel;

        % Number of pixels to move kernel through image @type numeric, step
        % size for overlapping image-patches
        stepSize;

        % Threshold for positive classification for simplification of
        % restructed map data @type numeric
        thresholdSimplifyGroupData;

        % Threshold for how much of object to use for training classifier
        % @type numeric
        thresholdPartialClassification;

        % Rows of full scale image to classify @type numeric,
        % @default single(1:size(image, 1))
        rowsForClassification;

        % Columns of full scale image to classify @type numeric
        % @default single(1:size(image, 2))
        columnsForClassification;

        % Number of features that sparse filter will extract @type numeric
        sparseFilteringFeatures;

        % Image rescaling factor @type numeric
        imageResizeFactor;

        % Filter to use on image data @type vsv.ai.imgFilter
        imageFilter;

        % Type of classification model to train
        % @type classreg.learning.classif.FullClassificationModel
        modelType;

        % Class used to extract the features from image data
        % @type vsv.ai.VSVFeatureExtraction
        featureExtraction;

        % True if models should be converted to compact models (do not save
        % training data) @type logical, @default true
        useCompactModels = true;

    end



    methods

        function this = VSVAIParameters(varargin)
        % VSVAIParameters Constructs model path from the file path
        %
        %  initialize the parameter with string value pairs
        %
        %  % Define a training parmaeters
        %  parameters = vsv.ai.VSVAIParameters('kernel', ones( 180, 80 ),...
        %                                      'stepSize', 1, ...
        %                                      'thresholdSimplifyGroupData', 0.25, ...
        %                                      'thresholdPartialClassification', 0.50, ...
        %                                      'rowsForClassification', [], ...
        %                                      'columnsForClassification', [], ...
        %                                      'sparseFilteringFeatures', 150, ...
        %                                      'imageResizeFactor', 0.125, ...
        %                                      'imageFilter', vsv.ai.imgfilter.BtwImageFilter, ...
        %                                      'modelType', 'ClassificationDiscriminant');
        %
        % @param String-value pairs of parameter sets
        %

            this.setParameters(varargin{:});
            this.initializeFeatureExtraction;
        end

        function [desiredRows, desiredCols, bn1] = desiredRowsCols(this, image)
        % desiredRowsCols Determines the appropriate rows and columns for
        % classification
        %
        %  the user can define a set of rows and columns in the image which
        %  are going to be used for classification. Other areas in the
        %  image are ignored. This can be used to speed up live
        %  classification.
        %
        %  If the user did not specified the  rowsForClassification and
        %  columnsForClassification the disiredRows and desiredCols will be
        %  single(1:size(image, 1)); and single(1:size(image, 2));
        %  respectively
        %
        % @param image - @type numeric the image that is going to be
        %                classified
        %
        % @return desiredRows - @type numeric, the rows as a vector giving
        %                       the rows in the image that are supposed to
        %                       be processed for classification
        % @return desiredCols - @type numeric, the columns as a vector
        %                       giving the columns in the image that are
        %                       supposed to be processed for classification

            if ~isempty(this.rowsForClassification) && ~isempty(this.columnsForClassification)
                desiredRows = single(this.rowsForClassification);
                desiredCols = single(this.columnsForClassification);
                % Define corners of search area for object
                fR          = desiredRows(1);
                lR          = desiredRows(end);
                fC          = desiredCols(1);
                lC          = desiredCols(end);
                % Define search area boundaries will be used for plotting
                b1          = [fR fC; lR fC; lR lC; fR lC; fR fC];
                bn1         = {b1};
            else
                % Define search area parameters
                desiredRows = single(1:size(image, 1));
                desiredCols = single(1:size(image, 2));
                % Define search area
                bn1 = {};
            end

        end

        function resizedObject = properResize(this, object)
        % properResize resizes image based on given imageResizeFactor
        %
        % The function uses imresize with the underlying imageResizeFactor
        % and nearest interpolation
        %
        % @param object - @type numeric, the object to be resized
        %
        % @return resizedObject - @type numeric, the resized object
        %

            resizedObject = imresize(object, this.imageResizeFactor, 'nearest');
        end

        function fullSizeObject = undoResize(this, object)
        % undoResize Takes a resized image and correctly resizes it back to
        % full size
        %
        %  use this to resize an image back to its original size that was
        %  resized before using the properResize function
        %
        % @param object - @type numeric, the image object
        % @return fullSizeObject - @type numeric, the resized image object

            fullSizeObject = imresize(object, 1/this.imageResizeFactor, 'nearest');
        end


    end

    methods(Access=private)

        function this = setParameterNameValue(this, name, value)
        % setParameterNameValue sets the properties of the class with the
        % appropriate names and values
        %
        % @param name @type char
        % @param value @type numeric

            switch name
                case 'modelPath'
                    this.setModelPath(value);
                case 'kernel'
                    this.kernel = value;
                case 'stepSize'
                    this.stepSize = value;
                case 'thresholdSimplifyGroupData'
                    this.thresholdSimplifyGroupData = value;
                case 'thresholdPartialClassification'
                    this.thresholdPartialClassification = value;
                case 'rowsForClassification'
                    this.rowsForClassification = value;
                case 'columnsForClassification'
                    this.columnsForClassification = value;
                case 'sparseFilteringFeatures'
                    this.sparseFilteringFeatures = value;
                case 'imageResizeFactor'
                    this.imageResizeFactor = value;
                case 'imageFilter'
                    this.imageFilter = value;
                case 'modelType'
                    this.modelType = value;
                otherwise
                    error( 'VSVAIParameters:setParameterNameValue:unknownInput', ...
                           'Unknown input parameter');
            end

        end

        function setParameters(this, varargin)
        % setParameters Sets properties from given inputs
        %
        % @param varargin @type alternation str and numeric
        %
            if nargin == 1 && isstruct(varargin{1})
                names  = fieldnames(varargin{1});
                values = struct2cell(varargin{1});
            else
                names  = varargin(1:2:end);
                values = varargin(2:2:end);
            end

            lengthValues = length(values);
            lengthNames  = length(names);

            if lengthValues ~= lengthNames
                error('Length of input names not equal to number of inputs')
            end

            for i = 1:lengthValues
                this.setParameterNameValue(names{i}, values{i});
            end

            this.properKernelSize;
        end

        function properKernelSize(this)
        % properKernelSize Determines the proper kernelSize if the
        % imageResize factor is less than 1
        %

            if this.imageResizeFactor < 1
                kernelTmp   = this.properResize(this.kernel);
                this.kernel = kernelTmp;
            end

        end


        function setModelPath(this, modelPath)
        % setModelPath sets model path given by user
        %
        % @param modelPath @type charing
            if ischar(modelPath) && exist(modelPath, 'file')
                this.modelPath = modelPath;
            else
                error('AIParameters:setModelPath:invalidFolder', ...
                      'Model path must be a char providing a valid file');
            end
        end

        function initializeFeatureExtraction(this)
        % intializeFeatureExtraction creates feature extraction class to be
        % used in training
        %
        % @type vsv.ai.VSVAIFeatureExtraction

            this.featureExtraction = vsv.ai.VSVAIFeatureExtraction(this.kernel,...
                                                                   this.stepSize);
        end

    end


end
