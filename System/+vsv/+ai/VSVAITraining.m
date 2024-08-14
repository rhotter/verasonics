classdef VSVAITraining < handle
% VSVAITraining trains classifier with given training data
%
% @type vsv.ai.VSVAITraining
%
% Version 1.0 | 2019-10-28
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    properties(SetAccess=private)

       % Parameters for training @type vsv.ai.VSVAIParameters
       parameters;

    end

    methods

        function this = VSVAITraining(parameters)
        % VSVAITraining sets parameters for training classifier
        %
        % @param trainingParameter @type struct

            this.setTrainingParameters(parameters);

        end

        function vsvAIModel = trainClassifier(this, filePath, verbosity)
        % trainClassifier trains the classifier with training data in the
        % given folder
        %
        %   This is the main function to a classification model and a
        %   sparse filtering model to be used in the classification
        %   process. The first input is the folder where the trainingData
        %   is located (trainingData). The second input controls whether
        %   the progress of training is displayed in the command window.
        %
        %   the function will automatically search for mat files in the
        %   training folder. those mat files should be files that were
        %   created using vsv.ai.VSVAICreateTrainingData and the
        %   createTrainingData function. the mat files will contain a image
        %   data (img) and a map (map) that is a logical map that indicates
        %   the objects in the image that are the targets.
        %
        %   the function may throw some errors if the files cannot be found
        %   or if there no valid files in the specified folder
        %
        % @param filePath - @type char, the file path where the training
        %                   data is located
        % @param verbosity - @type char, @optional, @default 'on', set to
        %                    'off' to limit warning messages
        %
        % @return vsvAIModel - @type vsv.ai.VSVAIModel, the trained
        %                     classificaiton model consisting of the
        %                     classificationModel and the sparseFiltering
        %                     model


            % Trains classifiers in the open folder and creates the model to be used to identify the features
            tic;

            if nargin < 3
                verbosity = 'on';
            end

            givenFileNames = vsv.file.FileTools.listFilesRegexp(filePath, '\.mat');
            if isempty(givenFileNames)
                error( 'trainClassifier:noTrainignData' ,...
                       'No training data found');
            end

            % get the training files
            trainingFileNames = this.checkTrainingFiles(givenFileNames, filePath);

            % @ToDo check whether its empty an throw an error if so
            st                            = load(trainingFileNames{1});
            image                         = single(st.img);
            [desiredRows, desiredCols, ~] = this.parameters.desiredRowsCols(image);

            imgResize    = this.parameters.properResize(single(image));
            [ks, si]     = this.parameters.featureExtraction.kernelImageSize(imgResize);
            [rows, cols] = this.parameters.featureExtraction.rowsCols(ks, si);
            fDSize       = size(zeros(rows*cols, ks(1)*ks(2)));
            featureData  = zeros(fDSize(1)*length(trainingFileNames), fDSize(2), class(image));
            groupData    = zeros(fDSize(1)*length(trainingFileNames), fDSize(2), class(image));

            index = 1;
            for i = 1:length(trainingFileNames)

                this.displayCounterControl(verbosity, trainingFileNames, i);
                st                  = load(trainingFileNames{i});
                [features, map, si] = this.prepareFeatureClassify(single(st.img), single(st.map),...
                                                                  desiredRows, desiredCols);
                featureData(index:index + si(1) - 1, :) = features;
                groupData(index:index + si(1) - 1, :)   = map;
                index                                   = index + si(1);
            end
%
            % Simplify group data
            [groupData, rowAverage] = this.simplifyGroupData(groupData);
            
            status = warning('query', 'stats:classreg:learning:fsutils:Solver:LBFGSUnableToConverge');
            warning('off', 'stats:classreg:learning:fsutils:Solver:LBFGSUnableToConverge');
            sparseFilteringModel    = sparsefilt(featureData, ...
                                              this.parameters.sparseFilteringFeatures, ...
                                              'IterationLimit', 10);
            warning(status.state, 'stats:classreg:learning:fsutils:Solver:LBFGSUnableToConverge');

            % Extract feature data from input data
            featureData = this.parameters.featureExtraction.transformSparse(sparseFilteringModel, featureData);

            % Define parts of the data to exclude from training
            excl                 = this.dataToExcludeInTraining(rowAverage);
            featureData(excl, :) = [];
            groupData(excl)      = [];

            % Train models
            classificationModel = this.trainModel(featureData, groupData);

            vsvAIModel = this.createModel( classificationModel, sparseFilteringModel);

            % make the model compact
            if this.parameters.useCompactModels
                 vsvAIModel.makeModelCompact();
            end
        end


    end

    methods(Access=private)

        function model = createModel(this, classificationModel, sparseFilteringModel)
        % createModel creates the model that contains the classification
        % model as well as any other information needed for classification
        %
        % @param classificationModel @type vsv.ai.VSVAIModel
        % @return model - @type vsv.ai.VSVAIModel the model to be created
        %

            model = vsv.ai.VSVAIModel(classificationModel, sparseFilteringModel, ...
                                      this.parameters);
        end

        function excl = dataToExcludeInTraining(this, rowAverage)
        % dataToExcludeInTraining Determines the data to be excluded from
        % training based on given
        % thresholdPartialClassification
        %
        % @param rowAverage - @type numeric, the rowAverage for removal
        %
        % @return excl - @type logical, true if sample should be excluded,
        %                false otherwise

            excl = (rowAverage > 0) & (rowAverage <= this.parameters.thresholdPartialClassification);
        end

        function saveTrainedClassifier(this, modelFileName, trainingDataFilePath, verbosity)
        % saveTrainedClassifier Saves a trained classifier and all of the
        % parameters needed for classification in a struct
        %
        % @param modelName @type string

            [classificationModel, sparseFilteringModel] = this.trainClassifier(trainingDataFilePath, verbosity);
            model = this.createModel(classificationModel, sparseFilteringModel);

            [fp,fn,~] = fileparts(modelFileName);
            save(fullfile(fp, [fn '.mat']), 'model')
        end

        function [data, rowAverage] = simplifyGroupData(this, groupData)
        % simplifyGroupData takes the given group data and takes the mean of the
        % data in each row. If this mean is less than the given threshold
        % the value of 0 is given to the row, if greater than one the
        % values of 1 is given to the row. This is done to not confuse the
        % classifier in training.
        %
        % @param groupData @type numeric
        %
        %

            groupDataSize = size(groupData);
            rowAverage    = zeros(groupDataSize(1), 1, class(groupData));
            data          = zeros(groupDataSize(1), 1, class(groupData));
            for row = 1:groupDataSize(1)
                rowAverage(row, 1) = mean(groupData( row, :));
            end
            thresholdValue                     = this.parameters.thresholdSimplifyGroupData * max(rowAverage);
            data(rowAverage >= thresholdValue) = 1;
            data(rowAverage < thresholdValue)  = 0;

        end

        function [features, groups, si] =  prepareFeatureClassify(this, image, map, desiredRows, desiredCols)
        % prepareFeatureClassify Prepares the input image for
        % classification by extracting the features from it
        %
        % @param image - @type numeric the input image
        % @param map - @type logical the map that indicates the target
        %               objects, the map is true for all pixels that are
        %               part of a target object
        % @param desiredRows - @type numeric the rows to be processed
        % @param desiredCols - @type numeric the columns to be processed
        %
        % @return features - @type numeric the feature vector
        % @return groups - @type numeric, the groups that assign a group to
        %                   each feature vector
        % @return si - @type numeric, the size of the features
            filtImg   = this.parameters.imageFilter.filterImage(image(desiredRows, desiredCols));
            imgResize = this.parameters.properResize(filtImg);
            mapResize = this.parameters.properResize(map(desiredRows, desiredCols));
            features  = this.parameters.featureExtraction.imagePatches(imgResize);
            si        = size(features);
            groups    = this.parameters.featureExtraction.imagePatches(single(mapResize));

        end

        function displayCounterControl(this, verbosity, trainingFileNames, i)
        % displayCounterControl Controls whether the display displays
        % current progress
        %
        % @param verbosity @type char

            switch verbosity
                case 'on'
                    disp([ 'i = ' num2str(i) ' of ' num2str(length(trainingFileNames))]);
                case 'off'

                otherwise
                    error('Incorrect input. Correct inputs are ''on'' or ''off''.');
            end

        end

        function fileName = isTrainingDataValid(~, fullFileName)
        % isTrainingDataValid Determines whether provided training data is
        % correct
        %
        % @param fullFileName @type char
        % @return fileName - @type char the file name if valid, {}
        %                    otherwise

            file = matfile(fullFileName);
            fileDetails = whos(file);
            filePartNames = {fileDetails.name};
            if strcmp(filePartNames, {'img', 'map'})
                fileName = fullFileName;
            else
                fileName = '';
            end
        end


        function trainingFiles = checkTrainingFiles(this, givenFileNames, filePath)
        % checkTrainingFiles Check files provided for training to see if
        % they are correct
        %
        % @param givenFileName @type str cell
        % @param filePath @type char where to check the files
        %
        % @return trainingFiles - @type cell array of valid files

            nFiles = length(givenFileNames);
            trainingFiles = cell(nFiles, 1);

            for counter = 1:nFiles
                fileName = this.isTrainingDataValid([filePath filesep givenFileNames{counter}]);
                trainingFiles{counter} = fileName;
            end

            trainingFiles = trainingFiles(~cellfun( 'isempty', trainingFiles));

            if isempty(trainingFiles)
                error( 'checkTrainingFiles:noFiles' , ...
                       'No training data found');
            end
        end

        function classificationModel = trainModel(this, featureData, groupData)
        % trainModel Trains the classification model. This model type is
        % determined based on the user model type input parameters
        %
        % The function is using the matlab machine learning toolbox to
        % train a classification model, which then can be used for
        % classification. See the machine learning toolbox for more details
        %
        % @param featureData - @type numeric, the feature vector, rows are
        %                     the samples and cols are the features for
        %                     each sample
        % @param groupData - @type numeric the grouping vector that assigns
        %                    a class (group) for each
        %
        % @return classificationModel - the classification model returned
        %                   either from fitcsvm, fitcknn, fitcnb, fitctree,
        %                   fitcdiscr
        %

            switch this.parameters.modelType
                case 'ClassificationSVM'
                    classificationModel = fitcsvm(featureData, groupData, 'KernelScale','auto', ...
                                                               'Standardize', true,  ...
                                                               'Solver',      'SMO',  ...
                                                               'KernelFunction','rbf',  ...
                                                               'BoxConstraint', 10);
                case 'ClassificationKNN'
                    classificationModel = fitcknn(featureData, groupData, 'NumNeighbors', 5, 'Standardize', 1);
                case 'ClassificationNaiveBayes'
                    classificationModel = fitcnb(featureData, groupData);
                case 'ClassificationTree'
                    classificationModel = fitctree(featureData, groupData);
                case 'ClassificationDiscriminant'
                    classificationModel = fitcdiscr(featureData, groupData);
            end

        end

        function setTrainingParameters(this, parameters)
        % Extracts and save training parameters from input
        %
        % @param trainingParameters @type vsv.ai.VSVAIParameters

            if isa(parameters, 'vsv.ai.VSVAIParameters')
                this.parameters = parameters;
            else
                error( 'VSVAITraining:setTrainingParameters:invalidArgument', ...
                       'Parameters must be a vsv.ai.VSVAIParameters' );
            end

        end
    end


end
