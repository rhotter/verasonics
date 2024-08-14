%% Demo of classifying objects in ultrasound images

%% Step 1: Acquire data

cd(vsv.file.getVSXDir)
% Set up for saving images
SetUpL11_5vFlashAngles_SaveImages

% Start Vantage Software
VSX

%% Step 2: Training Preparation

% Define training data folder
trainingDataPath = fullfile( vsv.file.getVSXDir, ['UserData' filesep 'AI'] );

% Contruct training data object with the training data file path
trainingData = vsv.ai.VSVAICreateTrainingData( trainingDataPath );

% Set number of images to be in training data set
numberOfTrainingImages = 2;

% Save training data in the folder 'AITraining' ** with a specified number of images
% for training.

trainingData.createTrainingData( 'AI', numberOfTrainingImages );

% Define a training parmaeters
parameters = vsv.ai.VSVAIParameters('kernel', ones( 180, 80 ),...
                                    'stepSize', 1, ...
                                    'thresholdSimplifyGroupData', 0.25, ...
                                    'thresholdPartialClassification', 0.50, ...
                                    'rowsForClassification', [], ...
                                    'columnsForClassification', [], ...
                                    'sparseFilteringFeatures', 150, ...
                                    'imageResizeFactor', 0.125, ...
                                    'imageFilter', vsv.ai.imgfilter.BtwImageFilter, ...
                                    'modelType', 'ClassificationDiscriminant');

%% Step 3: Train classifier

% Construct the classifier training object
training = vsv.ai.VSVAITraining( parameters );

% Train a classification model and a sparse filtering model to be used in
% the classification process. The first input is the folder where the
% trainingData is located (trainingData). The second input controls whether
% the progress of training is displayed in the command window.
model = training.trainClassifier( [trainingDataPath filesep 'AITraining'], 'on' );

% you have access to the the model components
classificationModel  = model.classificationModel;
sparseFilteringModel = model.sparseFilteringModel;

% Save trained classifier model
save( [trainingDataPath filesep 'ClassificationModel.mat'], 'model' );

%% Step 4: Classify Single Image

% Pick image to classify
fileNames = vsv.file.FileTools.listFilesRegexp( trainingDataPath, '\.vrs' );
image = vsv.file.readVRSFile([trainingDataPath filesep fileNames{2}]);
% Construct classifier object

modelPath = fullfile( vsv.file.getVSXDir, ['UserData' filesep 'AI' filesep 'brightCircleModel.mat'] );

% load the classifier from file
[classifier, msg] = vsv.ai.VSVClassifier.loadClassifier(modelPath);
% use the following to use a UI-based selection of the fiel
% [classifier, msg] = vsv.ai.VSVClassifier.loadClassifier(modelPath);

if isempty(classifier)
    error('Classifier is empty. valid file?')
end
% Set map threshold
mapThreshold = 20;

% create classification results
clResults = classifier.classifyImage(image, mapThreshold);

% Display classification of a single image
classifier.displayImageClassification( clResults, image );

subplot(1,2,1, gca);
subplot(1,2,2);
imagesc( clResults.scoreMapFull );


%% Classify Images LIVE!

cd(vsv.file.getVSXDir)
% Set up live classification
SetUpL11_5vFlashAngles_Classify

% Start Vantage software
VSX

