classdef (Sealed) VSVAICreateTrainingData < handle
% VSVAICreateTrainingData creates training data
%
% Usage
%  % Define training data folder
%  trainingDataPath = fullfile( vsv.file.getVSXDir, ['UserData' filesep 'AI'] );
%
%  % Contruct training data object with the training data file path
%  trainingData = vsv.ai.VSVAICreateTrainingData( trainingDataPath );
%
%  % Set number of images to be in training data set
%  numberOfTrainingImages = 2;
%
%  % Save training data in the folder 'AITraining' ** with a specified number of images
%  % for training.
%
%  trainingData.createTrainingData( 'AI', numberOfTrainingImages );
%
%
% Errors
%   - An error will be thrown during creation if the file path not exist.
%   - If the file path exist but does not contain any files (.vrs) another
%     error will be thrown
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    properties(SetAccess=private)

        % File path for acquired data for training @type char
        fullFilePath;

    end

    properties(Access={?utest.ai.VSVAICreateTraining, ?utest.ai.PackageTest})

        % this property manages the user interaction, i.e. selecting the
        % map and the question dialog. The property is private and and can
        % only be set by the testing framework
        userInteraction;
    end

    methods

        function this = VSVAICreateTrainingData(fullFilePath)
        % VSVAICreateTrainingData Constructs class
        %   A training data must be defined based on the training folder.
        %   The training folder should contain some training data. The
        %   trainign data should be .vrs files.
        %
        % @param fullFilePath - @type char the folder to the training data.
        %                       The training results will be saved in the
        %                       same folder in a subfolder defined by the
        %                       name provided to createTrainingData()


            if ~ischar(fullFilePath)
                error('VSVAICreateTrainingData:invalidArgument', ...
                      'Training file path must be a char');
            end
            if ~exist(fullFilePath, 'dir')
                error('VSVAICreateTrainingData:pathDoesNotExist', ...
                      'Given training path does not exist');
            end

            this.fullFilePath    = fullFilePath;
            this.userInteraction = vsv.ai.TrainUserInteraction();
        end

        function createTrainingData( this, name, numberOfTrainingFiles, doRamdom )
        % createTrainingData Picks the given number of random samples and
        % creates a training data set with the name 'name'
        %
        %  This will take randomly sample data (#numberOfTrainingFiles) in
        %  the folder defined by this.fullFilePath and create training data
        %  from it.
        %
        % @param name - @type char the name for the training data, this will
        %                     define the destination folder for storing the
        %                     training results. The destination folder is
        %                     this.fullFilePath/name
        % @param numberOfTrainingFiles - @type numeric


            import vsv.ai.VSVAIPlotClassifications

            if nargin < 4
                doRamdom = true;
            end

            if ~ischar(name)
                error('VSVAICreateTrainingData:createTrainingData:invalidArgument', ...
                      'Name must be a char giving the training folder name');
            end

            if ~isnumeric(numberOfTrainingFiles) || ~isscalar(numberOfTrainingFiles)
                error('VSVAICreateTrainingData:createTrainingData:invalidArgument', ...
                      'numberOfTrainingFiles must be a numeric scalar');
            end

            % list the files that
            ls    = vsv.file.FileTools.listFilesRegexp(this.fullFilePath, '\.vrs');

            if isempty(ls)
                error('VSVAICreateTrainingData:createTrainingData:noTrainingImages', ...
                      'There are no training images in the provided folder');
            end

            numFiles = length(ls);
            if doRamdom
                index    = randperm( numFiles, min( [ numFiles numberOfTrainingFiles ]));
            else
                index    = 1:min( [ numFiles numberOfTrainingFiles ]);
            end

            % create the training data folder if not exist
            folder = this.createTrainingDataFolder( name );

            % binary files
            bin = vsv.file.VRSBinaryFile;
            fi  = figure();

            % go over selected randomly selected files
            nIndex = length(index);
            for i = 1:nIndex

                filename = ls{ index(i) };
                file     = [this.fullFilePath filesep filename];
                this.userInteraction.currentFile = file;

                bin.readFile( file );

                img = bin.formatData();

                if isvalid(fi)
                    figure(fi);
                else
                    fi = figure;
                end

                clf('reset');
                imagesc( log(double( img ) ));
                colormap('gray');
                ax = gca( fi );

                answ = 'Yes';
                map  = false(size(img));

                while strcmp(answ, 'Yes')

                    BW = this.userInteraction.selectMap();

                    map = map | BW;
                    bn  = bwboundaries(map);
                    hold on;

                    VSVAIPlotClassifications.plotBoundaries(ax, bn, 'y' );

                    answ = this.userInteraction.questionRepeat();

                end

                [~, fn, ~] = fileparts(filename);
                save( fullfile(folder, [ fn '.mat' ] ), 'img', 'map' );
            end
        end

    end

    methods(Access=private)

        function folder = createTrainingDataFolder( this, name )
        % createTrainingDataFolder checks to see if given folder (name)
        % exists. If it doesn't exist then it will create the folder.
        %
        % @param name - @type char
        %
        % @return folder - @type char the folder that was created

            folder = [ this.fullFilePath filesep name 'Training' ];

            if ~exist( folder, 'dir' )
                mkdir( folder );
            end
        end

    end
end

