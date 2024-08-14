classdef (Sealed) VSVClassifier < handle
% VSVAICLASSIFIER classifies image data and displays the results
%
%  The vsv classifier uses a VSVAIModel to perform the classification and
%  the classifier needs to be constructed by providing such a model. The
%  model performs the actual classification and scoring of the image
%  patches where the classifier provides extra functions to perform image
%  display of the classification results.
%
%  see vsv.ai.VSVAIModel for more information about the actual
%  classification.
%
%  To construct a classifier use vsv.ai.VSVClassifier(model) where model
%  must be a vsv.ai.VSVAIModel. You can also use the load function:
%
%   [classifier, msg] = vsv.ai.VSVClassifier.loadClassifier(file);
%
%  To perform classification use
%  [classLabel, score, features] = classifier.classifyImage(image);
%
%  To display the classification you can use
%  classifier.displayImageClassification(image, mapThresh);
%
% @type vsv.ai.VSVClassifier
%
% Version 1.0 | 2019-10-24
% $Author: Ellyn Cashdollar & Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    properties(SetAccess=private)

        % Model for classification @type vsv.ai.VSVAIModel
        classificationModel;

    end

    properties(Access=private, Transient)

        % Figure handle @type figure
        fh;

        % Handle for classification plot @type graphics handle
        classificationPlotHandle;

        % Handle for detected objects plot @type graphics handle
        boundaryPlotHandle;

        % Handle for current image @type graphics handle
        imageHandle;

    end

    methods(Static)

        function [ vsvClassifier, msg ] = loadClassifier(file)
        % load and construct a classifier from a given file
        %
        %   use this function by using the class name
        %   [ vsvClassifier, msg ] = vsv.ai.VSVClassifier.loadClassifier();
        %
        % @return vsvClassifier - @type vsv.ai.VSVClassifier the classifier
        %                         constructed or [] if the classifer could
        %                         not be located
        % @return msg - @type char, a message which is empty if the
        %               classifier was loaded successfully, or providing
        %               some notes about what caused the classifier not to
        %               load. Reasons can be the file was not found, the
        %               user aborted file selection, the file does not
        %               contain a field called model, or the model in the
        %               mat file is not a vsv.ai.VSVAIModel

            vsvClassifier = [];
            msg = 'Classifier not selected';

            if nargin < 1
                file =[];
            end

            fileinput = file;
            file = vsv.file.FileTools.obtainFilename(file);
            if ~vsv.util.isNULL(file)
                if exist(file, 'file')
                    fc = load(file);
                    if isfield(fc, 'model')
                        model = fc.model;
                        if isa(model, 'vsv.ai.VSVAIModel' )
                            vsvClassifier = vsv.ai.VSVClassifier(model);
                            msg = '';
                        else
                            msg = 'Model in file is not an vsv.ai.VSVAIModel';
                        end
                    else
                        msg = 'File does not contain a model object';
                    end
                else
                    msg = 'File does not exist';
                end
            else
                if ischar( fileinput ) && ~exist(fileinput, 'file')
                    msg = 'File does not exist';
                end

            end

        end

    end

    methods(Access=public)

        function this = VSVClassifier(model)
        %VSVAICLASSIFIER construct an instance of this class
        %   Detailed explanation goes here
        %
        % @param model - @type vsv.ai.VSVAIModel

            this.setAIParameters(model);
        end


        function fh = displayImageClassification(this, results, image)
        % displayImageClassification Classify object trained for in given
        % image
        %
        % @param results - @type struct, containing the classification
        %                   results, use classifyImage to obtain the
        %                   classifciation results. The results struct
        %                   should at least contain a field classBoundary,
        %                   and classificationRegion



            bn  = results.classBoundary;
            bn1 = results.classificationRegion;
            fh  = this.plotImageClassification(bn, bn1, image);

        end

        function results = classifyImage(this, image, mapThresh)
        % classifyImage uses a trained classifier to detect and identify
        % objects in the current ultrasound image.
        %
        % The classification results:
        %       'scoreMapThresholded',  the map displaying the scores for
        %                               each pixel, the map is thresholded
        %                               from scoreMapFull using mapThreshold
        %       'scoreMapFull',         the map containing the scores of
        %                               the classification
        %       'classBoundary',        the boundaries highlighting the
        %                               classification
        %       'classificationRegion', the region that was selected for
        %                               classification, which is defined by
        %                               desired rows and columns
        %
        % @param image - @type numeric, the image to classify the objects in
        %
        % @return results - @type struct, containing the classification
        %                   results
            results = this.classificationModel.createClassificationResult(image, mapThresh);
        end


    end

    methods(Access=private)

        function fh = plotImageClassification(this, bn, bn1, image)
        % plotImage plot classification results on input image
        %
        % @param bn @type cell
        % @param bn1 @type cell
        % @param image @type numeric

            import vsv.ai.VSVAIPlotClassifications

            if isempty(this.fh) || ~isvalid(this.fh)
                this.fh = figure;
            end

            ax = gca(this.fh);

            if isempty(this.imageHandle) || ~isvalid(this.imageHandle)
                this.imageHandle = imagesc( log(image));
                colormap(ax, 'gray');
            else
                this.imageHandle.CData = log(image);
            end

            if ~isempty(this.classificationPlotHandle)
                ph = this.classificationPlotHandle;
                if iscell(ph )
                    nHandles = length(ph);
                    for i = 1:nHandles
                        delete(this.classificationPlotHandle{i});
                    end
                else
                    delete(ph);
                end
                this.classificationPlotHandle = [];
            end

            this.classificationPlotHandle = VSVAIPlotClassifications.plotBoundaries(ax, bn, 'y');

            if ~isempty( bn1 ) == 1
                if isempty( this.boundaryPlotHandle )
                    this.boundaryPlotHandle = VSVAIPlotClassifications.plotBoundaries(ax, bn1, 'c');
                end
            end
            fh = this.fh;
        end

        function setAIParameters(this, model)
        % setAIParameters Set AIParameters for the class
        %
        % @param model @type vsv.ai.VSVAIModel

            if isa( model, 'vsv.ai.VSVAIModel' )
                this.classificationModel = model;
            else
                error('VSVClassifier:invalidArgument', ...
                       'model must be a vsv.ai.VSVAIModel');
            end

        end

    end

end

