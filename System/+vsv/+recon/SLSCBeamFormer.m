classdef SLSCBeamFormer < vsv.recon.AbstractBeamFormer
% SLSCBeamforming 
% This class calculates an image using Short-Lag Spatial Coherence (SLSC).
% SLSC is a beamforming technique that uses coherence to calculate an 
% image. More information on this beamforming method can be found in the
% link below.
%
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3172134/ 
%
% The class SLSCBeamFormer was developed for use in the
% SetUpC5_2vWideBeam_ExtBeam example script, however it can be applied to 
% other example scripts as well.
%
% Version 1.0 | 2021-04-12
% Copyright 2001-2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    properties
        % Number of elements in the receive aperture
        NumElementsRcv (1,1) {mustBeNumeric,mustBeFinite,mustBeNonempty} = 128;
        
        % Short-lag values to use in coherence calculations in SLSC
        % beamforming. Specifies all of the distances between elements (in
        % number of elements) in spatial coherence calculations.
        ShortLagVals (1,:) {mustBeNumeric,mustBeFinite,mustBeNonempty} = (1:10);
    end
    
    properties(Dependent)
        % The maximum short-lag value used in the SLSC beamforming. 
        % Specifies the maximum distance between elements (in 
        % number of elements) in spatial coherence calculations.
        MaxShortLag (1,1){mustBeNumeric,mustBeFinite,mustBeNonempty};
    end
    
    % Public functions
    methods
        function this = SLSCBeamFormer(varargin)
        % SLSCBeamFormer This function sets the properties of the
        % class. 
        % @param varargin: A set of string labels and their corresponding 
        %                  values to be set as properties of the class. 
        %                  varargin = {'Param1',value,'Param2',value}
        % Examples:
        % * beamformer = vsv.recon.SLSCBeamFormer('NumElementsRcv',128,'MaxShortLag',7);
        % * beamformer = vsv.recon.SLSCBeamFormer('MaxShortLag',7);
        % * beamformer = vsv.recon.SLSCBeamFormer('NumElementsRcv',128,'ShortLagVals',(1:13));
        %  
            inputParameters = varargin;            
            this.setClassParameters(inputParameters);
        end
        
        function set.MaxShortLag(this,val)
            this.ShortLagVals = (1:val);
        end
        
        function val = get.MaxShortLag(this)
            val = this.ShortLagVals(end);
        end
        
        function Q = calculateQ(this)
        % calculateQ This function calculates the value of Q for
        % evaluation of the spatial coherence beamforming function. Q is
        % equal to the ratio of the maximum short-lag distance used (m) to 
        % the number of elements in the receive aperture (N). Mathematically,
        % it can be written as follows:
        %
        %   Q = m/N
        %
        % @result Q: @type numeric value used to evaluate the spatial
        %             coherence beamformer (Q is in %)
        %
            Q = this.MaxShortLag./this.NumElementsRcv * 100;
        end
    end
    
    % External Beamforming functions
    methods(Access=protected)
        function SLSCImage = computeImageProtected(this, signal1,~)
        % computeImageProtected computes the SLSC image from the IData,
        % inputted to the function and MaxShortLag defined in the class
        % @param signal1: @type numeric IData output from the system for each
        % channel
        % @result image: @type numeric Resulting SLSC image calculating
        %
            arguments
                this
                signal1 (:,:,:,:) {mustBeNumeric}
                ~
            end

            % Calculate the second signal used in SLSC calculations
            signal2 = this.calculateSignal2(signal1);
            
            % Calculate the channel scaling factor for the SLSC image
            channelScaling = this.calculateChannelScalingFactor(signal1);
            
            % Calculate the spatial coherence between signal1 and signal2
            R = this.calculateCoherence(signal1,signal2);
            
            % Calculate the resulting SLSC image from the spatial coherence
            % calculated between signal1 and signal2
            ImgData = this.calculateSLSCImage(R);
            
            % Scale image based on the number of channels that provided
            % data used in the reconstruction
            scaledImage = ImgData .* channelScaling;
            
            % Adjust the data for display
            SLSCImage = scaledImage./max(abs(scaledImage),[],'all') .*256;
            SLSCImage(SLSCImage <0) = 0;
        end
        
        function scalingFactor = calculateChannelScalingFactor(~,signal1)
        % calculateChannelScalingFactor Calculates the appropriate scaling 
        % factor for the SLSC image so that it displays correctly.
        %
            channelScale = sum(signal1 ~= 0,4);
            channelScale(channelScale == 0) = 1;
            scalingFactor =  1./channelScale; %magnitudeScale .* maxValScale;
        end
        
        function signal2 = calculateSignal2(this,signal1)
        % calculateSignal2 Calculates a second signal to use in coherence
        % calculations. The variable signal2 is signal1 shifted by 1 to m 
        % short-lags, where m is the maximum short-lag used. This  
        % information is stored in a five-dimensional  matrix. The
        % data shifted by the appropriate short-lag value (1:m) 
        % is concatenated along the fifth dimension of the matrix. This
        % results in the size signal2 being [size(signal1),m].
        %
            % Calculate the size of the data
            sizeOfData = size(signal1); 
            
            % Calculate the second signal used in SLSC calculations
            signal2 = zeros([sizeOfData,length(this.ShortLagVals)]); % Initialize signal2
            for m = this.ShortLagVals                                                                                                                                                                                                 
                % Define the different signal inputs for later calculations
                signal2(:,:,:,1:end-m,m) = signal1(:,:,:,m+1:end);
            end
        end
        
        function R = calculateCoherence(~,signal1,signal2)
        % calculateCoherence Calculates the coherence between signal1 and
        % signal2.
        %
            R = (signal1 .* signal2) ./ abs(signal1 .* signal2); % Assumes no kernel
            R(isnan(R)) = 0;
        end
        
        function ImgData = calculateSLSCImage(this,R)
        % calculateSLSCImage Calculates the SLSC image from the coherence
        % by weighting the sum of the coherences calculated for each of the
        % specified short-lag values (1:m).
        %
            mVector = reshape(1./(this.NumElementsRcv - this.ShortLagVals),...
                              [1,1,1,1,numel(this.ShortLagVals)]);
            ImgData = sum(mVector .* sum(R,4),5);
        end
                 
    end
    
    % Class setup functions
    methods(Access=protected,Sealed)
        function setClassParameters(this,inputParameters)
        % setClassParameters Sets the appropriate input parameters of the
        % class when initialized
        % @param inputParameters: @type cell Parameters input to the class
        %                         to initialize the class
        %
        
            % Sets the appropriate class parameters
            names = inputParameters(1:2:end);
            values = inputParameters(2:2:end);
            if numel(names) == numel(values)
                this.setupInputParameters(names,values);
            end
        end
        
        function setupInputParameters(this,names,values)
        % setupInputParameters This function makes sure that the correct
        % parameters are set.
        % @param names: @type string Names of the parameters to set
        % @param values: @type numeric Values to set with the
        %                corresponding parameter name
            
            for k = 1:length(names)
                name = names{k};
                value = values{k};
                switch name
                    case 'NumElementsRcv'
                        this.NumElementsRcv = value;
                    case 'MaxShortLag'
                        this.MaxShortLag = value;
                    case 'ShortLagVals'
                        this.ShortLagVals = value;
                    otherwise
                        error('setupInputParameters:IncorrectInput',...
                              'Incorrect inputs. Only input names NumElementsRcv, MaxShortLag, and ShortLagVals are accepted.')
                end
            end
        end
    end
end