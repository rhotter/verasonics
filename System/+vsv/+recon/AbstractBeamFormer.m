classdef AbstractBeamFormer < handle
% AbstractBeamFormer Base class for external beamforming classes
% Classes that inherit from AbstractBeamFormer must have the function
% computeImageProtected defined in the class. This allows for the
% computeImage function to operate correctly. This abstract class adds
% extra security that the data being processed is the correct type.
%  
% Version 1.0 | 2021-04-12
% Copyright 2001-2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

    methods(Abstract,Access=protected)
    % Defines the function that inheriting classes must have. All
    % inheriting classes must have IData and QData inputs to the function,
    % however it is not required that both be used in reconstruction.
    %
        image = computeImageProtected(this, IData,QData)
    end
    
    methods(Sealed)
        function image = computeImage(this, IData, QData)
        % computeImage This function calculates the image created with the
        % new beamforming method.
        % @param IData: @type numeric Reconstructed IData for each 
        %               transducer element
        % @param QData: @type numeric Reconstructed QData for each 
        %               transducer element
        % @result image: @type numeric Image created from external
        %                beamformer
        
            arguments
                this vsv.recon.AbstractBeamFormer
                IData {mustBeNumeric,mustBeFinite,mustBeNonempty}
                QData {mustBeNumeric,mustBeFinite}
            end
            
            image = this.computeImageProtected(IData,QData);
            mustBeNumeric(image);
            mustBeNonempty(image);
        end
    end
    
end


