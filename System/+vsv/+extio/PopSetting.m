classdef PopSetting < handle
%POPSETTING  settings for the pulse on position operation
%
%   Use this class as part of the Encoder input settings to specify a
%   scanning sequence for pulse on position operation. 'EncNumber'
%   specifies which encoder input should be read. 'StartTick' specifies the
%   encoder tick count to initiate the first Vantage acquisition.
%   'IncrementTick' specifies the number of ticks between acquisitions,
%   which may be negative for backwards position scanning (if negative,
%   'StartTick' will also be made negative). 'NumSteps' specifies the
%   number of position steps in the scanning sequence.
%
%   Example:
%
%   % creates pulse on position profile for encoder 1
%   popSetting = vsv.extio.PopSetting('EncNumber', 1, ...
%                                     'StartTick', 100, ...
%                                     'IncrementTick',100, ...
%                                     'NumSteps' 10;); 
%
% Version 1.0 | 2020-06-16
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    properties
        
        % Number of encoder to use for position must be > 0 
        % @type (1,1) int32  
        EncNumber (1,1) int32 {mustBeGreaterThan(EncNumber, 0)} = 1;
        
        % Measurement starting position in encoder ticks must be >= 0
        % @type (1,1) int32  
        StartTick (1,1) int32 {mustBeGreaterThanOrEqual(StartTick, 0)}  = 100;
        
        % Increment between position in encoder ticks @type (1,1) int32  
        IncrementTick (1,1) int32  = 100;
        
        % Numer of poistion steps (1,1) int32  must be > 0
        NumSteps (1,1) int32 {mustBeGreaterThan(NumSteps, 0)} = 10;
        
        % user definable time out when waiting for position to be reached
        % @type (1,1) int32  must be > 0
        TimeOut (1,1) int32 {mustBeGreaterThan(TimeOut, 0)} = 100;

    end
    
    properties(Hidden)       
        
        %index of current location
        PopIndex (1,1) int32 {mustBeGreaterThan(PopIndex, 0)} = 1;
        
    end
    
    %% public methods
    methods
        
        function this = PopSetting(varargin)
        % create a new pop setting using parameter value pairs as inputs
        %
        %   @signature this = PopSetting(param, value, ....)
        %
        % @Usage
        %
        %   popSetting = vsv.extio.PopSetting('EncNumber', 1,
        %   'StartTick', 100, 'IncrementTick',100, 'NumSteps' 10;); %
        %   creates pulse on position profile for encoder 1
        %
        %   popSetting = vsv.extio.PopSetting( 'EncNumber', 1);
        %
        % @param varargin - parameter value pairs
            
            % will apply the parameter value pairs to this object
            vsv.util.paramValPair(this, varargin{:});
            
        end
        
        function popArray = createPopArray( this )
        % Will create a 1D array of tick values
        %
            
            incrementTick = this.IncrementTick;
            startTick = this.StartTick;
            numSteps = this.NumSteps;
            
            if incrementTick <0
                startTick = startTick*-1;
            end
            
            popArray = startTick:incrementTick:(startTick+incrementTick*numSteps);
        end
        
    end
    
end

