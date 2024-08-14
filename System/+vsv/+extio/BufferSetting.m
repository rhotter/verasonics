classdef BufferSetting < handle
%BUFFERSETTING summarizes settings for the ExtIO result data buffers
%   
%   A ExtIO result-buffer data is defined by its BufferLength, and
%   BufferType.
%
%   Use this class as part of the Encoder, Analog, and Digital input
%   Settings. The resulting buffer will have a size of number of input
%   elements (e.g., number of encoders) times the BufferLength. 
% 
%   If BufferType is 'Fixed' then the buffer will fill in as many counts as
%   BufferLength then stop.  If 'Circular', the buffer will begin to
%   overwrite from the first sample once full.
%
%   Example:
%       
%   buffer = vsv.extio.BufferSetting('BufferType', 'Circular',
%   'BufferLength', 100); % create a circular buffer of size 100
%
%   buffer = vsv.extio.BufferSetting( 'BufferLength', 100); 
%
% Version 1.0 | 2020-06-16 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

% InternalDocumentation:
%
% @ToDo:
% - check whether int32 is sufficient for buffer length

    properties
        
        % the type of the buffer which can be "Circular" or "Fixed",
        % default is "Circular", @type string
        BufferType (1,:) string ...
                {mustBeMember(BufferType, ["Circular", "Fixed"] ) } = "Circular";
        
        % The size of the buffer, @default = 0
        % note if a double value is set this will be converted to int32
        BufferLength (1,1) int32 {mustBeGreaterThan(BufferLength, 0)} = 1;
        
        BufferDepth (1,1) int32 {mustBeGreaterThan(BufferDepth, 0)} = 1;

    end
    
    %% public methods
    methods
        
        function this = BufferSetting(varargin)
        % create a new buffer setting using parameter value pairs as inputs
        %
        %   @signature this = BufferSetting(param, value, ....)
        %
        % @Usage
        %
        %   buffer = vsv.extio.BufferSetting('BufferType', 'Circular',
        %   'BufferSize', 100); % create a circular buffer of size 100
        %
        %   buffer = vsv.extio.BufferSetting( 'BufferSize', 100); 
        %
        % @param varargin - parameter value pairs
        
            % will apply the parameter value pairs to this object, this
            % will call this.BufferType = buffer type if one of the
            % parameter value pairs is 'BufferType', buffer
            vsv.util.paramValPair(this, varargin{:});
            
        end
        
        function data = createBufferData( this, numElements)
        % Will create a 2D matrix filled with zeros 
        %
        %   The resulting matrix will be of size numElements X
        %   this.BufferLength
        %
        % @param numElements - @type numeric, the number of elements for
        %                      the buffer
        
            data = zeros( numElements, this.BufferLength);
        end
        
        
    end
    
end

