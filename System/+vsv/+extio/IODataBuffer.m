classdef IODataBuffer < handle
%IODATABUFFER The result data buffer
%  
% @ToDo Jac provide description here
%  The buffer data can be considered as an array of linear data vectors.
%  
%  dimension represents the number of data points ()
%
%   [ DataPoint, NumVectors * NumSamples] = size(this.Data,2) = this.DataSize
%
% Version 1.0 | 2020-11-05 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
    properties
        
        % Type of the buffer, either "Circular" or "Fixed". A circular
        % buffer means that data will be overwritten from the biginning of
        % the buffer once the buffer is full
        % @default = Circular
        Type (1,1) string {mustBeMember(Type, ["Circular", "Fixed"]) } ...
                    = "Circular";
    end
    
    % the user cannot modify the following parameters it will be defined by
    % other internal functions, this gurantees the correct state
    properties(SetAccess=private)
  
        % The current index in the buffer @type (1,1) uint32
        CurrentIndex    (1,1) uint32 = 0;
        
        % The data stored in the buffer @type (:,:) double, each column
        % represents a data vecotr
        Data            (:,:) double = 1;
        
        % Size of the internal data @type (1,2) uint32, 
        % size = [ numVectors]
        DataSize        (1,2) uint32 = [1 1];
        
        % number of data vecotrs/sets  
        NumDataVectors  (1,1) uint32  = 1;
        
        % number of samples, used in analog sampled input setting. In this 
        % case data will be written in data(:,
        % ((index-1)*NumSamples+1):(index*NumSamples)) = data; where the
        % size of data is [ numDataPoints, numSamples]
        NumSamples      (1,1) uint32  = 1;
        
        %  alogical that is true if this.Type == "Circular", default = true
        IsCircular      (1,1) logical = true;
    end
    
        
    methods
        
        function set.Type( this, val)
        % setter for the type
        %
            this.Type = val;
            this.updateIsCircular();
        end
        
        function this = IODataBuffer( bufferSize, numSamples, type)
        % creates a new data bufer based on buffer size and 
        %   will create a buffer with data size = [bufferSize(1), bufferSize(2)*numSamples]
        %
        % @param bufferSize - number of data points x number of data
        %                     vecotrs
        % @param numSamples - number of data samples 
        % @param type - @type char, string, the type of the buffer, eithe
        %               'circular' or 'fixed'.
        
            this.initializeBuffer(bufferSize, numSamples);
            if nargin > 2 && ~isempty(type)
                this.Type = type;
            end
            
        end
        
        
        function num = numDataVectors(this)
        %returns the number of data vectors that can be stored in this
        %buffer
        %
        %  Data vectors are defined as the second dimension of the buffer
        %  data
        
            num = this.NumDataVectors;
        end
        
        
        function num = numDataPoints(this)
        %returns the number of data vectors that can be stored in this
        %buffer
        %
        %  Data vectors are defined as the second dimension of the buffer
        %  data
        
            num = this.DataSize(1);
        end
        
        
        function isF = isBufferFull(this)
        % determines whether the buffer is full
        %
        % @return isF - @type logical, true if buffer is full, false
        %               otherwise. if the buffer is circular, this will
        %               always return false
        
            % circular buffer will never be full
            if this.IsCircular
                isF = false;
            else
                isF = this.CurrentIndex >= this.NumDataVectors;
            end
        end
        
        function data = getDataVectorAt( this, index)
        % add a data vector to the buffer, array should have the size  
        %
        %
        
            if index <= this.NumDataVectors 
                numSamples = this.NumSamples;
                data = this.Data;
                data = data(:, ((index-1)*numSamples+1):(index*numSamples));
            else
                 error('IODataBuffer:indexOutOfBounds', ...
                        'Given index is out of bounds');
            end
            
        end
        
        function addDataVector(this, array )
        % add a data vector to the buffer, array should have the size  
        %
        % 
        
            if this.IsCircular || ~this.isBufferFull()
          
                % this should make this more secure
                if this.isValidArraySize(array)
                
                    numSamples = this.NumSamples;

                    index = this.CurrentIndex;
                    index = index + 1;
                    
                    % Using the difinition of data vecotrs seperate of num
                    % samples saves us a division here.
                    if index > this.NumDataVectors 
                        index = 1;
                    end

                    % decide here how to handle if the vector does not fit
                    % if not checks this might throw an error if vector is
                    % smaller then the first dimension of data if vector is
                    % larger this will cause the data array to increase

                    % this function is time critical so you want this work as
                    % much efficeint as possible some comparisons may take some
                    % time

                    data = this.Data;             
                    data(:, ((index-1)*numSamples+1):(index*numSamples)) = array;              

                    this.Data = data;

                    this.CurrentIndex = index;
                    
                else
                    error('IODataBuffer:invalidDataSize', ...
                        'Given array has invalid size');
                end
                
            end
            
        end
       
    end
    
    methods(Access=private)
        
        function isValid = isValidArraySize(this, array)
            isValid =    size(array, 1) == this.DataSize(1) ...
                      && size(array, 2) == this.NumSamples;
        end
        
        function initializeBuffer(this, bufferSize, numSamples)
        % initialize the data buffer
            
            dataSize = [ bufferSize(1), bufferSize(2) * numSamples ];
            this.NumSamples = numSamples;
            
            data = zeros(dataSize);
            if ~isempty(data)
                this.Data           = data;
                this.DataSize       = size(data);
                this.NumDataVectors = bufferSize(2);
            else
                error('IODataBuffer:initializeBuffer:bufferMustNotBeEmpty', ....
                    'Buffer data must not be empty');
            end
         end
        
        
        function updateIsCircular(this)
            if this.Type == "Circular"
                this.IsCircular = true;
            else
                this.IsCircular = false;
            end
        end
        
    end
    
end

