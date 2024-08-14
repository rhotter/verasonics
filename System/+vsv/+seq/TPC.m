classdef TPC <  vsv.seq.AbstractSeqComponent ...
                & vsv.events.EnableEventsMixin ...
                & vsv.seq.storage.HasStorageParameter
% TPC class representation of the TPC struct for ultrasound sequences
%   
%
% Version 2 | 2021-05-22 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    
     
    properties(Constant)
        
        % @type double minimum TPC voltage in Volts
        minTpcVoltage = 1.6;
        
        % @type double Maximum TPC voltage in Volts
        maxTpcVoltage = 96;
        
    end
    
    properties
        
        % @type logical indicates whether the TPC profile is in use
        inUse (1,1) double = 1; %: 1
        
        % @type logical indicates whether the TPC profile is in use
        hv (1,1) double = vsv.seq.TPC.minTpcVoltage;               %: 1.6000
        
        % @type numeric the maximum high voltage in Volts
        maxHighVoltage (1,1) double = vsv.seq.TPC.minTpcVoltage;   %: 50
        
        % @type numeric high voltage limit
        highVoltageLimit (1,1) double = vsv.seq.TPC.minTpcVoltage; %: 50
        
        % @type char the name of the TPC object, which can be used in GUIs
        % e.g.: 'Doppler'
        name (1,:) char  = '';
    end
   
    
    methods(Static, Access=protected)
        
        function hv = getDefaultHV()
        % returns the default high voltage which is minTpcVoltage    
            hv = vsv.seq.TPC.minTpcVoltage;
        end
        
    end
    
    
    %% ninitialization
    methods
        
        function this = TPC()
        % constructs a new TPC object
        %
        
            this@vsv.seq.AbstractSeqComponent( 'TPC' );
        end
        
    end
    
    %% overriden and implemented vsv.seq.AbstractSeqComponent
    methods
        
        
        
        %@override vsv.seq.storage.HasStorageParameter
        function ID = getStorageID(this)
            
            tpcName = convertStringsToChars( this.name );
            if isempty(tpcName)
                ID = [ 'TPC_' num2str( this.getListIndex ) ];
            else
                if isempty( regexp( tpcName, 'TPC', 'ONCE' ) )
                    ID = [  'TPC_' tpcName];
                else
                    ID = tpcName;
                end
            end
            
        end
        
        %@override
        function [props, values] = listImportProperties(this)
        
            props = {   'inUse'; ...
                        'hv'; ...
                        'maxHighVoltage'; ...
                        'highVoltageLimit'; ...
                        'name'};
            if nargout > 1
                values = {  this.inUse; ...
                            this.hv; ...
                            this.maxHighVoltage; ...
                            this.highVoltageLimit; ...
                            this.name };
            end
        end
        
        %@override
        function new = newInstance(~)
            new = vsv.seq.TPC;
        end
        
        
       
    end
    
    %% setter methods
    methods
        
        function set.highVoltageLimit(this, value)
            
            if value >= this.minTpcVoltage
                this.highVoltageLimit = value;
                this.updateHv();
            else
                error('TPC:set:highVoltageLimit', 'hv value is out of range, must be larger or equal then minTpcVoltage');
            end
        end   
        
        function set.maxHighVoltage(this, value)
            
            if value >= this.minTpcVoltage
                this.maxHighVoltage = value;
                this.updateHv();
            else
                error('TPC:set:maxHighVoltage', 'hv value is out of range, must be larger or equal then minTpcVoltage');
            end
            
        end
        
        function set.hv(this, value)
        % @type logical indicates whether the TPC profile is in use
            
            if this.isHvValid(value)
                this.hv = value;
                this.notifyPropertyChanged('hv');
            else
                error('TPC:set:hv', 'hv value is out of range');
            end
        end
        
        
        %@override vsv.seq.AbstractSeqComponent make sure that we import hv at last
        function importStruct(this, str, fields)
        
            
            if nargin < 3
                fields = fieldnames(str);
            end
        
            ind = strcmp(fields, 'hv');
            if any(ind)
                hV = str.hv;
                fields(ind) = [];
                this.importStruct@vsv.seq.AbstractSeqComponent( str, fields);
                this.hv = hV;
            else
                this.importStruct@vsv.seq.AbstractSeqComponent( str, fields);
            end
                       
        end
        
                
        function isv = isHvValid(this, hv)
        % high voltage slider is valid if its larger or equal the min
        % voltage and smaller then the smaller of maxHighVoltage or
        % highVoltageLimit
        %
        
            limits = this.hvLimits();
            isv =    hv >= limits(1) ...
                  && hv <= limits(2);
        end
        
        function lim = hvLimits(this)
            
            minHV = this.minTpcVoltage;
            maxHV = min( this.maxHighVoltage, this.highVoltageLimit );
            lim = [minHV, maxHV];
            
        end
                
    end
    
    methods(Access=protected)
        
        %@override vsv.seq.storage.HasStorageParameter
        function [ params, vals ] = getStorageParameter(this)
            params = {'hv'};
            vals   = {this.hv};
        end
        
    end
    
    methods(Access=private)
        
        function updateHv(this)
            if this.hv > this.maxHighVoltage || this.hv > this.highVoltageLimit
                this.hv = min( [this.maxHighVoltage, this.highVoltageLimit] );             
            end
        end
        
    end
   
    
end

