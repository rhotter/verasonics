classdef ProcessDoppler < vsv.seq.AbstractSeqComponent...
                        & vsv.seq.storage.HasStorageParameter
% PROCESSDOPPLER
%
% Copyright (C) 2001-2022, Verasonics, Inc.  All worldwide rights and
% remedies under all intellectual property laws and industrial property
% laws are reserved.

    properties(SetObservable)
        
        % This is here for future support of doppler PRF
        % @type double PRF, range 0 – 99
        % prf (1,1) double { mustBeLessThanOrEqual(prf, 99),...
        %                    mustBeGreaterThanOrEqual(prf, 0)  } = 0;
    
        % @type double Power Threshold.
        pwrThreshold (1,1) double { mustBeLessThanOrEqual(pwrThreshold, 1.0),...
                                    mustBeGreaterThanOrEqual(pwrThreshold, 0.0)  } = 0.5;
        
        info (1, :) char = '';

    end
    
    methods
        function params = replaceParameters( this, params )
            ind = find( strcmp(params, 'pwrThreshold') );

            if ~isempty(ind)
                params{ind+1} = this.pwrThreshold;
            else
                params = [ params {'pwrThreshold', this.pwrThreshold} ];
            end
            
            % This is here for future support of doppler PRF
            % params = this.findAndReplace(params, 'prf');
            params = this.findAndReplace(params, 'pwrThreshold');     
        end
        
        function [list, values] = listImportProperties( this )
            list = this.listSupportedProperties();

            if nargout > 1
                values = { this.pwrThreshold };
                        %    this.prf };
            end
        end
        
        function list = listSupportedProperties( ~ )
            list =  { 'pwrThreshold' };
                    %   'prf' };
        end
        
        function setDefault( this )
            % This is here for future support of doppler PRF
            % this.prf          = 1; 
            this.pwrThreshold = 0.5; % Not sure what the defualt should be
            this.info         = '';        
        end
        
        function this = ProcessDoppler( varargin )
            this@vsv.seq.AbstractSeqComponent( 'Process' );

            if nargin == 1
                if isstruct(varargin{1})
                    this.importStruct(varargin{1});
                elseif iscell( varargin{1} )
                    params = varargin{1};
                    this.importParamValue(params(1:2:end), params(2:2:end));
                else
                    error('ProcessDoppler:invalidArgument', ...
                          'If one input provided it must be a struct or cell array');
                end
            elseif nargin > 1
                this.importParamValue(varargin(1:2:end), varargin(2:2:end) );
            end
        end
        
        function limits = getLimitForProperty( ~, property )
            switch property
                % case 'prf' % This is here for future support of doppler PRF
                %     limits = [0 99];
                case 'pwrThreshold'
                    limits = [0.0 1.0];
                otherwise
                    error('ProcessDoppler:invalidProperty', 'given property is unkown');
            end
        end
    
        function importParamValue( this, param, values )
            props = this.listSupportedProperties();
            
            nParams = numel(param);
            
            for i = 1:nParams
                if ismember( param{i}, props)
                    this.(param{i}) = values{i};
                end
            end 
        end
    end

    methods
        function li = getListIndex( this )
        % returns the list index of the sequence component in the lower
        % level sequence list
            li = this.ListIndex;
        end
        
        function isp = isPartOfList( ~ )
        % returns true if this sequence component is connected to a list. 
        % isp 
            isp = true;
        end
        
        function importStruct( this, str, fields )
        % returns true if this sequence component is connected to a list. 
        % isp 
        
            if nargin < 3
                fields = fieldnames(str);
            end
            
            nfields = numel(fields);
            for i = 1:nfields
                this.(fields{i}) = str.(fields{i});
            end
        
        end
        
        function ns = newInstance( ~ )
        % create a new instance
        %
            ns = vsv.seq.ProcessDoppler;
        end
    end
    
    %% getter and setter
    methods
        % This is here for future support of doppler PRF
        % function set.prf( this, value )
        %     this.prf = value;
        %     this.notifyPropertyChanged('prf');
        % end
        
        function set.pwrThreshold( this, value )
            this.pwrThreshold = value;
            this.notifyPropertyChanged('pwrThreshold');
        end
        
        function ID = getStorageID( this )
            %@override vsv.seq.storage.HasStorageParameter
            infoStr = convertStringsToChars(this.info);
            
            if isempty(infoStr)
                ID = [ 'DopplerProcess_' num2str( this.getListIndex ) ];
            else
                ID = [ 'DopplerProcess_' infoStr ];    
            end
        end
    end
    
    methods(Access=protected)
        function [params, values] = getStorageParameter( this )
            [params, values] = this.listImportProperties();
        end
        
        function params = findAndReplace( this, params, name )
            ind = find( strcmp(params, name) );

            if ~isempty(ind)
                params{ind+1} = this.(name);
            else
                params = [ params { name, this.(name) } ];
            end
        end
    end
end
