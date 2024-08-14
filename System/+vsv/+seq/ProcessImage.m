classdef ProcessImage  <  vsv.seq.AbstractSeqComponent...
                         & vsv.seq.storage.HasStorageParameter
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    

    properties(SetObservable)
        
        % @type char array 'intensity' 'signedColor', 'unsignedColor'
        srcData char {mustBeMember(srcData, ...
                        {'intensity' 'signedColor', 'unsignedColor'}) } ...
                        = 'intensity'; 
        
        % @type char array'none'(default),'low','medium','high'
        grainRemoval char {mustBeMember(grainRemoval, ...
                        { 'none', 'low','medium','high'}) } ...
                        = 'none';
        
        % @type char 'none'(default), 'simple', 'dynamic'
        persistMethod char {mustBeMember(persistMethod, ...
                        {'none' 'simple', 'dynamic'}) } ...
                        = 'none'; 
        
        % @type double Amount of persistence, range 0 – 99
        persistLevel (1,1) double { mustBeLessThanOrEqual(persistLevel, 99),...
                              mustBeGreaterThanOrEqual(persistLevel, 0)  } = 0;
        
        % @type char  Interpolation method ('4pt' interp(default)).
        interpMethod char {mustBeMember(interpMethod, {'4pt', '8pt'}) } ...
                     = '4pt'; 
        
        % @type char 'none'(default),'reduceSpeckle1','reduceSpeckle2'
        processMethod char {mustBeMember(processMethod, ...
                            {'none','reduceSpeckle1','reduceSpeckle2' }) } ...
                            = 'none'; 
        
        % @type char 'none'(default),'runAverage2','runAverage3'
        averageMethod char {mustBeMember(averageMethod, ...
                            {'none','runAverage2','runAverage3' }) } ...
                            = 'none'; 
    
        % @type double  Processing gain.
        pgain (1,1) double { mustBeLessThanOrEqual(pgain, 100),...
                        mustBeGreaterThanOrEqual(pgain, 0)  } = 1;
        
        % @type double  Low level reject (0-100)
        reject (1,1) double { mustBeLessThanOrEqual(reject, 100),...
                        mustBeGreaterThanOrEqual(reject, 0)  } = 0;
        
        % @type char 'none'(default),'runAverage2','runAverage3'
        compressMethod char {mustBeMember(compressMethod, ...
                            {'power','log'}) } ...
                            = 'power'; 
        
        % @type double  Low level reject (0-100)
        compressFactor (1,1) double { ...
                        mustBeLessThanOrEqual(compressFactor, 100),...
                        mustBeGreaterThanOrEqual(compressFactor, 0)  } = 40;
                
        % @type char  'full'(default),'lowerHalf','upperHalf'.
        mappingMethod char {mustBeMember(mappingMethod, ...
                            {'full','lowerHalf', 'upperHalf'}) } ...
                            = 'full'; 
        % @type double  0-255. Compare against pixel value to determine whether to overwrite.
        threshold (1,1) double { ...
                        mustBeLessThanOrEqual(threshold, 255),...
                        mustBeGreaterThanOrEqual(threshold, 0)  } = 255;
                    
        info (1, :) char = '';
    end
    
    methods
        
        
        function params = replaceParameters(this, params)
        
            ind = find( strcmp( params, 'srcData' ) );
            if ~isempty(ind)
                params{ind+1} = this.srcData;
            else
                params = [ params { 'srcData', this.srcData } ];
            end
            
            params = this.findAndReplace(  params, 'srcData');
            params = this.findAndReplace(  params, 'grainRemoval');
            params = this.findAndReplace(  params, 'persistMethod');
            params = this.findAndReplace(  params, 'persistLevel');
            params = this.findAndReplace(  params, 'interpMethod');
            params = this.findAndReplace(  params, 'processMethod');
            params = this.findAndReplace(  params, 'averageMethod');
            params = this.findAndReplace(  params, 'pgain');
            params = this.findAndReplace(  params, 'reject');
            params = this.findAndReplace(  params, 'compressMethod');
            params = this.findAndReplace(  params, 'compressFactor');
            params = this.findAndReplace(  params, 'mappingMethod');
            params = this.findAndReplace(  params, 'threshold');     
        end
        
        function [ list, values] = listImportProperties(this)
            list = this.listSupportedProperties();
            if nargout > 1
                values = {  this.srcData; ...
                            this.grainRemoval; ...
                            this.persistMethod; ...
                            this.persistLevel; ...
                            this.interpMethod; ...
                            this.processMethod; ...
                            this.averageMethod; ...
                            this.pgain; ...
                            this.reject; ...
                            this.compressMethod; ...
                            this.compressFactor; ...
                            this.mappingMethod; ...
                            this.threshold};
            end
        end
        
        function list = listSupportedProperties(~)
            list =  {   'srcData'; ...
                        'grainRemoval'; ...
                        'persistMethod'; ...
                        'persistLevel'; ...
                        'interpMethod'; ...
                        'processMethod'; ...
                        'averageMethod'; ...
                        'pgain'; ...
                        'reject'; ...
                        'compressMethod'; ...
                        'compressFactor'; ...
                        'mappingMethod'; ...
                        'threshold'};
        end
        
        function setDefault(this)
            
            % runAcq default values
            % Image->imgbufnum = 1;
            %             Image->framenum = -1;
            %             Image->sectionnum = 1;
            %             Image->srcData = Im_srcDataIntensity2D;
            %             Image->pdatanum = 1;
            %             Image->norm = 0;
            %             Image->pgain = 1;
            %             Image->reject = 0.0;
            %             Image->grainRemoval = Im_grainNone;
            %             Image->persistMethod = Im_persistNone;
            %             for (k=0; k<4; k++) Image->PersistLevel[k] = 0;
            %             Image->interpMethod = Im_interp4pt;
            %             Image->processMethod = Im_procMethodNone;
            %             Image->averageMethod = Im_averMethodNone;
            %             Image->compressMethod = Im_compMethodPower;
            %             Image->compressFactor = 40; // for sqrt compression with 'power'
            %             Image->mappingMethod = Im_mappingMethodFull;
            %             Image->threshold = 255;
            %             Image->displayWindow = 1;
            %             Image->display = 1;
            %             Image->extDisplay = 0;
            %             Image->InterpData = NULL;
            %             Image->interpLUTSize = 0; // for polar coord data
            %             Image->CalculateLUT = NULL;
          
            this.srcData        = 'intensity'; 
            this.grainRemoval   = 'none';
            this.persistMethod  = 'none'; 
            this.persistLevel   = 0;
            this.interpMethod   = '4pt'; 
            this.processMethod  = 'none'; 
            this.averageMethod  = 'none'; 
            this.pgain          = 1;
            this.reject         = 0;
            this.compressMethod = 'power'; 
            this.compressFactor = 40;
            this.mappingMethod  = 'full'; 
            this.threshold      = 255;
            this.info           = '';        
        end
        
        function this = ProcessImage(varargin)
            
            this@vsv.seq.AbstractSeqComponent( 'Process' );
            if nargin == 1
                if isstruct(varargin{1})
                    this.importStruct(varargin{1});
                elseif iscell( varargin{1} )
                    
                    params = varargin{1};
                    this.importParamValue(params(1:2:end), params(2:2:end));
                else
                    error('ProcessImage:invalidArgument', ...
                        'If one input provided it must be a struct or cell array');
                end
            elseif nargin > 1
                this.importParamValue(varargin(1:2:end), varargin(2:2:end) );
            end
        end
        
        
        function limits = getLimitForProperty(~, property )
            
            switch property
                case  'srcData'
                    limits = {'intensity' 'signedColor', 'unsignedColor'};    
                case 'grainRemoval' 
                    limits = { 'none', 'low','medium','high'};
                case 'persistMethod' 
                    limits = {'none' 'simple', 'dynamic'};
                case 'persistLevel' 
                    limits = [0 99];
                case 'interpMethod' 
                    limits = {'4pt'}; 
                case 'processMethod' 
                    limits = {'none','reduceSpeckle1','reduceSpeckle2' };
                case 'averageMethod' 
                    limits = {'none','runAverage2','runAverage3' };
                case { 'pgain', 'reject', 'compressFactor' }
                    limits = [0 100];
                case 'compressMethod'
                    limits = {'power','log'}; 
                case 'mappingMethod' 
                    limits = {'full','lowerHalf', 'upperHalf'};
                case 'threshold'
                    limits = [0 255];
                otherwise
                    error('ProcessImage:invalidProperty', 'given property is unkown');
            end
     
        end
    
        function importParamValue(this, param, values)
            
            props = this.listSupportedProperties();
            
            nParams = numel(param);
            
            for i = 1:nParams
                if ismember( param{i}, props)
                    this.(param{i}) = values{i};
                end
            end 
            
        end
                
    end
    
    %% 
    methods
        
        
        
        function li = getListIndex(this)
        % returns the list index of the sequence component in the lower
        % level sequence list
        %
            li = this.ListIndex;
        end
        
        function isp = isPartOfList(~)
        % returns true if this sequence component is connected to a list. 
        % isp 
            isp = true;
        end
        
        function importStruct(this, str, fields)
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
        
        function ns = newInstance(~)
        % create a new instance
        %
            ns = vsv.seq.ProcessImage;
        end
    end
    
    %% getter and setter
    methods
        % @type char array 'intensity' 'signedColor', 'unsignedColor'
        function  set.srcData(this, value)
            this.srcData = value;
            this.notifyPropertyChanged('srcData');
        end
        
        % @type char array'none'(default),'low','medium','high'
        function  set.grainRemoval (this, value)
            this.grainRemoval = value;
            this.notifyPropertyChanged('grainRemoval');
        end
        
        % @type char 'none'(default), 'simple', 'dynamic'
        function  set.persistMethod(this, value)
            this.persistMethod = value;
            this.notifyPropertyChanged('persistMethod');
        end
         
        % @type double Amount of persistence, range 0 – 99
        function  set.persistLevel(this, value)
            this.persistLevel = value;
            this.notifyPropertyChanged('persistLevel');
        end
         
        % @type char  Interpolation method ('4pt' interp(default)).
        function  set.interpMethod(this, value)
            this.interpMethod = value;
            this.notifyPropertyChanged('interpMethod');
        end
         
        % @type char 'none'(default),'reduceSpeckle1','reduceSpeckle2'
        function  set.processMethod(this, value)
            this.processMethod = value;
            this.notifyPropertyChanged('processMethod');
        end
         
        
        % @type char 'none'(default),'runAverage2','runAverage3'
        function  set.averageMethod(this, value)
            this.averageMethod = value;
            this.notifyPropertyChanged('averageMethod');
        end
         
    
        % @type double  Processing gain.
        function  set.pgain(this, value)
            this.pgain = value;
            this.notifyPropertyChanged('pgain');
        end
        
        % @type double  Low level reject (0-100)
        function  set.reject(this, value)
            this.reject = value;
            this.notifyPropertyChanged('reject');
        end
        
        % @type char 'none'(default),'runAverage2','runAverage3'
        function  set.compressMethod(this, value)
            this.compressMethod = value;
            this.notifyPropertyChanged('compressMethod');
        end
        
        
        % @type double  Low level reject (0-100)
        function  set.compressFactor(this, value)
            this.compressFactor = value;
            this.notifyPropertyChanged('compressFactor');
        end
        
        % @type char  'full'(default),'lowerHalf','upperHalf'.
        function  set.mappingMethod(this, value)
            this.mappingMethod = value;
            this.notifyPropertyChanged('mappingMethod');
        end
        
        % @type double  0-255. Compare against pixel value to determine whether to overwrite.
        function  set.threshold(this, value)
            this.threshold = value;
            this.notifyPropertyChanged('threshold');
        end
        
        function ID = getStorageID(this)
            
            %@override vsv.seq.storage.HasStorageParameter
            infoStr = convertStringsToChars(this.info);
            
            if isempty(infoStr)
                ID = [ 'ImageProcess_' num2str( this.getListIndex ) ];    
            else
                ID = [ 'ImageProcess_' infoStr ];    
            end
                
        end
        
    end
    
    
    methods(Access=protected)
        
        
        function [ params, values] = getStorageParameter(this)
            [ params, values] = this.listImportProperties();
        end
        
        
        function params = findAndReplace( this, params, name)
            
            ind = find( strcmp( params, name ) );
            if ~isempty(ind)
                params{ind+1} = this.(name);
            else
                params = [ params { name, this.(name) } ];
            end
        end
        
    end
end

