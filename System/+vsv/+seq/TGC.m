classdef (Sealed) TGC  <   vsv.seq.AbstractSeqComponent ...
                         & vsv.events.EnableEventsMixin ...
                         & vsv.seq.storage.HasStorageParameter
%TGC is a sequence component containing informaton about time gain
%compensation
%
%   TGC is mainly used to compute the time gain compensation wave form. The
%   waveform is calculated based on the parameters that are specified in
%   this class and will only be calculated when requested (i.e., when the 
%   user queries TGC.Waveform)
%
%   The function that computes the wave form is computTGCWaveform which
%   will be automatically called when the user asks for the Wavform. 
%   The computeTGCWaveform function first computes rangeMax in number of
%   800nsec intervals, which is the TGC generator sample period.  
%   It then divides this distance by 7 to find the depth of the control 
%   points (1st point is at 0 range).  These points, and a duplicate of 
%   the eighth point at max range (511) are used to compute a piecewise 
%   cubic Hermite waveform at the TGC sample period.
%
%
% Version 1.0 | 2020-03-31 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    
    %% Control events
    events
        
        % indicates that the control points have cahnged. It will fire an
        % vsv.ctrl.evt.TGCCtrlEvt providing the index of the 
        CntrlPtsChanged;
        
        % indicates that the wave form was recalculated
        WaveFormReset;
        
        % indicates that the gain has changed
        GainChanged;
        
    end
      
    %% main public properties
    properties
        
        % @type numeric the range max value. The range max value gives the
        % depth range max for computing the TGC value
        rangeMax (1,1) double = 128;
        
        % @type numeric arrays giving the control points. A control point
        % can be understand as a single point on the the TGC wave form. The 
        % TGC wave form is calculated b
        CntrlPts (1,:) double = vsv.seq.TGC.defaultCtrlPts();  
        
        % the center frequency is needed to calculate the TGC value and
        % need to be set
        CenterFreq (1,1) double {mustBeNonnegative} = 7e6;
        
        % @type char optional name of the TGC object to display in GUIs
        % 
        NameID (1,1) string ;
       
    end
    
    
    %% Properties Dependent on stucture data
    properties(SetAccess=private)
        
        % @type numeric array the waveform, the waveform is getting
        % calculated based on teh CntrlPts
        Waveform (1,:) double = [];
        
    end
    
    properties(SetAccess=private)
        
        % Gain factor for the contrl points @type numeric
        Gain = 1;
        
        % The gain for the TGC @type numeric 2 element vector
        TGCGainLimits = [0.1 4];
        
    end
    
    
    %% Private properties
    properties(Access=private)
        
        % to save the control points before gain was applied
        OriginCntrlPoints
        
        % maximum samples allowed 
        MaxSamples = 511;
        
        % maximum TGC val allowed for a single cntrl point
        MaxTGCVal  = 1023;
        
        % minimum TGC Val allowed for a single cntrl point
        MinTGCVal  = 0;
        
        % indicate whether wave form should be recomputed
        DoRecomputeWF = true;
        
        % @type logical is used to silence the updating of the origin
        DoUpdateOrigin (1,1) logical = true;
                
        % @type logical is used to silence the updating of control point
        % changed 
        DoNotifyCntrlPtsChanged (1,1) logical = true;
        
    end
    
    
       
    %% static methods 
    methods(Static)
        
        function default = defaultCtrlPts()
        % returns a default set of control points
        %   default = [300,450,575,675,750,800,850,900]
        %
        %@return default - @type numeric, 8 element array 
        
            default = [300,450,575,675,750,800,850,900];
        end
        
    end
    
    %% Abstract methods implementation
    methods
        
        %@override vsv.seq.storage.HasStorageParameter
        function params = listStorageParameter(~)
            params = {'CntrlPts', 'Gain' };
        end
        
        %@override vsv.seq.storage.HasStorageParameter
        function ID = getStorageID(this)
            ID = convertStringsToChars( strcat( 'TGC_', this.NameID ) );
        end
        
        %@override
        function new = newInstance(this)
        
            new = vsv.seq.TGC( this.CenterFreq, ...
                                'rangeMax', this.rangeMax, ...
                                'CntrlPts', this.CntrlPts);
        end
        
        %@override
        function [ props, values] = listImportProperties(this)
        
            props = {   'rangeMax'; ...
                        'CntrlPts'; ...
                        'CenterFreq'; ...
                        'NameID'};
                    
            if nargout > 1
                values = {  this.rangeMax; ...
                            this.CntrlPts; ...
                            this.CenterFreq; ...
                            this.NameID };
            end
        end
        
        %@overwrite to provide extra support for Waveform which is not part
        %of the listImport Properties
        function str = exportStruct(this)
        
            [ props, vals ] = this.listImportProperties(); 
            
            props{end+1} = 'Waveform';
            vals{end+1}  = this.Waveform;
            
            str = cell2struct( vals(:), props(:), 1);
        end
       
    end
    
    %% public general methods
    methods
        
        function this = TGC(centerFreq, varargin)
        % constructor that will create the TGC controll
        %
        %   this = TGC(TGC); will init the TGC control based on TGC
        %   which should be a struct containing rangeMax and CntrlPts
        %
        %   this = TGC('PropertyName', PropertyValue); creates a TGC
        %   control similar to a struct with name value pairs
        %
        % @param centerFreq - @type numeric the center freq 
            
            % this is a seq component with TGC as the component name
            this@vsv.seq.AbstractSeqComponent('TGC');
            
            if nargin >= 1 
                if ~this.isCenterFreqValid(centerFreq)
                    error('TGC:frequencyinvalid', 'invalid frequency (Freq. must be larger than 0)');
                end
                this.CenterFreq = centerFreq;                
            end
            if nargin >= 2
                if isstruct(varargin{1})
                    this.importStruct( varargin{1} );
                else
                    
                    params = varargin(1:2:end);
                    vals   = varargin(2:2:end);
                    str = cell2struct(vals(:), params(:), 1);
                    
                    this.importStruct(str);
                end
            end            
            this.CompName = 'TGC';
            
        end
        
        
        function isv = isCenterFreqValid(~, centerFreq)
        % checks whether the center freq is valid. 
        %
        % @param centerFreq - @type numeric && > 0 
        % @return isv - @type logical true if cenerFreq > 0 and numeric,
        %               false otherwise
            isv = centerFreq > 0;
        end
        
       
        function isv = isSetableAttributes(~, attribs)
            
           isv = numel(attribs) == 1 && strcmp( attribs{1}, 'Waveform' );
        end
        
        
    end
    
        
    
    %% Main methods to control TGC
    methods
        
        function setGain(this, gain )
        % Sets the TGC gain 
        %   The gain will increase all TGC Ctrl point values 
        %
        % @param gain - @type numeric 
        
        
            if this.isValidGain(gain)
                
                this.Gain = gain;
                this.notifyEnabled('GainChanged');
                
                % this are the new contrl points
                cntrPnts = this.OriginCntrlPoints .* gain;
                
                this.DoUpdateOrigin = false;
                try
                    this.CntrlPts = cntrPnts;
                catch err
                    this.DoUpdateOrigin = true;
                    rethrow(err);
                end
                this.DoUpdateOrigin = true;
                
            
            end
            
        end
        
        function isv = isValidGain(this, gain)
        % Gain is valid if its within its limits 
        %
        %   see TGCGainLimits, this.getTGCLimits
        % @param gain - @type numeric, the gain to check
        
            limits = this.TGCGainLimits;
            isv    = gain >= limits(1) < limits(2);
            
        end
        

        function limits = getTGCLimits(this)
        % returns the limits of the TGC min and max value that can be set
        % 
        % @return limits - @type numeric, two element array with min and
        %                   max value
            limits = [this.MinTGCVal this.MaxTGCVal];
        end
       
        function ctrl = validateCtrlPts(this, val)    
        % Function that validates the control points
        %
        % this will return a value in any case. If val is outside of the
        % range max and min values the returned value will be set to the
        % closest boundary value.
        %
        % @param val - @typeimp numeric, the value to be set
        % @return ctrl - @type numeric, either val, or this.MaxTGCVal if val
        %                > this.MaxTGCVal, or this.MinTGCVal if val <
        %                this.MinTGCVal 
        
            ctrl = min( this.MaxTGCVal, max( this.MinTGCVal, val ) );
        end
        
        function val = getTGCVal(this, ind)
        % returns the TGCVal at a given index ind
        % 
        % @param ind - @type numeric the index pointing to the control
        %              point in CntrPnts
        % @return val - @type numeric , this.CntrPnts(ind)
            if ind < 1 || ind > this.numControlPoints()
                error('TGC:setTGCVal:indexOutOfBounds', 'Index is out of bounds');
            else
                val = this.CntrlPts( ind );
            end
        end
        
        function setTGCVal(this, ind, val)
        % will set a single TGCVal at a given index
        %
        % This is the preffered way of setting single TGC values, because
        % this will be more efficient then setting the entire cntrl points
        %
        % @Error TGC:setTGCVal:indexOutOfBounds - if given index ind is out
        %        of bounds
        %
        % @param ind - @type numeric the index of the 
        % @param val - @type numeric the control point value to be set at
        %              given index
        
            if ind < 1 || ind > this.numControlPoints()
                error('TGC:setTGCVal:indexOutOfBounds', ...
                    'Index is out of bounds');
            else
                val = this.validateCtrlPts(val);
                % do only change if value is different from before
                if val ~= this.CntrlPts(ind)
                    
                    this.DoUpdateOrigin = false;
                    this.DoNotifyCntrlPtsChanged = false;
                    
                    try
                        this.CntrlPts( ind ) = val;
                    catch err
                        this.DoNotifyCntrlPtsChanged = true;
                        this.DoUpdateOrigin = true;
                        rethrow(err);
                    end
                    
                    this.DoNotifyCntrlPtsChanged = true;
                    this.DoUpdateOrigin = true;
                    this.updateOrigin(ind);
                    
                    this.notifyCntrlPtsChanged(ind);
                end
            end
            
        end
        
        
        function num = numControlPoints(this)
        % returns the number of control points
        %
        % @return num - @type numeric the number of control points which is
        %               usually 8
        
            num = numel(this.CntrlPts);
        end
        
        function importTGC(this, TGC)
        % Allows to import a struct TGC with CntrlPts and rangeMax fields
        % This is to support old struct data 
        %
           this.importStruct(TGC);
        end
        
        function samples = getNumberOfSamples(this)
        % returns the number of samples to represent the wave form
        % 
            samples = min(this.MaxSamples-0.1, ...
                          round(2*double(this.rangeMax)/800e-09/this.CenterFreq));
        end
        
        function computeTGCWaveform(this )
        % this is basic computeTGCWaveform function 
        %
        % This function does not need to be called, it will be called
        % automatically once the cntrl points are changing and someone
        % requests 
        
            if this.doRecomputeWF()
                nc = this.numControlPoints();

                samples = this.getNumberOfSamples();
                
                % Note the limit of 510.9 ensures that the eighth control point position
                % will always be less than the ninth one at 511, since the interpolation
                % routine requires that they are monotonic.

                % Construct P and C arrays to specify known points on waveform.
                %P =    zeros(1,nc+1);
                P = [ 0:(samples/(nc-1)):samples this.MaxSamples];

                cntrlPts = this.CntrlPts;
                C = double( [ cntrlPts(:)', cntrlPts(end)]);
                X = 0:this.MaxSamples;

                % Compute the interpolated points.
                this.Waveform =  round(interp1(P,C,X,'pchip'));
                this.DoRecomputeWF = false;
                
            end
        end
        
    end
    
    %% getter and setter
    methods
        
        function set.CntrlPts(this, value)
        % setter for teh contrl pts 
        % this will validate the value and set it to the underlying data 
        %
        % @param val - @type numeric the value supposed to set
        
            value = this.validateCtrlPts(value);
            this.CntrlPts = value;
            this.resetWaveform();
            this.notifyCntrlPtsChanged(-1);
            this.updateOrigin();
        end
        
        function set.rangeMax(this, val)
        % setter for range max
            this.rangeMax = val;
            this.resetWaveform();
        end
            
        function set.CenterFreq(this, val)
        % setter for frequency. Frequency must be in Hz
        %
        % 
        
            if val < 100e3
                warning('TGC:setCenterFreq:fcinHz', ... 
                    [ 'Frequency must must be in Hz, Given value ' ...
                       num2str(val) ' is most likely in MHz']);
            end
            this.CenterFreq = val;
            this.resetWaveform();
        end
        
        function val = get.Waveform(this)
        % getter for the waveform
        
            if this.doRecomputeWF()
                this.computeTGCWaveform();
            end
            val = this.Waveform ;
        end
                
    end
    
    %% supporting protected methods
    methods(Access=protected)
        
        %@override vsv.seq.storage.HasStorageParameter
        function [paramName, values] = getStorageParameter(this)
            values    = {this.CntrlPts, this.Gain};
            paramName = { 'CntrlPts', 'Gain'}; 
        end
        
        %@override vsv.seq.storage.HasStorageParameter
        function [ success, msg] = importStorageParameter_p(this, param)
        
            % find cntrl pts
            paramNames = param.Parameter;
            values     = param.Values;
            
            try
                
                % first apply the gain then the cntrl points
                indGain = strcmp(paramNames, 'Gain');
                if any(indGain)
                    gain = values(indGain);
                    this.setGain( gain{1} );
                end

                % first apply the gain then the cntrl points
                indCntrlPts = strcmp(paramNames, 'CntrlPts');
                if any(indCntrlPts)
                    cntrlPts = values(indCntrlPts);
                    this.CntrlPts = cntrlPts{1};
                end
                success = true;
                msg = '';
            catch err
                success = false;
                msg = {err.message}; 
            end
            
        end
        
        
        function resetWaveform(this)
        % this will cause the waveform to be recalculated the next time
        % it is requested
        %
        % This will also notify event WaveFormReset to indicate that the
        % waveform needs to be recalculated
        %
            this.DoRecomputeWF = true;
            this.notifyPropertyChanged('Waveform');
        end
        
        function doRecompute = doRecomputeWF(this)
        % returns true if the waveform needs to be recalculated
        %
        % @return doRecompute - @type logical, true if the waveform needs
        %                       to be recalculated
        
            doRecompute = this.DoRecomputeWF ;
        end
    end
    
    %% getter and setter functions to set properties to underlying struct 
    % data
    % these methods are setters for fast value assignment, and no value
    % evaluation
    methods(Access=private)
        
        function updateOrigin(this, index)
        % The origin contrl points are saved to restore them when the gain
        % is changing
        %
        % If the control points are changeing this will store that state.
        % This allows to restore the previous TGC contrl points when
        % changeing the Gain 
        
            if this.DoUpdateOrigin
                if nargin > 1
                    this.OriginCntrlPoints(index) = this.CntrlPts(index)/this.Gain;
                else

                    this.OriginCntrlPoints = this.CntrlPts./this.Gain;
                end
            end
            
        end
        
                
        function notifyCntrlPtsChanged(this, index)
        % 
        %   
        
        
            persistent events n all;
            
            if this.DoNotifyCntrlPtsChanged
            % cashing events rather then recreating it brings tremendous
            % speed improvements!
            
                

                if index < 1
                    if isempty(all)
                        all = vsv.seq.evt.TGCCtrlPtEvt(index);
                    end
                    this.notify( 'CntrlPtsChanged', all );
                else
                    if isempty(events)
                        events{index} = vsv.seq.evt.TGCCtrlPtEvt(index);
                        n = index;
                    elseif index > n
                        events{index} = vsv.seq.evt.TGCCtrlPtEvt(index);
                        n = index;
                    elseif isempty(events{index})
                        events{index} = vsv.seq.evt.TGCCtrlPtEvt(index);
                    end
                    this.notify( 'CntrlPtsChanged', events{index} );
                end
            end
            
        end
    end
    
   
end
    

