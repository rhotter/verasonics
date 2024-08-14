classdef VsCaseControl < vsv.seq.uicontrol.UserUiControl ...
                        & vsv.seq.uicontrol.mixin.HasTitle ...
                        & vsv.ctrl.mixin.CreatePropertyControlSetting ...
                        & vsv.events.EnableEventsMixin ...
                        & vsv.seq.storage.HasStorageParameter ...
%VsCaseControl A control that allows to select cases among a possible selection
%of cases
%
%   This is the base class for VsButtonGroups and for VsDropDown controls
%
%    
%
% Version 1.0 | 2020-01-28 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    

    properties(SetAccess=private)
        
        % @type char the case that is currently selected
        SelectedCase;
        
        % @type numeric the selected index
        SelectedIndex;
        
        % @type cell the possible cases 
        PossibleCases
        
        % @type the number of cases
        NumCases;
        
    end
    
    events
        
        % this indicates if the selection has changed
        SelectionChanged;
        
    end
    
    methods
        
        function strUI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
        
            strUI.Control  = {this.LocationCode, ...
                              'Style',          this.Style,     ...
                              'Title',          this.Title,     ...
                              'NumButtons',     this.NumCases,  ...
                              'Labels',         this.PossibleCases };
            strUI.Callback = this.Callback;
            if ~vsv.util.isNULL( this.Statement )
                strUI.Statement = this.Statement;
            end
        end
        
        
        % @override vsv.seq.storage.HasStorageParameter
        function ID = getStorageID(this)
            ID = this.Title;
        end
        
        function this = VsCaseControl(varargin)
        % Constructor for the VsButtonGroupControl class
        %
        %   the button group control must be initialized with the possible
        %   cases that can be selected
        %
        % Errors
        %   @error VSButtonGroupContro:setPossibleCases:invalidArgument if
        %          posc is not a cellstring
        %
        % @param possibleCases - possible cases @type cellstring optional
        % @param UI - old UI definition
        %
            this.initVsCaseControl(varargin{:});
        end
        
        function isvalid = isValidStyle(~, style)
            isvalid = ~isempty(style) && strcmp( style, 'VsCaseControl' );
        end
        
        function index = getIndexOfCase(this, caseName)
        % returns the index of the given case 
        %
        % @param caseName - @type char, the name of the case of which to
        %                   return the index
        %
        % @preturn index - @type the index of the caseName, or []
        
            index = [];
            if ischar(caseName) || ( isstring(caseName) && numel(caseName) == 1)
                index = this.getIndexOfCasePrivate(caseName);
            end
            
        end
        
        function isvalid = isValidCase(this, caseName)
        % Checks whether a case is valid, i.e., the selection is part of
        % the PossibleCases
        %
        % @param caseName - @type char the case to check for 
        % @return isvalid - @type logical true, if caseName is valid, false
        %                 otherwise
            
            isvalid = false;
            if ischar(caseName) || ( isstring(caseName) && numel(caseName) == 1)
                isvalid = ~isempty( this.getIndexOfCasePrivate(caseName) );
            end
        end
        
        function isval = isValidIndex( this, index )
        % Checks whether the given index is valid
        %
        % An index is valid if its > 0 and < this.NumCases
        % 
        % @param index - @type numeric, the index to check for
        % @return isval - @type logical true, if index is valid, false
        %                 otherwise
            isval = index > 0 && index <= this.NumCases;    
        end
        
        function success = setSelectedIndexCallback(this, index)
        % setter for the selected index
        %
        % Errors:  
        %   @error VsCaseControl:setSelectedIndex:invalidIndex - this will be
        %          thrown if index is invalid and no return argument was
        %          defined. The user can prevent the error to be thrown by
        %          checking the output argument
        %
        % @param index - @type numeric the name of the case
        % @return success - @type logical true, if index is valid, and the
        %                   given index was selected
            
            
            ov = this.SelectedIndex;
            success = this.setSelectedIndex(index);
            this.executeCallbackControlSave( this, 'SelectedIndex', ov);
            
        end
        
        function success = setSelectedIndex(this, index)
        % setter for the selected index
        %
        % Errors:  
        %   @error VsCaseControl:setSelectedIndex:invalidIndex - this will be
        %          thrown if index is invalid and no return argument was
        %          defined. The user can prevent the error to be thrown by
        %          checking the output argument
        %
        % @param index - @type numeric the name of the case
        % @return success - @type logical true, if index is valid, and the
        %                   given index was selected
            
            success = false;
            isvalid = this.isValidIndex(index);
            
            if isvalid
                % set the case
                this.setCasePrivate(index);
                success = true;
            else
                
                if nargout < 1
                    % the user is not asking whether this worked out, in
                    % this case we throw an error
                    error( 'VsCaseControl:setSelectedIndex:invalidIndex', ...
                        'Given Case in not in the list');
                end
                
            end
            
        end
        
        
        function success = setSelectedCaseCallback(this, caseName)
        % setter for the selected cases
        %
        % Errors:  
        %   @error VsCaseControl:setSelectedCase:invalidCase - this will be
        %          thrown if index is invalid and no return argument was
        %          defined. The user can prevent the error to be thrown by
        %          checking the output argument
        %
        % @param caseName - @type char the name of the case
        %
        % @return success - @type logical true, if caseName is valid, and the
        %                   given index was selected
        
            ov = this.SelectedCase;
            success = this.setSelectedCase(caseName);
            this.executeCallbackControlSave( this, 'SelectedCase', ov);
            
        end
        
        function success = setSelectedCase(this, caseName)
        % setter for the selected cases
        %
        % Errors:  
        %   @error VsCaseControl:setSelectedCase:invalidCase - this will be
        %          thrown if index is invalid and no return argument was
        %          defined. The user can prevent the error to be thrown by
        %          checking the output argument
        %
        % @param caseName - @type char the name of the case
        %
        % @return success - @type logical true, if caseName is valid, and the
        %                   given index was selected
        
            success = false;
            index = this.getIndexOfCase(caseName);
            if ~isempty(index)
                % set the case
                this.setCasePrivate(index);
                success = true;
            else
                
                if nargout < 1
                    % the user is not asking whether this worked out, in
                    % this case we throw an error
                    error( 'VsCaseControl:setSelectedCase:invalidCase', ...
                        'Given Case in not in the list');
                end
                
            end
            
        end
        
        function settings = getSelectedCaseControlSettings(this)
            settings = this.createPropertyControlSetting('SelectedCase');
        end
        
        function settings = getSelectedIndexControlSettings(this)
            settings = this.createPropertyControlSetting('SelectedIndex');
        end
    end
    
    methods(Access=private)
        
        function setCasePrivate(this, index)
        % private setter for index without checking for errors
        %
        % only for internal use
        
            this.SelectedIndex = index;
            this.SelectedCase  = this.PossibleCases{index};
            this.notifyEnabled( 'SelectionChanged' );
            
        end
        
        function index = getIndexOfCasePrivate(this, caseName)
        % returns the index of the given case 
        %
        % @param caseName - @type char, the name of the case of which to
        %                   return the index
        %
        % @preturn index - @type the index of the caseName, or []
        
            index = find( strcmp( caseName, this.PossibleCases ) );
            if isempty(index)
                index = [];
            end
        end
        
    end
    
    methods(Access ={?vsv.ctrl.mixin.CreatePropertyControlSetting, ...
                               ?vsv.ctrl.AbstractPropertyController})
        
        function setting = createPropertyControlSetting(this, property)
        % create a vsv.ctrl.PropertyControllerSettings for the given
        % property name
        %
        %   the default setting to describe access rights for a given
        %   property is direct property access with post set properties. In
        %   case the inheriting class uses different getters and setters
        %   this can be used to create the settings object that best
        %   describes the access rights for the given property
        %
        % @param property - @type char the name of the property
        % @return setting - @type vsv.ctrl.PropertyControllerSettings
            
            
            persistent selectedCaseTemplate selectedIndexTemplate
            
            switch property
                case 'SelectedCase'
                    
                    if isempty(selectedCaseTemplate)
                    % the settings object
                        selectedCaseTemplate = vsv.ctrl.PropertyControllerSettings(...
                                        'PropertyName',        property , ...
                                        'PropertyType',        'char', ...
                                        'PropertySetter',      'setSelectedCase', ...
                                        'PropertyChangeEvent', 'SelectionChanged', ...
                                        'PropertyLimits',       this.PossibleCases);
                    end
                    setting = selectedCaseTemplate.newCopy();         
                case 'SelectedIndex'
                    if isempty(selectedIndexTemplate)
                    % the settings object
                        selectedIndexTemplate = vsv.ctrl.PropertyControllerSettings(...
                                        'PropertyName',        property , ...
                                        'PropertyType',        'numeric', ...
                                        'PropertySetter',      'setSelectedIndex', ...
                                        'PropertyChangeEvent', 'SelectionChanged', ...
                                        'PropertyLimits',       this.PossibleCases);
                    end
                    setting = selectedIndexTemplate.newCopy();   
                    
                otherwise
                    error('VsToggleButtonControl:createPropertyControlSetting:unknownProperty', ...
                        'Given Property is unknown');
            end
            setting.setPropertyLimits( this.PossibleCases );
            
        end
        
    end
    
    
    methods(Access=protected)
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function success = importSingleStorageParameter(this, parameterName, value) 
            
            
            switch parameterName
                case 'SelectedCase'
                    index = this.getIndexOfCase(value);
                case 'SlectedIndex'
                    index = value;
                otherwise
                    success = false;
            end
            
            % user callbacks can only deal with index!!!
            if index ~= this.SelectedIndex
                success = this.setSelectedIndexCallback( index );
            end
            
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function [paramName, values] = getStorageParameter(this)
        
            paramName = {'SelectedCase'};
            values    = {this.SelectedCase};
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function Class = getSupportedClass(~)
            Class = 'vsv.seq.uicontrol.VsCaseControl'; 
        end
        
        function initVsCaseControl(this, varargin)
        % need to find the possible cases first 
            
            if nargin >= 3
                % nargin >= 3 means we have parameter value 
               index = vsv.util.cellfind( varargin, 'PossibleCases'); 
               if any(index)
                   
                  index = find(index);
                  param = varargin(index);
                  value = varargin(index+1);
                  
                  if numel(param) == 1 
                      
                      this.setPossibleCases( value{1} );
                      varargin( index:index+1 ) = [];
                  else
                      error('VsCaseControl:duplicatePossibleCasesParameter', ...
                        'Missing PossibleCases argument, needed to initialize the object');
                  end
               else
                  % couldn't find PossibleCases nor Labels need to see
                  % whether the user provided NumButtons
                  error('VsCaseControl:missingPossibleCases', ...
                        'Missing PossibleCases argument, needed to initialize the object');
               end
            end
        
            if ~isempty(varargin)
                this.initUserUiControl(varargin{:});
            else
                this.initUserUiControl( );
            end
            
        end
        
        function setPossibleCases(this, posc)
        % setter for possible cases
        %
        % this is private because the possible cases should not change
        % after definition
        %
        % Errors
        %   @error VSButtonGroupContro:setPossibleCases:invalidArgument if
        %          posc is not a cellstring
        %
        % @param posc - possible cases @type cellstring
        
            if iscellstr(posc)  && ~isempty(posc) %#ok
                if numel( unique(posc) ) == numel(posc)
                    this.PossibleCases = posc(:);
                    this.NumCases      = numel(posc);
                    % select the first element to start with
                    this.SelectedIndex = 1;
                    this.SelectedCase  = this.PossibleCases{1};
                else
                    error('VSButtonGroupContro:setPossibleCases:duplicateCases', ...
                    'No duplicate cases allowed as the possible cases')
                end
            else
                error('VSButtonGroupContro:setPossibleCases:invalidArgument', ...
                    'Possible cases must be a cell string')
            end
            
        end
        
        
        function success = applyParamValue(this, param, value)
        %
        %
        
            success = true;
            switch param
                case 'Title'
                    this.setTitle(value);
                case 'SelectedCase'
                    this.setSelectedCase(value);
                case 'SelectedIndex'
                    this.setSelectedIndex(value);
                case 'PossibleCases'
                    this.setPossibleCases(value);
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.UserUiControl(param, value);
            end
            
        end
    
    end
    
    
end

