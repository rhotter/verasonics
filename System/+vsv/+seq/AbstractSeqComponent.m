classdef AbstractSeqComponent < vsv.events.EnableEventsMixin
%ABSTRACTSCRIPTCOMPONENT a component in a ultrasound sequence
%   
% A seq component is an object with settings that define a component of
% an ultraound sequence. Examples are PData, TX, and other objects. 
%
% Each seq component is usally saved in an array at the lower seq levels.
% The index of the object location in the list (i.e., the pointer to that
% object) is stored in this object with the ListIndex. 
%
% The class defines a PropertyChanged event that will indicate that a
% property of the class has changed. This will send an event with the
% property name. The subclass can use notifyPropertyChanged('PropertyName')
% to notify listeners. 
% Although, matlab provides a mechanism for listening to property events
% this method can have some advantages:
%   - Matlabs Property Events are "one Property -> many listener" approach
%   - This approach is "all Properties -> many listeners" 
%
% The components need to communicate changes of their properties to other
% VSX components. For example if a property is changing, runAcq needs to be
% informed. To give the user the experience to work with the seq.
% components directly but keeping the sequence components light weight, the
% listener approach seems to be suitable. In this case another application
% component can register a listener to changes of properties and perform
% certain actions. 
%
% However, with the original approach for each property, we would need to
% install a seperate listener that may has to call the same function. This
% can cause an overhead if you need to dynamically register listeners
% during runtime. In such cases the "all property -> many listener" is a
% better approach because it only requires to register one listener for
% changes of all proeprties. 
%
%
% Version 1.0 | 2020-04-22 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
       
    properties(SetAccess=protected)
        
        % the name of the component @type char, this will match the name as
        % used in VSX for defining the sequence component such as TX, TW, 
        CompName (1,:) char = '';
        
    end
    
    properties(SetAccess={?vsv.seq.sets.AbstractSeqCompSet})
    
        % @type int16 this is the reference index in the sequence component
        % array. The index is used as sort of a pointer and could also be
        % understood as an ID. Currently, this can only be assigned by the
        % sequence component set, which is the container for storing
        % sequence components and should define the index in the container
        ListIndex (1,1) int16 = 1;
        
        % @type logical inidcates whether this component is part of a list 
        IsPartOfList (1,1) logical = false;
        
    end
    
    properties(Access=private)
        
        % @type vsv.events.PropertyChangeEvent a cache for event objects
        % recreating events can be slow, its much faster when reusing the
        % events, then 
        CachedEvent (:,1) vsv.events.PropertyChangeEvent {vsv.util.isScalarOrEmpty} = vsv.events.PropertyChangeEvent.empty;
    end
    
    events
        
        % Indicates that one of the properties has changed, sends a event 
        % @type vsv.events.PropertyChangeEvent
        PropertyChanged;
        
    end
    
         
    %% abstract methods that need to be implemented by the sub class
    methods(Abstract)
        
        % create a new instance of this component
        %  this will return a new pointer to a new data object
        %
        % @return @type vsv.seq.AbstractSeqComponent
        ns = newInstance(this);
        
        % @signature ns = newInstance(this);
        % List the properties that can be imported
        %
        %   This are mostly the properties that are defined in the original
        %   structs used for the sequence components.
        %   
        % @signature [ props, values] = listImportProperties(this);  
        %
        % @return props - @type cell string a list of properties that can
        %                 be imported
        % @return values - @type cell array with values that fit to the
        %                 properties
        [ props, values] = listImportProperties(this);
        
    end
    
    %% public methods preimplemented for all subclasses
    methods
        
        function this = AbstractSeqComponent(compName)
        % creates a new seq component based on the component name
        %   The component name is required to define the component. For
        %   each subclass this should be unique. 
        %
        % @param compName - @type char, the name of the component, e.g.,
        %                   'TX', 'Pdata'
        
            this.setCompName(compName);
        end
        
        function str = exportStruct(this)
        % Will export this object to a struct that matches the struct
        % component in the main workspace
        %   
        %   This uses this.listImportProperties to export the properties to
        %   the fields listed by that function.
        %
        % @return str - @type struct, a represenation of this object as a
        %               struct, which will match the struct 
        
            [ props, vals ] = this.listImportProperties(); 
            str = cell2struct( vals(:), props(:), 1);
        end
        
        function importStruct(this, str, fields)
        % imports a given struct, str, to this object
        % 
        % @param str - @type struct the struct to import
        % @oaram fields - @type cell string with field names to import.
        %                 @optional, @default fieldnames(str)
        
            props = this.listImportProperties();
            
            if nargin < 3
                fields = fieldnames(str);
            end
        
            % only get common fields
            setFields = intersect( props, fields);
            nFields = numel(setFields);
            % assign all fields and data from struct to this object
            for i = 1:nFields
                fieldi = setFields{i};
                this.(fieldi) = str.(fieldi);
            end
            
        end
        
        function isp = isPartOfList(this)
        % returns true if this sequence component is connected to a list. 
        % @return isp - @type logical true if this component is part of a
        %               list.
            isp = this.IsPartOfList;
        end
        
        function li = getListIndex(this)
        % returns the list index of the sequence component in the lower
        % level sequence list
        % @return li - @type numeric the list index 
            li = this.ListIndex;
        end
        
        function name = getCompName(this)
        % return the name of the component
        % @return name - @type char the name of the components
            name = this.CompName;
        end
        
    end
    
    
    %% protected methods
    methods(Access=protected)
        
        function evt = getPropertyEvent(this, property)
        % returns the property event 
        % @return evt - @type vsv.events.PropertyChangeEvent
        
            if isempty(this.CachedEvent)
                this.CachedEvent = vsv.events.PropertyChangeEvent;
            end
            this.CachedEvent.PropertyName = property;
            evt = this.CachedEvent;
        end
        
        function notifyPropertyChanged(this, property)
        % notify that a property has changed
        %
        %
        % @return evt - @type vsv.events.PropertyChangeEvent
            if this.EventsEnabled
                evt = this.getPropertyEvent(property);
                this.notifyEnabled( 'PropertyChanged', evt );
            end
            
        end
                
        function setCompName(this, name)
        % setter for the name
        % @param name - @type char the name of the component
        
            [name, msg] = vsv.mixin.HasName.validateName(name);
            if isempty(msg)
                this.CompName = name;
            else
                error( 'AbstractSeqComponent:setCompName:invalidArgument', ...
                    [ 'Invalid name: ' name]);
            end
        end
        
        
    end
    
end

