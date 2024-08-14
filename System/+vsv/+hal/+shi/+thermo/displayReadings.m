function displayReadings(labelHandles, selectedConnector, thermistorValues, resistances)
% Display thermistor values read from the hardware in a human-readable way
%
% Parameters:
%   @param labelHandles - an array of previously-created uilabel handles
%   @param selectedConnector - the integer scanhead connector number
%   @param thermistorValues - readings previously collected by the
%                                readProbeThermistorValues function
%   @param resistances - resistance values corresponding to the thermistorValues
%
% Return values: none
%
% Copyright 2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

    numThermistors = numel(thermistorValues);
    
    
    if ~isempty(labelHandles) && numel(labelHandles) >= numThermistors && numThermistors == numel(resistances)
        
        formatStr = '[Transducer %d]  temperature sensor %d:   measured A/D value: %.0f, resistance %.0f Ohms';
        for i=1:numThermistors
            if ishandle(labelHandles(i))
                s = sprintf(formatStr, selectedConnector, i, thermistorValues(i), resistances(i) );
                set(labelHandles(i),'Text', s);
            end
        end
        
    end
    
end