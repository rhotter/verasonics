function thermistorValues = readProbeThermistorValues(selectedConnector, numThermistors, doUseHardware)
% Retrieve thermistor values from the probe on the indicated connector
%
% Parameters:
%   @param selectedConnector - the integer scanhead connector number
%
% Return value:
%   @return thermistorValues - @type integer array containing thermistor A/D readings
%               (values of -1 indicate read errors)
%
% Copyright 2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
    
    if nargin < 3 || isempty(doUseHardware)
        doUseHardware = true;
    end

    import com.verasonics.hal.shi.Shi
    thermistorValues = zeros(1,numThermistors); % preallocate an arry of the right size
    
    for thermistor = 1:numThermistors
        
        if doUseHardware
            thermistorInfo = Shi.getProbeThermistorInfo(selectedConnector, thermistor-1);
            if isempty(thermistorInfo)
                fprintf(2, 'Error from getProbeThermistorInfo() for thermistor #%d\n',thermistor)
                tempSenseAD = -1;
            else
                tempSenseAD = thermistorInfo.currentValue;
            end
        else
            % generate fake A/D value that fall into the legal range
            tempSenseAD = randi([0,4094]);
        end
        thermistorValues(thermistor) = tempSenseAD;
    end

end