function resistanceValues = convertThermistorValues(V)
% Convert thermistor A/D values to resistance values in Ohms
%
% Parameters:
%   @param V - an array of numeric values
%
% Return value:
%   @return resistanceValues - @type double, a matrix the same size
%               as the input, containing resistance values in Ohms
%
% Copyright 2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

    % convert A/D readings into resistance values using the formula
    %     Ohms = (v * 2490) / 4095 - v)
    resistanceValues = zeros(1,length(V));
    for i=1:length(V)
        clippedV = min(V(i), 4094); % clip at 4094 so we don't get infinity
        resistanceValues(i) = clippedV * 2490 / (4095 - clippedV);
    end
end