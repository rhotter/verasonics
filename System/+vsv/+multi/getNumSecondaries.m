function numSecondaries = getNumSecondaries()
% get the number of secondary nodes connected based on the static IP addresses
[~,resultStr]=unix('ip -4 addr');
out = regexp(resultStr,'(?<=inet\s10\.10\.)\d+(\.\d+){1}','tokens');
numSecondaries = length(out);
