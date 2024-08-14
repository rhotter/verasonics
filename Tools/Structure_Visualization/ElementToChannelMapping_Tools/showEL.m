function [ ELlist ] = showEL( elNum, apIndex )
% showEL: display mapping of elements to channels, connector pins, and Apod
% arrays

% check for Trans.Connector or Trans.HVMux.Aperture; return with nothing if
% they are not present in the workspace
if evalin('base',  'exist(''Trans'', ''var'') && isfield(Trans, ''Connector'')')
    Trans = evalin('base', 'Trans');
    El2Ch = Trans.Connector;
    aperstr = ''; % not an HVMux probe, so aperture value is irrelevant
elseif evalin('base',  'exist(''Trans'', ''var'') && isfield(Trans, ''HVMux'') && isfield(Trans.HVMux, ''Aperture'')')
    Trans = evalin('base', 'Trans');
    % this is an HVMux probe, so need to find aperture index requested or
    % default to 1 otherwise
    if nargin == 2
        % limit the requested value to the range defined by the size of the
        % Aperture array
        apIndex = max(1, min(size(Trans.HVMux.Aperture, 2), apIndex));
    else
        apIndex = 1;
    end
    El2Ch = Trans.HVMux.Aperture(:, apIndex);
    aperstr = ['  aperture ', num2str(apIndex)];
else
    ELlist = 'Both Trans.Connector and Trans.HVMux.Aperture not found';
    return
end

apod2El = find(El2Ch);
El2apod = zeros(Trans.numelements, 1);
El2apod(apod2El) = 1:nnz(El2Ch);

if evalin('base',  'exist(''UTA'', ''var'') && isfield(UTA, ''ChPinNames'')')
    UTA = evalin('base', 'UTA');
    pinNames = 1;
else
    pinNames = 0;
end


if nargin < 1 || isempty(elNum)
    elNum = (1:Trans.numelements);
end
elStart = min(elNum);
elStart = max(elStart, 1);
elEnd = max(elNum);
elEnd = min(elEnd, Trans.numelements);
if elStart>elEnd
    elStart = 1;
    elEnd = Trans.numelements;
end

ELlist = cell((elEnd-elStart+1), 1);
for i = elStart:elEnd
    if El2Ch(i) == 0
        chStr = '  not used';
        pinStr = '';
        apodstr = '';
    else
        chStr = ['   CH ', num2str(El2Ch(i), '%03d')];
        apodstr = ['  APOD(', num2str(El2apod(i),'%03d'), ')', aperstr];
        if pinNames
            pinStr = ['   pin ', UTA.ChPinNames{El2Ch(i), 2}];
        else
            pinStr = '';
        end
    end
    ELlist{i-elStart+1, 1} = ['EL ', num2str(i, '%03d'), chStr, pinStr, apodstr];
end

end

