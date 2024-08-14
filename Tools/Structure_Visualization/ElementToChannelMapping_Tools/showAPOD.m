function [ Apodlist ] = showAPOD( apodNum, apIndex )
% showAPOD: display mapping of Apod array entries to elements, channels, and connector pins

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
        % limit the requested aperture index value to the range defined by
        % the size of the Aperture array
        apIndex = max(1, min(size(Trans.HVMux.Aperture, 2), apIndex));
    else
        apIndex = 1;
    end
    El2Ch = Trans.HVMux.Aperture(:, apIndex);
    aperstr = ['  aperture ', num2str(apIndex)];
else
    Apodlist = 'Both Trans.Connector and Trans.HVMux.Aperture not found';
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

sizeApod = nnz(El2Ch);

if nargin < 1 || isempty(apodNum)
    apodNum = (1:sizeApod);
end
apodStart = min(apodNum);
apodStart = max(apodStart, 1);
apodEnd = max(apodNum);
apodEnd = min(apodEnd, sizeApod);
if apodStart>apodEnd
    apodStart = 1;
    apodEnd = sizeApod;
end

Apodlist = cell((apodEnd-apodStart+1), 1);
for i = apodStart:apodEnd
    chStr = ['   CH ', num2str(El2Ch(apod2El(i)), '%03d')];
    if pinNames
        pinStr = ['   pin ', UTA.ChPinNames{El2Ch(apod2El(i)), 2}];
    else
        pinStr = '';
    end
    Apodlist{i, 1} = ['APOD(', num2str(i, '%03d'), ')   EL ', num2str(apod2El(i), '%03d'), chStr, pinStr, aperstr];
end

end

