function [ CHlist ] = showCH( chNum, apIndex )
% showCH: display mapping of channels to elements and connector pins
%
% Syntax:  CHlist = showCH( chNum, apIndex );
%
%       [chNum]     (optional) a range of channels to display [default is 1:UTA.numch]
%       [apIndex]   (optional) aperture index for HVMux probes [default is 1]
%
%       CHlist      cell array containing rows of strings listing 'CH Ch# EL El# pin PinLabel APOD(index)'
%
% Example: showCH(1:2)
%          'CH 001  EL 076  pin U5  APOD(076)'
%          'CH 002  EL 045  pin U4  APOD(045)'
%


if evalin('base',  '~exist(''UTA'', ''var'') || ~isfield(UTA, ''ChPinNames'')')
    CHlist = 'UTA.chmap not found';
    return
end
UTA = evalin('base', 'UTA');
if nargin < 1 || isempty(chNum)
    chNum = (1:UTA.numCh);
end
chStart = min(chNum);
chStart = max(chStart, 1);
chEnd = max(chNum);
chEnd = min(chEnd, UTA.numCh);
if chStart>chEnd
    chStart = 1;
    chEnd = UTA.numCh;
end

if evalin('base',  'exist(''Trans'', ''var'') && isfield(Trans, ''Connector'')')
    Trans = evalin('base', 'Trans');
    El2Ch = Trans.Connector;
    [apod2El, ~, apod2Ch] = find(El2Ch);
    Ch2El = zeros(UTA.numCh, 1);
    Ch2El(apod2Ch) = apod2El;
    ELmap = 1;
    aperstr = '';
elseif evalin('base',  'exist(''Trans'', ''var'') && isfield(Trans, ''HVMux'') && isfield(Trans.HVMux, ''Aperture'')')
    Trans = evalin('base', 'Trans');
    if nargin == 2
        apIndex = max(1, min(size(Trans.HVMux.Aperture, 2), apIndex));
    else
        apIndex = 1;
    end
    El2Ch = Trans.HVMux.Aperture(:, apIndex);
    [apod2El, ~, apod2Ch] = find(El2Ch);
    Ch2El = zeros(UTA.numCh, 1);
    Ch2El(apod2Ch) = apod2El;
    ELmap = 1;
    aperstr = ['  aperture ', num2str(apIndex)];
else
    ELmap = 0;
end

apodL = length(apod2Ch);
Ch2apod = zeros(UTA.numCh, 1);
Ch2apod(apod2Ch) = 1:apodL;



CHlist = cell((chEnd-chStart+1), 1);

if ELmap
    for i = chStart:chEnd
        if Ch2El(i) == 0
            Elstr = ' not used';
            apodstr = '';
        else
            Elstr = ['  EL ', num2str(Ch2El(i),'%03d'), ' '];
            apodstr = ['  APOD(', num2str(Ch2apod(i),'%03d'), ')', aperstr];
        end
        CHlist{i-chStart+1, 1} = [UTA.ChPinNames{i, 1}, Elstr, ' pin ', UTA.ChPinNames{i, 2}, apodstr];
    end
else
    for i = chStart:chEnd
        CHlist{i-chStart+1, 1} = [UTA.ChPinNames{i, 1}, ' pin ', UTA.ChPinNames{i, 2}];
    end
end

end

