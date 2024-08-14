function generateSecondaryEcho(mainGUIFigure)
% vsv.multi.generateSecondaryEcho
%
% Generates echo commands for the Secondary GUIs
%

% Copyright 2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.

objList = mainGUIFigure.Children;

for a =1:length(objList)
    if isa(objList(a),'matlab.ui.control.UIControl')
        callback = objList(a).Callback;
        if strcmp(objList(a).Tag,'FreezeButton')
            objList(a).Callback = @(src,event)cellfun(@(x)feval(x,src,event),...
                {callback,...
                @(src,event)doEchoFreeze(objList(a))});
        else
            if iscell(callback) %Any Control requiring input arguments (vsx_gui sliders for example)
                newCallback  = callback{1};
                arguments = callback(2:end);
                doEchoCbk = @(src,event,varargin)doEcho(src,event,objList(a));
                callbackList = {newCallback, doEchoCbk };
                newCallback = @(src,evt)combinedCallback(src, evt, callbackList, arguments{:});
                objList(a).Callback = newCallback;
                
            else %all other controls
                objList(a).Callback = @(src,event)cellfun(@(x)feval(x,src,event), ...
                    {callback, ...
                    @(src,event)doEcho(src,event,objList(a))});
            end
        end
    end
    % Explicitly setup the closeRequest Function
    
    mainGUIFigure.CloseRequestFcn = @(src,event)doEchoclose(mainGUIFigure);
end
return

function doEcho(~,~,obj)
% Sometimes it is defined by tag and sometimes it is defined by string
if ~isempty(obj.Tag)
    identifier = 'tag';
    identifierString = obj.Tag;
else
    identifier = 'string';
    identifierString = obj.String;
end

if strcmp(obj.Style,'edit')
    SecondaryCmd = sprintf('set(findobj(''%s'',''%s''),''String'',%s)', identifier, identifierString, obj.String);
else
    SecondaryCmd = sprintf('set(findobj(''%s'',''%s''),''Value'',%i)', identifier, identifierString, obj.Value);
end
SecondaryCmd2 = sprintf('CB = get(findobj(''%s'',''%s''),''Callback'')', identifier, identifierString);

if contains(func2str(obj.Callback),'combinedCallback') %hack
    SecondaryCmd3=sprintf('feval(CB{1},findobj(''%s'',''%s''), 0, CB{2})', identifier, identifierString);
else
    SecondaryCmd3=sprintf('feval(CB,findobj(''%s'',''%s''), 0)', identifier, identifierString);
end

remoteCmd = ['  ' SecondaryCmd ';' SecondaryCmd2 ';' SecondaryCmd3]; %concatenate command
fprintf('remoteCmd: %s\n',remoteCmd);   % for debugging
assignin('base','remoteCmd',remoteCmd); % command gets evaluated in VSX.m
com.verasonics.common.util.Prop.setProperty('remoteCmd',remoteCmd);
return

function doEchoFreeze(obj)
syncSocketId=evalin('base','syncSocketId');
RDMAconfig = evalin('base','RDMAconfig');
remoteCmd = sprintf('freeze=%i',obj.Value);
for a = 1:RDMAconfig.numSecondaryNodes
    socket(a) = com.verasonics.common.socket.Net.getSocketById(syncSocketId(a));
    if obj.Value
        socket(a).writeString('F');
    else
        socket(a).writeString('C');
    end
    socket(a).readString();%response from Secondary that it is ready to receive a command
    socket(a).writeString(remoteCmd);
end
return

function doEchoclose(obj)
try
    syncSocketId=evalin('base','syncSocketId');
    RDMAconfig = evalin('base','RDMAconfig');
    remoteCmd = 'close(findobj(''Tag'',''UI''));freeze=0;';
    for a = 1:RDMAconfig.numSecondaryNodes
        socket(a) = com.verasonics.common.socket.Net.getSocketById(syncSocketId(a));
        socket(a).writeString('C');
        socket(a).readString();%response from Secondary that it is ready to receive a command
        socket(a).writeString(remoteCmd); %send the close GUI command
    end
    disp('Closed remote systems')
    assignin('base', 'vsExit', 1);
    if exist('hvtmr','var'), stop(hvtmr); delete(hvtmr); end
    delete(obj);
catch
    delete(obj);
    disp('Cannot Close Remote Systems');
end
return

% Execute combined callbacks
function combinedCallback(src, evt, callbacks, varargin)
nCallbacks = numel(callbacks);
for i = 1:nCallbacks
    callbacki = callbacks{i};
    callbacki( src, evt, varargin{:});
end
return