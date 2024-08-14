function ColorMapTool(varargin)

persistent originalMap gammaCurveMenu

% ColorMapTool is used to adjust the gamma curve for colormap adjustment.
% Please note that only the gray scale (lower half) of the doppler script
% will be adjusted

% Create the figure window.
ScrnSize = get(0,'ScreenSize');
width = 350; % width
height = 200; % height
fWidth = width + 100;
fHeight = height + 100; % figure height, including border.
bkgrnd = [0.9 0.9 0.9];

hg = findobj('tag','ColorMapTool');
if ishandle(hg)
    pos = get(hg,'Position');
    delete(hg);
else
    pos = [(ScrnSize(3)-fWidth)/2 (ScrnSize(4)-fHeight)/2 fWidth fHeight];
end

handle.Figure = figure('Visible','on',...
    'Position',pos,...
    'Name','Color Map Tool',...
    'Color',bkgrnd,...
    'NumberTitle','off',...
    'MenuBar','none',...
    'Resize','on',...
    'tag','ColorMapTool');
set(handle.Figure,'CloseRequestFcn',{@closefunc});
set(handle.Figure,'DefaultUicontrolBackgroundColor',bkgrnd)

% default values
winNum = 1;
xCoor = linspace(0,1,256)';
yCoor = xCoor;

passCheck = 0;
%% Create dropdown menu for display window, curve selection, and save as..
if evalin('base','exist(''Resource'',''var'')')
    Resource = evalin('base','Resource');
    if isfield(Resource,'DisplayWindow') && isfield(Resource.DisplayWindow,'figureHandle')
        if ishandle(Resource.DisplayWindow(winNum).figureHandle)
            passCheck = 1;
        end
    end
end

if passCheck
    nWindow = length(Resource.DisplayWindow);
    if ~isfield(Resource.DisplayWindow(winNum),'splitPalette')||isempty(Resource.DisplayWindow(winNum).splitPalette)
        sp = 0;
    else
        sp = Resource.DisplayWindow(winNum).splitPalette;
    end
else
    disp('No display window or Resource is not shown in the workspace.');
    delete(handle.Figure);
    return
end

for j = 1:nWindow
    WinStr{j} = [blanks(10),num2str(j,1)];
end

dispWindowTxt = uicontrol('Style','text',...
    'String','Display Window',...
    'Units','normalized',...
    'Position',[0.7 0.75 0.25 0.05],...
    'FontUnits','normalized',...
    'FontSize',0.85,...
    'FontWeight','bold');

dispWindowMenu = uicontrol('Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.7 0.65 0.25 0.07],...
    'String',WinStr,...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Callback',@windowSel);

% load default curves
setting = load('defaultGamma');
Gamma = setting.Gamma;

nCurve = length(Gamma);
CurveStr = cell(nCurve+1,1);
CurveStr{1} = '   Original ';
for j = 1:nCurve
    CurveStr{j+1} = ['   Curve ',num2str(j,1)];
end

CurveStrCustom = CurveStr;
CurveStrCustom{end+1} = '   Custom ';

gammaCurveTxt = uicontrol('Style','text',...
    'String','Gamma Curve',...
    'Units','normalized',...
    'Position',[0.7 0.55 0.25 0.05],...
    'FontUnits','normalized',...
    'FontSize',0.85,...
    'FontWeight','bold');

gammaCurveMenu = uicontrol('Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.7 0.45 0.25 0.07],...
    'String',CurveStr,...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Callback',@CurveSel);

% save new curve as default
SaveStr = CurveStr(2:6);

saveCurveTxt = uicontrol('Style','text',...
    'String','Save Curve to...',...
    'Units','normalized',...
    'Position',[0.7 0.35 0.25 0.05],...
    'FontUnits','normalized',...
    'FontSize',0.85,...
    'FontWeight','bold');

saveCurveMenu = uicontrol('Style','popupmenu',...
    'Units','normalized',...
    'Position',[0.7 0.25 0.25 0.07],...
    'String',SaveStr,...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Callback',@SaveSel);

%% Creat axes for gamma curve and 3 markers
handle.Axes = axes(...
    'Parent',handle.Figure,...
    'Units','normalized',...
    'Position',[0.1,0.15,0.55,0.725],...
    'NextPlot','replacechildren',...
    'FontSize',14,...
    'Color',[0.95 0.95 0.95],...
    'XLim',[0 1],...
    'YLim',[0 1],...
    'Tag','gammaAxes');

% ColorMapTool will be opened with the following conditions
% 1) If the customGamma is not the workspace -> Original Colormap
% 2) If the customGamma is in the workspace, either the same as one of default
%    curves or after previous adjustment. Load it back!

cMap = Resource.DisplayWindow(winNum).Colormap;
finalMap = cMap; % make finalMap be seen by other nested functions

if evalin('base','exist(''customGamma'',''var'')')
    customGamma = evalin('base','customGamma');
    for j = 1:5
        if isequal(Gamma(j).Curve,customGamma.Curve)
            curveNum = j;
            break;
        else
            curveNum = 6; % if not the default curvers, use custom
            set(gammaCurveMenu,'String',CurveStrCustom);
        end
    end
    curveApplied = 1;
    set(saveCurveMenu,'enable','on');
else
    curveNum = 0; % no customGamma? --> Original Colormap
    originalMap = cMap;
    curveApplied = 0;
    set(saveCurveMenu,'enable','off');
end
set(gammaCurveMenu,'Value',curveNum+1);

% if calling after loading preSet, a customGamma might exist but original
% Map was not be saved. In this case, the originalMap is set to gray(256)
if isempty(originalMap)
    originalMap = gray(256);
    if isequal(cMap(128),1) % Doppler script
        originalMap(1:128,:) = originalMap(1:2:256,:);
        originalMap(129:256,:) = cMap(129:256,:);
    end
end

if isequal(cMap(1:128,1),cMap(1:128,2),cMap(1:128,3))
    % if yes, colorMap is gray
    isgray = 1;
else
    % if not, cMap is in other colorMap
    isgray = 0;
end

% set marker limit
xlim = 0;
ylim = 1;
delta = 0.01;

applyCurve;

%% all callbacks

    function windowSel(hObject,~)
        winNum = get(hObject,'Value');
        if ~isfield(Resource.DisplayWindow(winNum),'splitPalette')||isempty(Resource.DisplayWindow(winNum).splitPalette)
            sp = 0;
        else
            sp = Resource.DisplayWindow(winNum).splitPalette;
        end
        dispUpdate;
    end

    function CurveSel(hObject,~)
        curveNum = get(hObject,'Value')-1;
        if curveNum > 0
            set(saveCurveMenu,'enable','on');
        else
            set(saveCurveMenu,'enable','off');
        end
        applyCurve;
        dispUpdate;
    end

    function SaveSel(hObject,~)
        newCurveNum = get(hObject,'Value');
        location = which('defaultGamma.mat');
        Gamma(newCurveNum).Curve = get(handle.Line,'YData')';
        Gamma(newCurveNum).Marker(:,1) = cell2mat(get(handle.Marker,'XData'));
        Gamma(newCurveNum).Marker(:,2) = cell2mat(get(handle.Marker,'YData'));
        save(location,'Gamma');
        set(gammaCurveMenu,'Value',newCurveNum+1);
    end

    function applyCurve
        switch curveNum
            case 0 % back to original map
                cla
                finalMap = originalMap;
                if evalin('base','exist(''customGamma'',''var'')')
                    evalin('base','clear customGamma;');
                end
                text(0.1,0.5,'Please select a Gamma Curve','FontSize',12);
                set(gammaCurveMenu,'String',CurveStr);
                curveApplied = 0;
            case {1,2,3,4,5}
                yCoor = Gamma(curveNum).Curve;
                markerX = Gamma(curveNum).Marker(:,1);
                markerY = Gamma(curveNum).Marker(:,2);
                set(gammaCurveMenu,'String',CurveStr);
                curveApplied = 1;
            case 6
                if evalin('base','exist(''customGamma'',''var'')')
                    customGamma = evalin('base','customGamma');
                    yCoor = customGamma.Curve;
                    markerX = customGamma.Marker(:,1);
                    markerY = customGamma.Marker(:,2);
                else
                    text(0.1,0.5,'Please select a Gamma Curve','FontSize',12);
                    return
                end
                curveApplied = 1;
        end

        if curveApplied
            % set marker limit
            xlim = markerX(1);
            ylim = markerY(end);

            handle.Line = plot(handle.Axes,xCoor,yCoor,'-b','LineWidth',1.5); grid on; hold on;
            for n = 1:5
                handle.Marker(n) = plot(markerX(n),markerY(n),'o',...
                    'MarkeredgeColor','k',...
                    'MarkerfaceColor','k',...
                    'MarkerSize',6);
                set(handle.Marker(n),'ButtonDownFcn',{@MarkerButtonDown,n});
            end
            hold off;

            % have customGamma curve in the workspace
            customGamma.Curve = yCoor;
            customGamma.Marker(:,1) = markerX;
            customGamma.Marker(:,2) = markerY;
            assignin('base','customGamma',customGamma);
        end

    end

    function MarkerButtonDown(~,~,nMarker)
        % Set the figure motion function accordingly
        set(handle.Figure,'WindowButtonMotionFcn',{@moveMarker,nMarker});
        set(handle.Figure,'WindowButtonUpFcn',@MarkerButtonUp);
    end

    function MarkerButtonUp(~,~)
        set(handle.Figure,'WindowButtonMotionFcn',[]);
        set(handle.Figure,'WindowButtonUpFcn',[]);
    end

    function moveMarker(~,~,nMarker)
        hMarker = handle.Marker(nMarker);
        X(:,1) = get(handle.Marker,'XData'); X = cell2mat(X);
        Y(:,1) = get(handle.Marker,'YData'); Y = cell2mat(Y);

        % update marker
        cp = get(handle.Axes, 'CurrentPoint');
        x  = cp(1,1);
        y  = cp(1,2);

        % Is marker outside the handle.Axeses limits?
        switch nMarker
            case 1 % y value will be always zero
                y = 0;
                if x <= 0 || x >= X(2)
                    x =  min(max([0+delta x]),X(2)-delta);
                end
                xlim = x;
            case 2
                if y <= Y(1) || y >= Y(3) || x <= X(1) || x >= X(3)
                    x = min(max([X(1)+delta x]),X(3)-delta);
                    y = min(max([Y(1)+delta y]),Y(3)-delta);
                end
            case 3
                if y <= Y(2) || y >= Y(4) || x <= X(2) || x >= X(4)
                    x = min(max([X(2)+delta x]),X(4)-delta);
                    y = min(max([Y(2)+delta y]),Y(4)-delta);
                end
            case 4
                if y <= Y(3) || y >= Y(5) || x <= X(3) || x >= X(5)
                    x = min(max([X(3)+delta x]),X(5)-delta);
                    y = min(max([Y(3)+delta y]),Y(5)-delta);
                end
            case 5  % x value will be always 1
                x = 1;
                if y <= Y(4) || y >= 1
                    y = min(max([Y(4)+delta y]),1-delta);
                end
                ylim = y;
        end
        set(hMarker,'XData',x,'YData',y);

        % update display window
        markerX = cell2mat(get(handle.Marker,'XData'));
        markerY = cell2mat(get(handle.Marker,'YData'));

        % update curve based on new marker
        yCoor = pchip(markerX,markerY,xCoor);
        yCoor(xCoor<xlim) = 0; yCoor(yCoor>ylim)=ylim;
        set(handle.Line,'XData',xCoor,'YData',yCoor);

        % have customGamma curve in the workspace once the curved is moved
        set(gammaCurveMenu,'String',CurveStrCustom);
        set(gammaCurveMenu,'Value',7);
        customGamma.Curve = yCoor;
        customGamma.Marker(:,1) = markerX;
        customGamma.Marker(:,2) = markerY;
        assignin('base','customGamma',customGamma);
        dispUpdate;
    end

    function dispUpdate
        if curveApplied
            if isgray
                if sp == 1
                    finalMap(1:128,:) = [yCoor(1:2:256) yCoor(1:2:256) yCoor(1:2:256)];
                else
                    finalMap = [yCoor yCoor yCoor];
                end
            else % if not graymap, RGB has to be modified
                if evalin('base','exist(''customGamma'',''var'')') % if no customGamma, just apply finalMap
                    customGamma = evalin('base','customGamma');
                    markerX = customGamma.Marker(:,1);
                    markerY = customGamma.Marker(:,2);
                    for n = 1:3
                        buffer = pchip(markerX,markerY,cMap(:,n));
                        buffer(cMap(:,n)<xlim) = 0;
                        buffer(buffer>ylim) = ylim;
                        finalMap(:,n) = buffer;
                    end
                end
            end
        end
        set(handle.Axes,'FontSize',14,'Color',[0.95 0.95 0.95],'XLim',[0 1],'YLim',[0 1]);
        assignin('base','newMap',finalMap);
        evalin('base',['Resource.DisplayWindow(',num2str(winNum),').Colormap = newMap;']);
        Control.Command = 'set&Run';
        Control.Parameters = {'DisplayWindow',winNum,'colormap',finalMap};
        assignin('base','Control',Control);
    end

    function closefunc(~,~)
        delete(handle.Figure);
        set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
    end

end
