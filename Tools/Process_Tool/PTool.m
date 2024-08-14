function PTool
%IPTool - Tool for manipulating image processing parameters.
%   Call before running an imaging script. When an attribute is changed,
%   the new value is set for the first Process structure found with
%   the imageDisplay method.

% Create the figure window.
ScrnSize = get(0,'ScreenSize');
width = 350; % width
height = 400; % height
fWidth = width + 100;
fHeight = height + 100; % figure height, including border.
bkgrnd = [0.8 0.8 0.8];

hf = findobj('tag','ProcessTool');
if ishandle(hf)
    pos = get(hf,'Position');
    delete(hf);
else
    pos = [(ScrnSize(3)-fWidth)/2 (ScrnSize(4)-fHeight)/2 fWidth fHeight];
end

hf = figure('Visible','off',...
           'Position',pos,...
           'Name','Process Tool',...
           'NumberTitle','off',...
           'MenuBar','none',...
           'Color',bkgrnd,...
           'Resize','on',...
           'tag','ProcessTool');
set(hf,'CloseRequestFcn',{@closefunc});
set(hf,'DefaultUicontrolBackgroundColor',bkgrnd)
axis([0 1 0 1]);
axis off;
line([0 1],[0.92 0.92],'Color','b','LineWidth',2);
% Define attribute sizes.
AT = struct('AS',[0.3 0.057],...
            'AFS',0.50);
% Define slider group offsets and sizes. All units are normalized.
SG = struct('TO',[0.0 0.0],...     % title offset
            'TS',[0.3,0.057],...   % title size
            'TF',0.55,...           % title font size
            'SO',[0.34,0.015],...   % slider offset
            'SS',[0.25,0.04],...    % slider size
            'EO',[0.25,0.02],...    % edit box offset
            'ES',[0.08,0.037]);     % edit box size

% Create the Update button, Process selector and Process info.
uicontrol(hf,'Style','pushbutton',...
            'String','Update',...
            'Units','normalized',...
            'Position',[0.05 0.91 0.15 0.054],...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'BackgroundColor',bkgrnd+0.07,...
            'Callback',{@update_Callback});
uicontrol(hf,'Style','text',...
            'String','Process',...
            'Units','normalized',...
            'Position',[0.23 0.89 0.17 0.06],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'FontWeight','bold');
hPNum = uicontrol(hf, 'Style','popupmenu',...
            'String',{'0'},...
            'Units','normalized',...
            'Position',[0.4 0.90 0.13 0.05],...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','processNum',...
            'Callback',{@processNum_Callback});
uicontrol(hf,'Style','text',...
            'String','classname:',...
            'Units','normalized',...
            'Position',[0.55 0.91 AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold');
hClassImage = uicontrol(hf,'Style','text',...
            'String','Image',...
            'Units','normalized',...
            'Position',[0.78 0.91 AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS);
uicontrol(hf,'Style','text',...
            'String','method:',...
            'Units','normalized',...
            'Position',[0.55 0.87 AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold');
hMethodID = uicontrol(hf,'Style','text',...
            'String','imageDisplay',...
            'Units','normalized',...
            'Position',[0.69 0.87 AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS);
uicontrol(hf,'Style','text',...
            'String','Parameters:',...
            'Units','normalized',...
            'Position',[0.05 0.77 0.25 0.06],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'FontWeight','bold');

% Create the UIControl objects for class Image, method imageDisplay.
% -- imgbufnum
Pos = [0.3 0.77];
hImgdsply(1) = uicontrol(hf,'Style','text',...
            'String','imgbufnum',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(2) = uicontrol(hf,'Style','text',...
            'String','1',...
            'Units','normalized',...
            'Position',[Pos(1)+0.25 Pos(2) AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'Visible','off');
% -- framenum
Pos = [0.3 0.73];
hImgdsply(3) = uicontrol(hf,'Style','text',...
            'String','framenum',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(4) = uicontrol(hf,'Style','text',...
            'String','-1',...
            'Units','normalized',...
            'Position',[Pos(1)+0.24 Pos(2) AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'Visible','off');
% -- pdatanum
Pos = [0.3 0.69];
hImgdsply(5) = uicontrol(hf,'Style','text',...
            'String','pdatanum',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(6) = uicontrol(hf,'Style','text',...
            'String','1',...
            'Units','normalized',...
            'Position',[Pos(1)+0.25 Pos(2) AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'Visible','off');
% -- srcData
Pos = [0.3 0.65];
hImgdsply(7) = uicontrol(hf,'Style','text',...
            'String','srcData',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(8) = uicontrol(hf,'Style','text',...
            'String','''intensity''',...
            'Units','normalized',...
            'Position',[Pos(1)+0.25 Pos(2) AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'Visible','off');
% -- pgain attribute, value and slider
Pos = [0.30 0.60];
hImgdsply(9) = uicontrol(hf,'Style','text',...
            'String','pgain',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(10) = uicontrol('Style','edit','String',num2str(1.0,'%3.1f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','pgainValue',...
            'Callback',{@pgainValue_Callback},...
            'Visible','off');
hImgdsply(11) = uicontrol(hf,'Style','slider',...
            'Max',25,'Min',0.1,'Value',1.0,...
            'SliderStep',[0.1/(25-0.1), 1/(25-0.1)],...
            'Units','normalized',...
            'Position',[Pos+SG.SO,SG.SS],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','pgainSlider',...
            'Callback',{@pgainSlider_Callback},...
            'Visible','off');
% -- reject attribute, value and slider
Pos = [0.30 0.55];
hImgdsply(12) = uicontrol(hf,'Style','text',...
            'String','reject',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(13) = uicontrol('Style','edit','String',num2str(0,'%3.0f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','rejectValue',...
            'Callback',{@rejectValue_Callback},...
            'Visible','off');
hImgdsply(14) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[0.0101, 0.0505],...
            'Units','normalized',...
            'Position',[Pos+SG.SO,SG.SS],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','rejectSlider',...
            'Callback',{@rejectSlider_Callback},...
            'Visible','off');
% persist method text and popupmenu
Pos = [0.30,0.50];
hImgdsply(15) = uicontrol(hf,'Style','text',...
            'String','persistMethod',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(16) = uicontrol(hf, 'Style','popupmenu',...
            'String',{'none';'simple';'dynamic'},...
            'Units','normalized',...
            'Position',[Pos(1)+0.25 Pos(2)+0.01,0.25,0.05],...
            'FontUnits','normalized',...
            'FontSize',0.55,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','persistM',...
            'Callback',{@persistMethod_Callback},...
            'Visible','off');
% persist level text and slider(s)
Pos = [0.30,0.40];
nPLevels = 0; % no. of persist levels to display
persistLevelLength = 0; % Define here so available to sub functions
hImgdsply(17) = uicontrol(hf,'Style','text',...
            'String','persistLevel',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(18) = uicontrol('Style','edit','String',num2str(0,'%2.0f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO .062 .037],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','persistValue1',...
            'Callback',{@pLevel1Value_Callback},...
            'Visible','off');
hImgdsply(19) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[0.0101, 0.0505],...
            'Units','normalized',...
            'Position',[Pos(1)+0.30 Pos(2)-0.038 0.043 0.15],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','persistSlider1',...
            'Callback',{@pLevel1Slider_Callback},...
            'Visible','off');
hImgdsply(20) = uicontrol('Style','edit','String',num2str(0,'%2.0f'),...
            'Units','normalized',...
            'Position',[Pos(1)+0.34 Pos(2)+0.02 .062 .037],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','persistValue2',...
            'Callback',{@pLevel2Value_Callback},...
            'Visible','off');
hImgdsply(21) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[0.0101, 0.0505],...
            'Units','normalized',...
            'Position',[Pos(1)+0.39 Pos(2)-0.038 0.043 0.15],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','persistSlider2',...
            'Callback',{@pLevel2Slider_Callback},...
            'Visible','off');
hImgdsply(22) = uicontrol('Style','edit','String',num2str(0,'%2.0f'),...
            'Units','normalized',...
            'Position',[Pos(1)+0.43 Pos(2)+0.02 .062 .037],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','persistValue3',...
            'Callback',{@pLevel3Value_Callback},...
            'Visible','off');
hImgdsply(23) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[0.0101, 0.0505],...
            'Units','normalized',...
            'Position',[Pos(1)+0.48 Pos(2)-0.038 0.043 0.15],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','persistSlider3',...
            'Callback',{@pLevel3Slider_Callback},...
            'Visible','off');
hImgdsply(24) = uicontrol('Style','edit','String',num2str(0,'%2.0f'),...
            'Units','normalized',...
            'Position',[Pos(1)+0.52 Pos(2)+0.02 .062 .037],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','persistValue4',...
            'Callback',{@pLevel4Value_Callback},...
            'Visible','off');
hImgdsply(25) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[0.0101, 0.0505],...
            'Units','normalized',...
            'Position',[Pos(1)+0.57 Pos(2)-0.038 0.043 0.15],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','persistSlider4',...
            'Callback',{@pLevel4Slider_Callback},...
            'Visible','off');
% remove grain text and popupmenu
Pos = [0.30,0.30];
hImgdsply(26) = uicontrol(hf,'Style','text',...
            'String','removeGrain',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(27) = uicontrol(hf, 'Style','popupmenu',...
    'String',{'none';'low';'medium';'high'},...
    'Units','normalized',...
    'Position',[Pos(1)+0.25 Pos(2)+0.01,0.25,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','rmGrain',...
    'Callback',{@removeGrain_Callback},...
    'Visible','off');
% process method text and popupmenu
Pos = [0.30,0.25];
hImgdsply(28) = uicontrol(hf,'Style','text',...
            'String','processMethod',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(29) = uicontrol(hf, 'Style','popupmenu',...
            'String',{'none';'reduceSpeckle1';'reduceSpeckle2'},...
            'Units','normalized',...
            'Position',[Pos(1)+0.25 Pos(2)+0.01,0.34,0.05],...
            'FontUnits','normalized',...
            'FontSize',0.55,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','processM',...
            'Callback',{@processMethod_Callback},...
            'Visible','off');
% average method text and popupmenu
Pos = [0.30,0.20];
hImgdsply(30) = uicontrol(hf,'Style','text',...
            'String','averageMethod',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(31) = uicontrol(hf, 'Style','popupmenu',...
    'String',{'none';'runAverage2';'runAverage3'},...
    'Units','normalized',...
    'Position',[Pos(1)+0.25 Pos(2)+0.01,0.34,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','averageM',...
    'Callback',{@averageMethod_Callback},...
    'Visible','off');
% compress method text and popupmenu
Pos = [0.30,0.15];
hImgdsply(32) = uicontrol(hf,'Style','text',...
            'String','compressMethod',...
            'Units','normalized',...
            'Position',[Pos AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(33) = uicontrol(hf, 'Style','popupmenu',...
    'String',{'power';'log'},...
    'Units','normalized',...
    'Position',[Pos(1)+0.28 Pos(2)+0.01,0.22,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.55,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','compressM',...
    'Callback',{@compressMethod_Callback},...
    'Visible','off');
% -- compress factor attribute, value and slider
Pos = [0.30 0.10];
hImgdsply(34) = uicontrol(hf,'Style','text',...
            'String','compressFactor',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,AT.AS],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');
hImgdsply(35) = uicontrol('Style','edit','String',num2str(0,'%3.0f'),...
            'Units','normalized',...
            'Position',[Pos(1)+0.26 Pos(2)+0.02 SG.ES],...
            'BackgroundColor',bkgrnd+0.1,...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'tag','compFactorValue',...
            'Callback',{@compressFactorValue_Callback},...
            'Visible','off');
hImgdsply(36) = uicontrol(hf,'Style','slider',...
            'Max',99,'Min',0,'Value',0,...
            'SliderStep',[1/99, 5/99],...
            'Units','normalized',...
            'Position',[Pos(1)+0.35 Pos(2)+0.015 SG.SS],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'tag','compFactorSlider',...
            'Callback',{@compressFactorSlider_Callback},...
            'Visible','off');

DopplerTxt = uicontrol(hf,'Style','text',...
            'String','Doppler process is not supported in current release',...
            'Units','normalized',...
            'Position',[0.1 0.5 0.9 0.06],...
            'HorizontalAlignment','left',...
            'FontUnits','normalized',...
            'FontSize',AT.AFS,...
            'FontWeight','bold',...
            'Visible','off');

% List Process pushbutton
uicontrol(hf,'Style','pushbutton',...
            'String','List Process',...
            'Units','normalized',...
            'Position',[0.05 0.05 0.20 0.054],...
            'FontUnits','normalized',...
            'FontSize',0.5,...
            'BackgroundColor',bkgrnd+0.07,...
            'Callback',{@listProcess_Callback});


% Call update callback to check for Process structures in workspace.
procNum = 0;
Process = [];
update_Callback;

% Make the GUI visible.
set(hf,'Visible','on');

% Nested functions
    function update_Callback(~,~)
        % Look for Process structures in base workspace
        if ~evalin('base','exist(''Process'',''var'')')
            % reset Process no.; hide Class and Method text and any uicontrols
            set(hPNum,'String',{'0'});
            set(hClassImage,'Visible','off');
            set(hMethodID,'Visible','off');
            for j = 1:length(hImgdsply)
                set(hImgdsply,'Visible','off')
            end
            procNum = 0;
            return
        end
        Process = evalin('base','Process');
        nproc = size(Process,2);
        PString = cell(nproc,1);
        for j=1:nproc
            PString{j} = num2str(j,1);
        end
        set(hPNum,'String',PString);
        procNum = 1;
        processNum_Callback;
    end

    function processNum_Callback(~,~)
        if procNum == 0, return, end
        procNum = get(hPNum,'Value');
        Process = evalin('base','Process');
        set(hClassImage,'String',Process(procNum).classname);
        set(hMethodID,'String',Process(procNum).method);
        % set and make visible attributes of selected Process
        classname = Process(procNum).classname;
        switch classname
            case 'Image'

            if strcmp(Process(procNum).method,'imageDisplay')

                % Set default parameters for all attributes
                set(hImgdsply(2),'String','1'); % imgbufnum
                set(hImgdsply(4),'String','1'); % framenum
                set(hImgdsply(6),'String','1'); % pdatanum
                set(hImgdsply(8),'String','''intensity''');  % default value is intensity
                set(hImgdsply(10),'String','1.0'); % pgain edit value
                set(hImgdsply(11),'Value',1.0); % pgain slider value
                set(hImgdsply(13),'String','0'); % reject edit value
                set(hImgdsply(14),'Value',0); % reject slider value
                set(hImgdsply(16),'Value',1); % persist method
                set(hImgdsply(18),'String','0'); % persist level 1 edit value
                set(hImgdsply(19),'Value',0); % persist level 1 slider value
                set(hImgdsply(20),'String','0'); % persist level 2 edit value
                set(hImgdsply(21),'Value',0); % persist level 2 slider value
                set(hImgdsply(22),'String','0'); % persist level 3 edit value
                set(hImgdsply(23),'Value',0); % persist level 3 slider value
                set(hImgdsply(24),'String','0'); % persist level 4 edit value
                set(hImgdsply(25),'Value',0); % persist level 4 slider value
                set(hImgdsply(27),'Value',1); % remove grain
                set(hImgdsply(29),'Value',1); % process method
                set(hImgdsply(31),'Value',1); % average method
                set(hImgdsply(33),'Value',1); % compress method
                set(hImgdsply(35),'String','40'); % compress factor edit value
                set(hImgdsply(36),'Value',40); % compress factor slider value

                % Initialize compressFactor to WS
                if ~evalin('base','exist(''compFactorAll'',''var'')')
                    compFactorAll.power = ones(1,length(Process))*40;
                    compFactorAll.log = ones(1,length(Process))*40;
                    for n = 1:length(Process)
                        factorInd = find(strcmp(Process(n).Parameters,'compressFactor'));
                        if factorInd
                            powerInd = find(strcmp(Process(n).Parameters,'power'));
                            if powerInd
                                compFactorAll.power(n) = Process(n).Parameters{factorInd+1};
                            else
                                compFactorAll.log(n) = Process(n).Parameters{factorInd+1};
                            end
                        end
                    end
                    assignin('base','compFactorAll',compFactorAll);
                end

                % Overwrite defaults with attributes from Process structure

                SrcData = 'intensity';
                for k = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{k},'srcData')
                        SrcData = Process(procNum).Parameters{k+1};
                        set(hImgdsply(8),'String',SrcData);
                    end
                end

                switch SrcData
                    case {'intensity','intensity2D','intensity3D'}

                        for k = 1:2:size(Process(procNum).Parameters,2)
                            switch Process(procNum).Parameters{k}
                                case 'imgbufnum'
                                    set(hImgdsply(2),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'framenum'
                                    set(hImgdsply(4),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'pdatanum'
                                    set(hImgdsply(6),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'pgain'
                                    set(hImgdsply(10),'String',num2str(Process(procNum).Parameters{k+1},'%3.1f'));
                                    set(hImgdsply(11),'Value',Process(procNum).Parameters{k+1});
                                case 'reject'
                                    set(hImgdsply(13),'String',num2str(Process(procNum).Parameters{k+1},3));
                                    set(hImgdsply(14),'Value',Process(procNum).Parameters{k+1});
                                case 'persistMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(16),'Value',1);
                                            nPLevels = 0;
                                        case 'simple'
                                            set(hImgdsply(16),'Value',2);
                                            nPLevels = 1;
                                        case 'dynamic'
                                            set(hImgdsply(16),'Value',3);
                                            nPLevels = 4;
                                    end
                                case 'persistLevel'
                                    persistLevelLength = length(Process(procNum).Parameters{k+1});
                                    for j = 1:persistLevelLength
                                        n = 18 + 2*(j-1);
                                        set(hImgdsply(n),'String',num2str(Process(procNum).Parameters{k+1}(j),2));
                                        set(hImgdsply(n+1),'Value',Process(procNum).Parameters{k+1}(j));
                                    end
                                case 'grainRemoval'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(27),'Value',1);
                                        case 'low'
                                            set(hImgdsply(27),'Value',2);
                                        case 'medium'
                                            set(hImgdsply(27),'Value',3);
                                        case 'high'
                                            set(hImgdsply(27),'Value',4);
                                    end
                                case 'processMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(29),'Value',1);
                                        case 'reduceSpeckle1'
                                            set(hImgdsply(29),'Value',2);
                                        case 'reduceSpeckle2'
                                            set(hImgdsply(29),'Value',3);
                                    end
                                case 'averageMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(31),'Value',1);
                                        case 'runAverage2'
                                            set(hImgdsply(31),'Value',2);
                                        case 'runAverage3'
                                            set(hImgdsply(31),'Value',3);
                                    end
                                case 'compressMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'power'
                                            set(hImgdsply(33),'Value',1);
                                        case 'log'
                                            set(hImgdsply(33),'Value',2);
                                    end
                                case 'compressFactor'
                                    compFactorAll = evalin('base','compFactorAll');
                                    if isequal(get(hImgdsply(33),'Value'),1)
                                        compFactorAll.power(procNum) = Process(procNum).Parameters{k+1};
                                    else
                                        compFactorAll.log(procNum) = Process(procNum).Parameters{k+1};
                                    end
                                    set(hImgdsply(35),'String',num2str(Process(procNum).Parameters{k+1},2));
                                    set(hImgdsply(36),'Value',Process(procNum).Parameters{k+1});
                                    assignin('base','compFactorAll',compFactorAll);
                                otherwise
                            end
                        end
                        % Make all controls visible.
                        for k = 1:17
                            set(hImgdsply(k),'Visible','on');
                        end
                        for k = 1:4  % turn off all persist levels controls
                            set(hImgdsply(18+2*(k-1)),'Visible','off');
                            set(hImgdsply(18+2*(k-1)+1),'Visible','off');
                        end

                        % turn on persist levels that are active
                        for k = 1:nPLevels
                            set(hImgdsply(18+2*(k-1)),'Visible','on');
                            set(hImgdsply(18+2*(k-1)+1),'Visible','on');
                        end
                        for k = 26:36
                            set(hImgdsply(k),'Visible','on');
                        end

                    case 'unsignedColor'

                        % turn off all sliders first
                        set(hImgdsply,'Visible','off');

                        for k = 1:2:size(Process(procNum).Parameters,2)
                            switch Process(procNum).Parameters{k}
                                case 'imgbufnum'
                                    set(hImgdsply(2),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'framenum'
                                    set(hImgdsply(4),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'pdatanum'
                                    set(hImgdsply(6),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'persistMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(16),'Value',1);
                                            nPLevels = 0;
                                        case 'simple'
                                            set(hImgdsply(16),'Value',2);
                                            nPLevels = 1;
                                        case 'dynamic'
                                            set(hImgdsply(16),'Value',3);
                                            nPLevels = 4;
                                    end
                                case 'persistLevel'
                                    persistLevelLength = length(Process(procNum).Parameters{k+1});
                                    for j = 1:persistLevelLength
                                        n = 18 + 2*(j-1);
                                        set(hImgdsply(n),'String',num2str(Process(procNum).Parameters{k+1}(j),2));
                                        set(hImgdsply(n+1),'Value',Process(procNum).Parameters{k+1}(j));
                                    end
                                case 'grainRemoval'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(27),'Value',1);
                                        case 'low'
                                            set(hImgdsply(27),'Value',2);
                                        case 'medium'
                                            set(hImgdsply(27),'Value',3);
                                        case 'high'
                                            set(hImgdsply(27),'Value',4);
                                    end
                                otherwise
                            end
                        end

                        % Make all controls visible.
                        for k = [(1:8),(15:17),26,27]
                            set(hImgdsply(k),'Visible','on');
                        end

                        if nPLevels > 1, nPLevels = persistLevelLength; end
                        % turn on persist levels that are active
                        for k = 1:nPLevels
                            set(hImgdsply(18+2*(k-1)),'Visible','on');
                            set(hImgdsply(18+2*(k-1)+1),'Visible','on');
                        end

                    case 'signedColor'
                        % turn off all sliders first
                        set(hImgdsply,'Visible','off');

                        for k = 1:2:size(Process(procNum).Parameters,2)
                            switch Process(procNum).Parameters{k}
                                case 'imgbufnum'
                                    set(hImgdsply(2),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'framenum'
                                    set(hImgdsply(4),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'pdatanum'
                                    set(hImgdsply(6),'String',num2str(Process(procNum).Parameters{k+1},2));
                                case 'persistMethod'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(16),'Value',1);
                                            nPLevels = 0;
                                        case 'simple'
                                            set(hImgdsply(16),'Value',2);
                                            nPLevels = 1;
                                        case 'dynamic'
                                            set(hImgdsply(16),'Value',3);
                                            nPLevels = 4;
                                    end
                                case 'persistLevel'
                                    persistLevelLength = length(Process(procNum).Parameters{k+1});
                                    for j = 1:persistLevelLength
                                        n = 18 + 2*(j-1);
                                        set(hImgdsply(n),'String',num2str(Process(procNum).Parameters{k+1}(j),2));
                                        set(hImgdsply(n+1),'Value',Process(procNum).Parameters{k+1}(j));
                                    end
                                case 'grainRemoval'
                                    Strng = Process(procNum).Parameters{k+1};
                                    switch Strng
                                        case 'none'
                                            set(hImgdsply(27),'Value',1);
                                        case 'low'
                                            set(hImgdsply(27),'Value',2);
                                        case 'medium'
                                            set(hImgdsply(27),'Value',3);
                                        case 'high'
                                            set(hImgdsply(27),'Value',4);
                                    end
                                otherwise
                            end
                        end
                        % Make all controls visible.
                        for k = [(1:8),(15:17),26,27]
                            set(hImgdsply(k),'Visible','on');
                        end

                        % turn on persist levels that are active
                        for k = 1:nPLevels
                            set(hImgdsply(18+2*(k-1)),'Visible','on');
                            set(hImgdsply(18+2*(k-1)+1),'Visible','on');
                        end
                end

                set(DopplerTxt,'Visible','off');
            end

            case 'Doppler'
                % Make all controls invisible.
                set(hImgdsply,'Visible','off');
                set(DopplerTxt,'Visible','on');

            case 'External'
                % Make all controls invisible.
                set(hImgdsply,'Visible','off');

        end
    end

% PGain Slider Callback
    function pgainSlider_Callback(src,~)
        pv = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'pgain' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'pgain')
                    Process(procNum).Parameters{j+1} = pv;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'pgain' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'pgain';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = pv;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                                   'pgain',pv};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);

        end
        set(findobj('tag','pgainValue'),'String',num2str(pv,'%2.1f'));
    end

% pgainValue Callback
    function pgainValue_Callback(src,~)
        pv = str2double(get(src,'String'));
        if (0.1<=pv) && (pv<=16.0)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'pgain' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'pgain')
                        Process(procNum).Parameters{j+1} = pv;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'pgain' attribute to Process structure
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'pgain';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = pv;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                                       'pgain',pv};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','pgainSlider'),'Value',pv);
            end
        else
            pv = get(findobj('tag','pgainSlider'),'Value');
            set(findobj('tag','pgainValue'),'String',num2str(pv,'%2.1f'));
        end
    end

    % reject Slider Callback
    function rejectSlider_Callback(src,~)
        rej = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'reject' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'reject')
                    Process(procNum).Parameters{j+1} = rej;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'reject' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'reject';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = rej;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                                   'reject',rej};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','rejectValue'),'String',num2str(rej,'%2.0f'));
    end

% reject Value Callback
    function rejectValue_Callback(src,~)
        rej = str2double(get(src,'String'));
        if (rej>=0) && (rej<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'reject' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'reject')
                        Process(procNum).Parameters{j+1} = rej;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'reject' attribute to Process structure
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'reject';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = rej;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'reject',rej};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','rejectSlider'),'Value',rej);
            end
        else
            rej = get(findobj('tag','rejectSlider'),'Value');
            set(findoj('tag','rejectValue'),'String',num2str(rej,'%2d'));
        end
    end

% persist method callback
    function persistMethod_Callback(src,~)
        persistVal = get(src,'Value');
        switch persistVal
            case 1
                persistStrng = 'none';
                for k = 1:4
                    set(hImgdsply(18+2*(k-1)),'Visible','off');
                    set(hImgdsply(18+2*(k-1)+1),'Visible','off');
                end
            case 2
                persistStrng = 'simple';
                set(hImgdsply(18),'Visible','on');
                set(hImgdsply(19),'Visible','on');
                for k = 2:4
                    set(hImgdsply(18+2*(k-1)),'Visible','off');
                    set(hImgdsply(18+2*(k-1)+1),'Visible','off');
                end
            case 3
                persistStrng = 'dynamic';
                if persistLevelLength > 1
                    for k = 1:persistLevelLength
                        set(hImgdsply(18+2*(k-1)),'Visible','on');
                        set(hImgdsply(18+2*(k-1)+1),'Visible','on');
                    end
                else % turn on all sliders
                    for k = 1:4
                        set(hImgdsply(18+2*(k-1)),'Visible','on');
                        set(hImgdsply(18+2*(k-1)+1),'Visible','on');
                    end
                end
        end
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'processMethod' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'persistMethod')
                    Process(procNum).Parameters{j+1} = persistStrng;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'processMethod' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistMethod';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = persistStrng;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                                   'persistMethod',persistStrng};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
    end

% persist Slider 1 Callback
    function pLevel1Slider_Callback(src,~)
        pl = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'persistLevel' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'persistLevel')
                    PL = Process(procNum).Parameters{j+1};
                    PL(1) = pl;
                    Process(procNum).Parameters{j+1} = PL;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'persistLevel' attribute to Process structure
                PL = pl;
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'persistLevel',PL};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','persistValue1'),'String',num2str(pl,'%2.0f'));
        % modify persistence in Doppler script
        if isequal(procNum,3)
            h = findobj('String','Color Persistence');
            if ishandle(h)
                UI = evalin('base','UI');
                for i = 1:length(UI)
                    if find(h==UI(i).handle)
                        set(UI(i).handle(2),'Value',pl);
                        set(UI(i).handle(3),'String',num2str(pl,'%3.0f'));
                    end
                end
            end
        end
    end

% persist Value 1 Callback
    function pLevel1Value_Callback(src,~)
        pl = str2double(get(src,'String'));
        if (pl>=0) && (pl<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'persistLevel' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'persistLevel')
                        PL = Process(procNum).Parameters{j+1};
                        PL(1) = pl;
                        Process(procNum).Parameters{j+1} = PL;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'persistLevel' attribute to Process structure
                    PL = pl;
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'persistLevel',PL};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','persistSlider1'),'Value',pl);
            end
        else
            pl = get(findobj('tag','persistSlider1'),'Value');
            set(findoj('tag','persistValue1'),'String',num2str(pl,'%2.0f'));
        end
        % modify persistence in Doppler script
        if isequal(procNum,3)
            h = findobj('String','Color Persistence');
            if ishandle(h)
                UI = evalin('base','UI');
                for i = 1:length(UI)
                    if find(h==UI(i).handle)
                        set(UI(i).handle(2),'Value',pl);
                        set(UI(i).handle(3),'String',num2str(pl,'%3.0f'));
                    end
                end
            end
        end
    end

% persist Slider 2 Callback
    function pLevel2Slider_Callback(src,~)
        pl = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'persistLevel' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'persistLevel')
                    PL = Process(procNum).Parameters{j+1};
                    PL(2) = pl;
                    Process(procNum).Parameters{j+1} = PL;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'persistLevel' attribute to Process structure
                PL = [0 pl 0 0];
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'persistLevel',PL};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','persistValue2'),'String',num2str(pl,'%2.0f'));
    end

% persist Value 2 Callback
    function pLevel2Value_Callback(src,~)
        pl = str2double(get(src,'String'));
        if (pl>=0) && (pl<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'persistLevel' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'persistLevel')
                        PL = Process(procNum).Parameters{j+1};
                        PL(2) = pl;
                        Process(procNum).Parameters{j+1} = PL;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'persistLevel' attribute to Process structure
                    PL = [0 pl 0 0];
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'persistLevel',PL};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','persistSlider2'),'Value',pl);
            end
        else
            pl = get(findobj('tag','persistSlider2'),'Value');
            set(findoj('tag','persistValue2'),'String',num2str(pl,'%2.0f'));
        end
    end

% persist Slider 3 Callback
    function pLevel3Slider_Callback(src,~)
        pl = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'persistLevel' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'persistLevel')
                    PL = Process(procNum).Parameters{j+1};
                    PL(3) = pl;
                    Process(procNum).Parameters{j+1} = PL;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'persistLevel' attribute to Process structure
                PL = [0 0 pl 0];
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'persistLevel',PL};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','persistValue3'),'String',num2str(pl,'%2.0f'));
    end

% persist Value 3 Callback
    function pLevel3Value_Callback(src,~)
        pl = str2double(get(src,'String'));
        if (pl>=0) && (pl<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'persistLevel' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'persistLevel')
                        PL = Process(procNum).Parameters{j+1};
                        PL(3) = pl;
                        Process(procNum).Parameters{j+1} = PL;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'persistLevel' attribute to Process structure
                    PL = [0 0 pl 0];
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'persistLevel',PL};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','persistSlider3'),'Value',pl);
            end
        else
            pl = get(findobj('tag','persistSlider3'),'Value');
            set(findoj('tag','persistValue3'),'String',num2str(pl,'%2.0f'));
        end
    end

% persist Slider 4 Callback
    function pLevel4Slider_Callback(src,~)
        pl = get(src,'Value');
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'persistLevel' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'persistLevel')
                    PL = Process(procNum).Parameters{j+1};
                    PL(4) = pl;
                    Process(procNum).Parameters{j+1} = PL;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'persistLevel' attribute to Process structure
                PL = [0 0 0 pl];
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'persistLevel',PL};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','persistValue4'),'String',num2str(pl,'%2.0f'));
    end

% persist Value 4 Callback
    function pLevel4Value_Callback(src,~)
        pl = str2double(get(src,'String'));
        if (pl>=0) && (pl<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'persistLevel' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'persistLevel')
                        PL = Process(procNum).Parameters{j+1};
                        PL(4) = pl;
                        Process(procNum).Parameters{j+1} = PL;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'persistLevel' attribute to Process structure
                    PL = [0 0 0 pl];
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'persistLevel';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = PL;
                end
                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'persistLevel',PL};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','persistSlider4'),'Value',pl);
            end
        else
            pl = get(findobj('tag','persistSlider4'),'Value');
            set(findoj('tag','persistValue4'),'String',num2str(pl,'%2.0f'));
        end
    end

    function removeGrain_Callback(src,~)
        rmGVal = get(src,'Value');
        switch rmGVal
            case 1
                rmGStrng = 'none';
            case 2
                rmGStrng = 'low';
            case 3
                rmGStrng = 'medium';
            case 4
                rmGStrng = 'high';
        end
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'grainRemoval' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'grainRemoval')
                    Process(procNum).Parameters{j+1} = rmGStrng;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'grainRemoval' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'grainRemoval';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = rmGStrng;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'grainRemoval',rmGStrng};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
    end

    function processMethod_Callback(src,~)
        procMVal = get(src,'Value');
        switch procMVal
            case 1
                procMStrng = 'none';
            case 2
                procMStrng = 'reduceSpeckle1';
            case 3
                procMStrng = 'reduceSpeckle2';
        end
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'processMethod' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'processMethod')
                    Process(procNum).Parameters{j+1} = procMStrng;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'processMethod' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'processMethod';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = procMStrng;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'processMethod',procMStrng};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
    end

    function averageMethod_Callback(src,~)
        avgVal = get(src,'Value');
        switch avgVal
            case 1
                avgStrng = 'none';
            case 2
                avgStrng = 'runAverage2';
            case 3
                avgStrng = 'runAverage3';
        end
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'averageMethod' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'averageMethod')
                    Process(procNum).Parameters{j+1} = avgStrng;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'averageMethod' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'averageMethod';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = avgStrng;
            end
            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'averageMethod',avgStrng};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
    end

    function compressMethod_Callback(src,~)
        compVal = get(src,'Value');
        switch compVal
            case 1
                compStrng = 'power';
            case 2
                compStrng = 'log';
        end
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'compressMethod' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'compressMethod')
                    Process(procNum).Parameters{j+1} = compStrng;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'compressMethod' attribute to Process structure
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'compressMethod';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = compStrng;
            end

            compFactorAll = evalin('base','compFactorAll');
            if isequal(compVal,1)
                compFactor = compFactorAll.power(procNum);
            else
                compFactor = compFactorAll.log(procNum);
            end
            set(hImgdsply(35),'String',num2str(compFactor));
            set(hImgdsply(36),'Value',compFactor);

            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'compressMethod',compStrng};
            Cntrl(n+1).Command = 'set&Run';
            Cntrl(n+1).Parameters = {'Process',procNum,...
                'compressFactor',compFactor};

            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
    end

    function compressFactorSlider_Callback(src,~)
        compFactor = round(get(src,'Value'));
        if exist('Process','var')&&~isempty(Process)
            flag = 0; % flag for 'compressFactor' set in Process
            for j = 1:2:size(Process(procNum).Parameters,2)
                if strcmp(Process(procNum).Parameters{j},'compressFactor')
                    Process(procNum).Parameters{j+1} = compFactor;
                    flag = 1;
                    break
                end
            end
            if flag==0 % add 'compressFactor' attribute to Process structure and compressFactor variable to workspace
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'compressFactor';
                Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = compFactor;
            end

            compFactorAll = evalin('base','compFactorAll');
            if isequal(get(hImgdsply(33),'Value'),1)
                compFactorAll.power(procNum) = compFactor;
            else
                compFactorAll.log(procNum) = compFactor;
            end
            assignin('base','compFactorAll',compFactorAll);

            if evalin('base','exist(''Control'',''var'')')
                Cntrl = evalin('base','Control');
            else
                Cntrl = struct('Command',[],'Parameters',[]);
            end
            if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
            Cntrl(n).Command = 'set&Run';
            Cntrl(n).Parameters = {'Process',procNum,...
                'compressFactor',compFactor};
            assignin('base','Process',Process);
            assignin('base','Control',Cntrl);
        end
        set(findobj('tag','compFactorValue'),'String',num2str(compFactor,'%2.0f'));
    end

    function compressFactorValue_Callback(src,~)
        compFactor = str2double(get(src,'String'));
        if (compFactor>=0) && (compFactor<=99)
            if exist('Process','var')&&~isempty(Process)
                flag = 0; % flag for 'compressFactor' set in Process
                for j = 1:2:size(Process(procNum).Parameters,2)
                    if strcmp(Process(procNum).Parameters{j},'compressFactor')
                        Process(procNum).Parameters{j+1} = compFactor;
                        flag = 1;
                        break
                    end
                end
                if flag==0 % add 'compressFactor' attribute to Process structure and compressFactor variable to workspace
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = 'compressFactor';
                    Process(procNum).Parameters{size(Process(procNum).Parameters,2)+1} = compFactor;
                end

                compFactorAll = evalin('base','compFactorAll');
                if isequal(get(hImgdsply(33),'Value'),1)
                    compFactorAll.power(procNum) = compFactor;
                else
                    compFactorAll.log(procNum) = compFactor;
                end
                assignin('base','compFactorAll',compFactorAll);

                if evalin('base','exist(''Control'',''var'')')
                    Cntrl = evalin('base','Control');
                else
                    Cntrl = struct('Command',[],'Parameters',[]);
                end
                if isempty(Cntrl(1).Command), n=1; else n=length(Cntrl)+1; end
                Cntrl(n).Command = 'set&Run';
                Cntrl(n).Parameters = {'Process',procNum,...
                    'compressFactor',compFactor};
                assignin('base','Process',Process);
                assignin('base','Control',Cntrl);
                set(findobj('tag','compFactorSlider'),'Value',compFactor);
            end
        else
            compFactor = get(findobj('tag','compFactorSlider'),'Value');
            set(findoj('tag','compFactorValue'),'String',num2str(compFactor,'%2d'));
        end
    end
    function closefunc(~,~)
        delete(hf);
        set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
    end

    function listProcess_Callback(~,~)
        % List the Process Object on the command line.
        fprintf('Process(%u).classname = ''%s'';\n',procNum,Process(procNum).classname);
        fprintf('Process(%u).method = ''%s'';\n',procNum,Process(procNum).method);
        fprintf('Process(%u).Parameters = {',procNum);
        n = size(Process(procNum).Parameters,2);
        for j=1:2:n
            if j~=1, fprintf('                         '); end
            fprintf('''%s'',',Process(procNum).Parameters{j});
            switch class(Process(procNum).Parameters{j+1})
                case 'double'
                    v = Process(procNum).Parameters{j+1};
                    if length(v) == 1
                        if (v-fix(v))<0.001, fprintf('%2d',fix(v));
                        else fprintf('%2.1f',v);
                        end
                    else
                        fprintf('[%2.1f',v(1));
                        for k = 2:length(v)
                            fprintf(' %2.1f',v(k));
                        end
                        fprintf(']');
                    end
                case 'char'
                    fprintf('''%s''',Process(procNum).Parameters{j+1});
            end
            if j<(n-1), fprintf(',...\n');
            else fprintf('};\n');
            end
        end
    end


end

