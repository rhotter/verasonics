function vsx_gui
%
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% VSX_GUI  This gui is opened by the main VSX program and allows control of various
% acquisition and processing parameters.

% Revised
%   April 23, 2020  VTS-1695 Delete pan and zoom controls
%   July 17, 2019 for 4.1.0 release

% Initialize variables that will be available in nested scopes.

% Variables that do not change state after VSX initialization.
trackP5 = 0;     % profile 1 tracks profile 5 when true.
hvps2 = 0;       % hv profile to use for 2nd slider.
nfrms = 0;       % number of DisplayWindow frames.
vpfRoot = '';    % root directory of Vantage installation.
vsTempDir = '';  % location of temporary directory.

% State variables for TGC controls.
TGCAll = [];     % all slider position for each TGC profile.
SPAll = [];      % all slider gain factor for each TGC profile.
SP = [];         % active TGC slider positions.
tgcSldrs   = [];   % tgcSliders
tgcAllSldr = [];
tgcSelect  = [];

hv1Sldr    = [];
hv1Value   = [];

hv2Sldr    = [];
hv2Value   = [];

speedSldr  = [];
speedValue = [];

initialize()

    function initialize()
        % Are we running with hardware? Get state of VDAS variable.
        VDAS = evalin('base','VDAS');
        Resource = evalin('base', 'Resource');
        vpfRoot = Resource.SysConfig.vpfRoot;
        vsTempDir = evalin('base','vsTempDir');

        % Close any previously opened GUI windows.
        delete(findobj('tag','UI'));
        % Initialize and hide the GUI as it is being constructed.
        ScrnSize = get(0,'ScreenSize');
        Bkgrnd = [0.8 0.8 0.8];
        f = figure('Visible','off',...  %'Units','normalized',...
            'Position',[ScrnSize(3)-500,(ScrnSize(4)-620)/2,450,620],... %'Position',[0.7,0.25,0.25,0.50],...
            'Name','VSX Control',...
            'Color',Bkgrnd,...
            'NumberTitle','off',...
            'MenuBar','none', ...
            'Resize','on', ...
            'DefaultUicontrolBackgroundColor',Bkgrnd,...
            'tag','UI');

        % Determine number of DisplayWindow frames.
        dsplywin = 0;
        if isfield(Resource,'DisplayWindow')
            dsplywin = 1;
            if isfield(Resource.DisplayWindow,'numFrames')
                nfrms = Resource.DisplayWindow(1).numFrames;
            else
                nfrms = 1;
            end
        end

        % ***** Create the GUI components *****
        % Define UIPos, which contains the default GUI positions - three columns of 10 controls. The x,y
        %    locations increment up columns, with each column being a separate page. The origin
        %    specified by UIPos is the lower left corner of a virtual box that encloses the control.
        UIPos = zeros(10,2,3);
        UIPos(:,1,1) = 0.0625;
        UIPos(:,1,2) = 0.375;
        UIPos(:,1,3) = 0.6875;
        UIPos(:,2,1) = 0.0:0.1:0.9;
        UIPos(:,2,2) = 0.0:0.1:0.9;
        UIPos(:,2,3) = 0.0:0.1:0.9;
        assignin('base','UIPos',UIPos);
        % Define slider group offsets and sizes. All units are normalized.
        SG = struct('TO',[0.0,0.0975],...   % title offset
            'TS',[0.25,0.025],...   % title size
            'TF',0.8,...            % title font size
            'SO',[0.0,0.06],...     % slider offset
            'SS',[0.25,0.031],...   % slider size
            'EO',[0.075,0.031],...  % edit box offset
            'ES',[0.11,0.031]);     % edit box size
        assignin('base','SG',SG);
        % Define pushbutton offsets and sizes.
        PB = struct('FS',0.3,...    % font size
            'BO',[0.025,0.04],...   % button offset
            'BS',[0.2,0.07]);       % button size
        assignin('base','PB',PB);
        % Define buttonGroup offsets and sizes for a single button.
        BG = struct('TFS',0.25,...  % button group title font size
            'BGO',[0,0.025],...     % button group offset (BGO(2) is 0.1 - BGS(2))
            'BGS',[0.25,0.075],...  % button group size
            'BI',0.030,...          % button increment (in units of full window)
            'BO',[0.1,0.2],...      % button offset (units relative to BGS box)
            'BS',[0.9,0.4],...      % button size      "    "    "    "     "
            'BFS',0.8);             % button font size
        assignin('base','BG',BG);

        % Titles
        Pos = UIPos(10,:,1);
        uicontrol('Style','text','String',...
            'Front End',...
            'Units','normalized',...
            'Position',[Pos+[0.0 0.06],0.25,0.03],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','title1',...
            'FontWeight','bold');
        uicontrol('Style','text','String',...
            'Processing',...
            'Units','normalized',...
            'Position',[Pos+[0.3175 0.06],0.25,0.03],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','title2',...
            'FontWeight','bold');
        uicontrol('Style','text','String',...
            'Display',...
            'Units','normalized',...
            'Position',[Pos+[0.625 0.06],0.25,0.03],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','title3',...
            'FontWeight','bold');

        % Transducer connector
        utaType = Resource.SysConfig.UTAtype;
        if ~isempty(utaType) && utaType(3) > 1
            tcx = 0.09;
            tcy = 0.94;
            tcw = 0.2;
            Conn = Resource.Parameters.Connector;
            conStr = strjoin(string(Conn), ', ');
            if length(Conn) == 1
                TcString = "Using Connector " + conStr;
            else
                TcString = "Using Connectors " + conStr;
                tcx = tcx - 0.005*length(Conn);
                tcw = tcw + 0.01*length(Conn);
            end
            uicontrol('Style','text',...
                'String',TcString,...
                'Units','normalized',...
                'Position',[tcx tcy tcw 0.018],...
                'FontUnits','normalized',...
                'FontSize',0.8,...
                'BackgroundColor',[1.0,0.8,0.4]);
        end

        % TGC Controls
        Pos = UIPos(9,:,1);
        tgcTxt = uicontrol('Style','text','String','TGC',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,SG.TS],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','tgctxt',...
            'FontWeight','bold');

        % Set sliders to TGC.CntlPts positions, if TGC exists.
        if evalin('base','exist(''TGC'',''var'')')
            TGC = evalin('base','TGC');
            SP = TGC(1).CntrlPts/1023;
            numTGCs = length(TGC);
        else
            SP = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
            numTGCs = 1;
        end

        % if more than 1 TGC, have a dropdown menu for TGC numbers
        tgcSelect = gobjects(1); % Initialize with empty graphics object
        if numTGCs > 1
            TGCAll = zeros(numTGCs, 1); % TGCAll is used to save the slider value
            SPAll = ones(numTGCs, 1);
            set(tgcTxt,'Position',[Pos+SG.TO-[0.05,-0.003],SG.TS]);
            tgcSelect = uicontrol('Style','popupmenu',...
                'Units','normalized',...
                'Position',[Pos+SG.TO+[0.125,0.01],0.125,0.02],...
                'String',arrayfun(@num2str, 1:numTGCs, 'UniformOutput', false),...
                'FontUnits','normalized',...
                'FontSize',0.9,...
                'Tag','TGCnum');
        else
            TGCAll = 0; % TGCAll is used to save the slider value
            SPAll = 1;
        end

        tgcSldrInc = 0.0325; % increment between TGC sliders
        Pos = Pos+SG.SO;
        tgcSldrs = gobjects(length(SP), 1); % Initialize with empty graphics object
        for i=1:length(SP)
            tgcSldrs(i) = uicontrol(f,'Style','slider',...
                'Max',1.0,'Min',0,'Value',SP(i),...
                'SliderStep',[0.05 0.2],...
                'Units','normalized',...
                'Position',[Pos+[0 -(i-1)*tgcSldrInc],SG.SS],...
                'BackgroundColor',Bkgrnd-0.05,...
                'Tag',sprintf('tgc%dSldr',i),...
                'Callback',{@tgcSldr_Callback, i});
        end

        Pos = UIPos(6,:,1);
        uicontrol(f,'Style','text','String','TGC All Gain',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,SG.TS],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','tgcAllTxt',...
            'FontWeight','bold');
        tgcAllSldr = uicontrol(f,'Style','slider',...
            'Max',1.0,'Min',-1.0,'Value',0,...
            'SliderStep',[0.01,0.04],...
            'Units','normalized',...
            'Position',[Pos+SG.SO,SG.SS],...
            'Tag','tgcAllSldr',...
            'BackgroundColor',Bkgrnd-0.05,...
            'Callback',{@tgcAllSldr_Callback, tgcSldrs});
        if isgraphics(tgcSelect)
            tgcSelect.Callback = {@tgcSelect_Callback,tgcSldrs,tgcAllSldr};
        end

        % High Voltage slider.
        %   The slider's min and max range is determined by the limits of the Verasonics
        %   TPC and a max voltage limit that can be set in the user's setup script.
        %
        %   The user script may override the "Max" attribute of this slider to impose
        %   a high voltage maximum limit that is less than the capability of the detected
        %   Verasonics TPC.
        TPC = evalin('base','TPC');

        % - For more than one active profile, determine the number to use for the 2nd slider. Profile 5
        %   has priority over profiles 2-4.
        if TPC(5).inUse
            hvps2 = 5;
        else
            for i = 2:4
                if TPC(i).inUse ~= 0
                    hvps2 = i;
                    break,
                end
            end
        end
        assignin('base', 'hv2GUIprofile', hvps2);

        rangeMin = evalin('base', 'minTpcVoltage');
        hv1Sldr = gobjects(1); % Initialize with empty graphics object
        hv1Value = gobjects(1); % Initialize with empty graphics object
        if TPC(1).inUse
            % - Render HV control number 1
            Pos = UIPos(4,:,1);
            uicontrol('Style','text','String','High Voltage P1',...
                'Units','normalized',...
                'Position',[Pos+SG.TO,SG.TS],...
                'FontUnits','normalized',...
                'FontSize',0.8,...
                'Tag','hv1txt',...
                'FontWeight','bold');
            hv1Sldr = uicontrol(f,'Style','slider',...
                'Max',TPC(1).highVoltageLimit,...
                'Min',rangeMin,...
                'Value',TPC(1).hv,...
                'SliderStep',[0.05 0.2],...
                'Units','normalized',...
                'Position',[Pos+SG.SO,SG.SS],...
                'BackgroundColor',Bkgrnd-0.05,...
                'Interruptible','off',...
                'BusyAction','cancel',...
                'Tag','hv1Sldr');
            hv1Value = uicontrol('Style','edit','String',num2str(TPC(1).hv,'%.1f'),...
                'Units','normalized',...
                'Position',[Pos+SG.EO,SG.ES],...
                'Tag','hv1Value',...
                'Callback',{@hv1Value_Callback,hv1Sldr},...
                'BackgroundColor',Bkgrnd+0.1);
            hv1Sldr.Callback = {@hv1Sldr_Callback, hv1Value};
        end
        % - If more than one profile, create 2nd high voltage control
        hvTimer = gobjects(1); % Initialize with empty graphics object
        if hvps2 ~= 0
            Pos = UIPos(3,:,1);
            uicontrol('Style','text',...
                'String',['High Voltage P' num2str(hvps2)],...
                'Units','normalized',...
                'Position',[Pos+SG.TO,SG.TS],...
                'FontUnits','normalized',...
                'FontSize',0.8,...
                'Tag','hv2txt',...
                'FontWeight','bold');
            hv2Sldr = uicontrol(f,'Style','slider',...
                'Max',TPC(hvps2).highVoltageLimit,...
                'Min',rangeMin,...
                'Value',TPC(hvps2).hv,...
                'SliderStep',[0.05 0.2],...
                'Units','normalized',...
                'Position',[Pos+SG.SO,SG.SS],...
                'BackgroundColor',Bkgrnd-0.05,...
                'Interruptible','off',...
                'BusyAction','cancel',...
                'Tag','hv2Sldr');
            hv2Value = uicontrol('Style','edit','String',num2str(TPC(hvps2).hv,'%.1f'),...
                'Units','normalized',...
                'Position',[Pos+SG.EO,SG.ES],...
                'Tag','hv2Value',...
                'Callback',{@hv2Value_Callback,hv2Sldr,hv1Sldr,hv1Value},...
                'BackgroundColor',Bkgrnd+0.1);
            hv2Sldr.Callback = {@hv2Sldr_Callback, hv2Value, hv1Sldr, hv1Value};
            if hvps2 == 5
                % check for profile voltage tracking parameter
                % value has already been initialized to default of 'off'
                if isfield(Resource,'HIFU') && ...
                        isfield(Resource.HIFU,'voltageTrackP5') && ...
                        ~isempty(Resource.HIFU.voltageTrackP5)
                    trackP5 = Resource.HIFU.voltageTrackP5;
                end
                hv2Value.Position = [Pos+[0.045 0.03],SG.ES];
                hv2Actual = uicontrol('Style','edit','String',num2str(TPC(hvps2).hv,'%.1f'),...
                    'Units','normalized',...
                    'Position',[Pos+[0.125 0.03],SG.ES],...
                    'Enable','inactive',...
                    'Tag','hv2Actual',...
                    'BackgroundColor',Bkgrnd+0.05);
                if VDAS
                    % Create a timer object that will update the actual push capacitor voltage string.
                    hvTimer = timer('TimerFcn',{@hvTimer_Callback,hv2Sldr,hv2Actual},...
                        'Period', 1.5,'ExecutionMode','fixedSpacing');
                    start(hvTimer);
                end
            end
        end

        % - RcvData loop control
        uicontrol('Style','togglebutton',...
            'String','Rcv Data Loop',...
            'Units','normalized',...
            'Position',[UIPos(1,:,1)+[0.0175 0.077],0.20,0.05],...
            'FontUnits','normalized',...
            'FontSize',0.4,...
            'Tag','RcvLoop',...
            'BackgroundColor',Bkgrnd+0.05,...
            'Callback',{@rcvdataloop_Callback});

        % - Simulate control
        if VDAS
            uicontrol('Style','togglebutton',...
                'String','Simulate',...
                'Units','normalized',...
                'Position',[UIPos(1,:,1)+[0.0175 0.027],0.20,0.05],...
                'FontUnits','normalized',...
                'FontSize',0.4,...
                'Tag','Simulation',...
                'BackgroundColor',Bkgrnd+0.05,...
                'Callback',{@simulate_Callback});
        end
        % - Add Speed Correction slider only if a Recon is defined.
        if evalin('base','exist(''Recon'',''var'')')
            Pos = UIPos(9,:,2);
            uicontrol('Style','text','String','Speed Of Sound',...
                'Units','normalized',...
                'Position',[Pos+SG.TO,SG.TS],...
                'FontUnits','normalized',...
                'FontSize',0.8,...
                'FontWeight','bold',...
                'Tag','speedtxt');
            speedSldr = uicontrol(f,'Style','slider',...
                'Max',1.4,'Min',0.6,'Value',1.0,...
                'SliderStep',[0.00125 0.0125],...
                'Units','normalized',...
                'Position',[Pos+SG.SO,SG.SS],...
                'BackgroundColor',Bkgrnd-0.05,...
                'Interruptible','off',...
                'BusyAction','cancel',...
                'Tag','speedSldr');
            speedValue = uicontrol('Style','edit','String',num2str(1.0,'%1.3f'),...
                'Units','normalized',...
                'Position',[Pos+SG.EO,SG.ES],...
                'Tag','speedValue',...
                'Callback',{@speedValue_Callback,speedSldr},...
                'BackgroundColor',Bkgrnd+0.1);
            speedSldr.Callback = {@speedSldr_Callback,speedValue};
        end

        % - Freeze control
        freezeBtn = uicontrol('Style','togglebutton',...
            'String','Freeze',...
            'Units','normalized',...
            'Position',[UIPos(1,:,2)+PB.BO,PB.BS],...
            'FontUnits','normalized',...
            'FontSize',PB.FS,...
            'BackgroundColor',Bkgrnd+0.05,...
            'Tag','FreezeButton',...
            'Callback',{@freeze_Callback});

        toolStr = {'none';'filterTool';'showTXPD';'saveRF'};

        % Add the following GUI controls only if we have a DisplayWindow specification.
        if dsplywin ~= 0
            if evalin('base','exist(''Process'',''var'')')
                Process = evalin('base','Process');
                % Find values used in the first 'Image/imageDisplay' Process structure for displayWindow 1
                for i = 1:size(Process,2)
                    if strcmp(Process(i).classname,'Image') && strcmp(Process(i).method,'imageDisplay')
                        toolStr = {'none';'filterTool';'showTXPD';'saveRF';'PTool';'ColorMapTool'};
                        break
                    end
                end
            end
        end

        % Tools popup menu
        Pos = UIPos(7,:,3);
        uicontrol('Style','text',...
            'String','Tools',...
            'Units','normalized',...
            'Position',[Pos+[0.015 0.05],0.22,0.06],...
            'HorizontalAlignment','Center',...
            'FontUnits','normalized',...
            'FontSize',0.35,...
            'FontWeight','bold',...
            'Tag','toolsTxt');

        uicontrol('Style','popupmenu',...
            'Units','normalized',...
            'Position',[Pos+[0.01,-0.05],0.25,0.13],...
            'String',toolStr,...
            'FontUnits','normalized',...
            'FontSize',0.15,...
            'Tag','toolsMenu',...
            'Callback',@toolSelect_Callback);

        Pos = UIPos(6,:,3);
        uicontrol('Style','text',...
            'String','PreSet',...
            'Units','normalized',...
            'Position',[Pos+[0.015 0.07],0.22,0.06],...
            'HorizontalAlignment','Center',...
            'FontUnits','normalized',...
            'FontSize',0.35,...
            'FontWeight','bold',...
            'Tag','preSetTxt');
        uicontrol('Style','pushbutton',...
            'String','Save',...
            'Units','normalized',...
            'Position',[Pos+[0.02 0.06],0.1,0.04],...
            'FontUnits','normalized',...
            'FontSize',0.5,...
            'Tag','preSetSave',...
            'Callback',{@savePreSet_Callback});
        uicontrol('Style','pushbutton',...
            'String','Load',...
            'Units','normalized',...
            'Position',[Pos+[0.13 0.06],0.1,0.04],...
            'FontUnits','normalized',...
            'FontSize',0.5,...
            'Tag','preSetLoad',...
            'Callback',{@loadPreSet_Callback});

        % - Add Cineloop control if DisplayWindow contains more than one frame.
        if nfrms > 1
            Pos = UIPos(1,:,3);
            uicontrol('Style','text','String','CineLoop',...
                'Units','normalized',...
                'Position',[Pos+SG.TO,SG.TS],...
                'FontUnits','normalized',...
                'FontSize',0.8,...
                'FontWeight','bold',...
                'Tag','cinetxt');
            cineSldr = uicontrol('Style','slider',...
                'Max',nfrms,'Min',1,'Value',nfrms,...
                'SliderStep',[1/(nfrms-1) 1/(nfrms-1)],...
                'Units','normalized',...
                'Position',[Pos+SG.SO,SG.SS],...
                'BackgroundColor',Bkgrnd-0.05,...
                'Interruptible','off',...
                'BusyAction','cancel',...
                'Tag','CLSlider');
            % -- Cineloop number
            cineValue = uicontrol('Style','edit','String',num2str(nfrms,'%2.0f'),...
                'Units','normalized',...
                'Position',[Pos+[0.02 0.03],0.08,0.03],...
                'Tag','CLValue',...
                'Callback',{@cineValue_Callback, cineSldr, freezeBtn},...
                'BackgroundColor',Bkgrnd+0.1);
            % -- Cineloop save button
            cineSave = uicontrol('Style','togglebutton',...
                'String','Save',...
                'Units','normalized',...
                'Position',[Pos+[0.12 0.023],0.10,0.04],...
                'FontUnits','normalized',...
                'FontSize',0.5,...
                'Enable','off',...
                'Tag','cineSave',...
                'Callback',{@cineSave_Callback, freezeBtn});
            cineSldr.Callback = {@cineSldr_Callback, cineValue, freezeBtn};
            freezeBtn.Callback = {@freeze_Callback, cineSldr, cineValue, cineSave};
        end

        set(f,'CloseRequestFcn',{@closefunc,hvTimer});

        % Make the GUI visible, unless the call has requested that it be hidden.
        visibility = 'on';
        if evalin('base', 'exist(''Mcr_GuiHide'', ''var'')') && ...
                evalin('base', 'Mcr_GuiHide')
            visibility = 'off';
        end
        set(f,'Visible', visibility);
    end


    function success = applySliderValue(slider, value)
        if ~isempty(slider) && isvalid(slider)
            slider.Value = value;
            success = true;
        else
            success = false;
        end
    end


    function success = applyEditValue(edit, value)
        if ~isempty(edit) && isvalid(edit)
            edit.String = num2str( value );
            success = true;
        else
            success = false;
        end
    end



    function closefunc(source,~,hvTimer)
        assignin('base', 'vsExit', 1);
        % VTS-1306 If hvTimer has been initialized as a timer (script using TPC
        % Profile 5), it will be an object with class name timer so the
        % lines below will stop and delete it.  For all other scripts,
        % "hvTimer" will have been created as a placeholder graphics object
        % with no properties so the expression "isa(hvTimer, 'timer')" will
        % evaluate to false and the stop, delete commands will be skipped.
        if isa(hvTimer, 'timer')
            stop(hvTimer);
            delete(hvTimer);
        end
        delete(source);
        hFigure = findobj('Type','figure');
        if ~isempty(hFigure)
            close(hFigure(contains(get(hFigure,'Tag'),'Msgbox')));
        end
    end

    function updateTGCSliders( nTGC, tgcSldrs, tgcAllSldr)
        TGC = evalin('base','TGC');
        TGCValues = TGC(nTGC).CntrlPts/1023;
        SP = TGCValues/SPAll(nTGC);
        for i=1:length(SP)
            set(tgcSldrs(i), 'Value', TGCValues(i));
            assignin('base', sprintf('tgc%d',i), TGC(nTGC).CntrlPts(i));
        end
        set(tgcAllSldr,'Value',TGCAll(nTGC));
        assignin('base','nTGC',nTGC);
        assignin('base','action','tgc');
    end

% TGC Callback functions
%   The array SP keeps track of the slider positions before applying the SPAll(nTGC) gain factor.
%   This allows returning saturated TGC sliders to original gain curve values if the SPAll(nTGC)
%   gain is lowered.
    function tgcSelect_Callback(hObject,~,tgcSldrs,tgcAllSldr)
        nTGC = get(hObject,'Value');
        updateTGCSliders( nTGC, tgcSldrs, tgcAllSldr);
    end
    function tgcSldr_Callback(tgcSldr,~,tgcSldrNum)
        nTGC = evalin('base', 'nTGC');
        tgcValue = get(tgcSldr,'Value');
        SP(tgcSldrNum) = tgcValue/SPAll(nTGC);
        assignin('base', sprintf('tgc%d', tgcSldrNum), min(1023,1023*tgcValue));
        assignin('base', 'action', 'tgc');
    end
    function tgcAllSldr_Callback(tgcAllSldr,~,tgcSldrs)
        nTGC = evalin('base', 'nTGC');
        TGCAll(nTGC) = get(tgcAllSldr,'Value');
        SPAll(nTGC) = 1.05*TGCAll(nTGC)*TGCAll(nTGC) + 1.95*TGCAll(nTGC) + 1;  % Convert to gain factor between 0.25 and 4.0
        for i=1:length(SP)
            tgcValue = SP(i)*SPAll(nTGC);
            set(tgcSldrs(i), 'Value', min(1.0,tgcValue));
            assignin('base', sprintf('tgc%d',i), min(1023,1023*tgcValue));
        end
        assignin('base', 'action', 'tgc');
    end

% RcvData Loop Callback
    function rcvdataloop_Callback(source,~)
        assignin('base','rloopButton',get(source,'Value'));
        assignin('base', 'action', 'rcvloop');
    end

% Simulate Callback
    function simulate_Callback(source,~)
        assignin('base','simButton',get(source,'Value'));
        assignin('base', 'action', 'simulate');
    end

% Freeze Callback
    function freeze_Callback(source,~,varargin)
        if evalin('base','isequal(initialized,0)')
            return
        end
        frzstate = get(source,'Value');
        assignin('base','freeze',frzstate);
        if nargin == 5
            cineSldr = varargin{1};
            cineValue = varargin{2};
            cineSave = varargin{3};
            if frzstate == 0
                set(cineSldr, 'Value', nfrms);
                set(cineValue,'String',num2str(nfrms,'%2.0f'));
                set(cineSave,'Enable','off');
            else
                set(cineSave,'Enable','on');
            end
        end
    end

% Tool Select Callback
    function toolSelect_Callback(hObject,~)
        toolValue = get(hObject,'Value');
        switch toolValue
            case 1
                close(findobj('tag','filterTool'));
                close(findobj('tag','TXPD'));
                close(findobj('tag','ColorMapTool'));
                close(findobj('tag','ProcessTool'));
            case 2
                filterTool;
            case 3
                showTXPD
            case 4
                saveRF;
            case 5
                PTool;
            case 6
                ColorMapTool;
        end
    end

% High Voltage 1 Slider Callback
    function hv1Sldr_Callback(hv1Sldr,~,hv1Value)
        if trackP5
            hv = str2double(get(hv1Value,'String'));
            set(hv1Sldr,'Value',hv);
            return
        end
        hv = get(hv1Sldr,'Value');
        hvset = setTpcVoltage(hv, 1);
        % Since requested value is not necessarily the value that
        % was obtained, we set the slider to the resulting value.
        set(hv1Sldr,'Value',hvset);
        set(hv1Value,'String',num2str(hvset,'%.1f'));
        evalin('base','tStartHvSldr = tic;'); % set time the slider was moved for error suppression.
    end

% High Voltage 1 Value Callback
    function hv1Value_Callback(hv1Value,~,hv1Sldr)
        if trackP5
            hv = get(hv1Sldr,'Value');
            set(hv1Value,'String',num2str(hv,'%.1f'));
            return
        end
        hv = str2double(get(hv1Value,'String'));
        % Protect against bad user input (e.g."1.6.6").
        if isnan(hv)
            hv = get(hv1Sldr,'Value');
        end
        % Don't allow setting hv outside slider's Min/Max range.
        sliderMin = get(hv1Sldr, 'Min');
        sliderMax = get(hv1Sldr, 'Max');
        if hv < sliderMin, hv = sliderMin; end
        if hv > sliderMax, hv = sliderMax; end
        hvset = setTpcVoltage(hv, 1);
        % Since requested value is not necessarily the value that
        % was obtained, we set the slider to the resulting value.
        set(hv1Sldr,'Value',hvset);
        set(hv1Value,'String',num2str(hvset,'%.1f'));
        evalin('base','tStartHvSldr = tic;'); % set time the slider was moved for error suppression.
    end

% High Voltage 2 Slider Callback
    function hv2Sldr_Callback(hv2Sldr,~,hv2Value,hv1Sldr,hv1Value)
        if hvps2 == 0, return, end
        hv = get(hv2Sldr,'Value');
        TPC = evalin('base', 'TPC');
        hvset = setTpcVoltage(hv, hvps2);
        % Since requested value is not necessarily the value that
        % was obtained, we set the slider to the resulting value.
        set(hv2Sldr,'Value',hvset);
        set(hv2Value,'String',num2str(hvset,'%.1f'));
        if hvps2==5 && trackP5 % If profile 5 active check for tracking feature enabled
            hv = max(get(hv1Sldr, 'Min'), hvset); % don't go below slider minimum
            hvset = setTpcVoltage(hv, trackP5);
            if trackP5 % only update hv1 gui if trackP5 points to it and not some other profile
                set(hv1Sldr,'Value',hvset);
                set(hv1Value,'String',num2str(hvset,'%.1f'));
            end
        end
        evalin('base','tStartHvSldr = tic;'); % set time the slider was moved for error suppression.
    end

% High Voltage 2 Value Callback
    function hv2Value_Callback(hv2Value,~,hv2Sldr,hv1Sldr,hv1Value)
        if hvps2 == 0, return, end
        hv = str2double(get(hv2Value,'String'));
        % Protect against bad user input (e.g."1.6.6").
        if isnan(hv)
            hv = get(hv2Sldr,'Value');
        end
        % Don't allow setting hv outside slider's Min/Max range.
        sliderMin = get(hv2Sldr, 'Min');
        sliderMax = get(hv2Sldr, 'Max');
        if hv < sliderMin, hv = sliderMin; end
        if hv > sliderMax, hv = sliderMax; end
        TPC = evalin('base', 'TPC');
        hvset = setTpcVoltage(hv, hvps2);
        % Since requested value is not necessarily the value that
        % was obtained, we set the slider to the resulting value.
        set(hv2Sldr,'Value',hvset);
        set(hv2Value,'String',num2str(hvset,'%.1f'));
        if hvps2==5 && trackP5 % If profile 5 active check for tracking feature enabled
            hv = max(get(hv1Sldr, 'Min'), hvset); % don't go below slider minimum
            hvset = setTpcVoltage(hv, trackP5);
            if trackP5 == 1 % only update hv1 gui if trackP5 points to it and not some other profile
                set(hv1Sldr,'Value',hvset);
                set(hv1Value,'String',num2str(hvset,'%.1f'));
            end
        end
        evalin('base','tStartHvSldr = tic;'); % set time the slider was moved for error suppression.
    end

% Speed Slider Callback
    function speedSldr_Callback(speedSldr,~,speedValue)
        sv = get(speedSldr,'Value');
        assignin('base', 'speedCorrect', sv);
        assignin('base', 'action', 'speed');
        set(speedValue,'String',num2str(sv,'%1.3f'));
    end

% Speed Value Callback
    function speedValue_Callback(speedValue,~,speedSldr)
        sv = str2double(get(speedValue,'String'));
        svMin = get(speedValue, 'Min');
        svMax = get(speedValue, 'Max');
        if svMin<=sv && sv<=svMax
            assignin('base', 'speedCorrect', sv);
            assignin('base', 'action', 'speed');
            set(speedSldr,'Value', sv);
        else
            sv = get(speedSldr,'Value');
            set(speedValue,'String',num2str(sv,'%1.3f'));
        end
    end

% Cineloop Callback
    function cineSldr_Callback(cineSldr,~,cineValue,freezeBtn)
        if get(freezeBtn,'Value')==0    % no action if not in freeze
            set(cineSldr,'Value',nfrms);  % reset slider to end
            return
        end
        cv = round(get(cineSldr,'Value'));
        Control.Command = 'cineDisplay';
        Control.Parameters = {'displayWindow', 1, 'frameNumber', cv};
        runAcq(Control);
        set(cineValue,'String',num2str(cv,'%2.0f'));
    end

% cineValue Callback
    function cineValue_Callback(cineValue,~,cineSldr,freezeBtn)
        if get(freezeBtn,'Value')==0  % no action if not in freeze
            return
        end
        cv = round(str2double(get(cineValue,'String')));
        if 1<=cv && cv<=nfrms
            set(cineSldr, 'Value', cv);
            Control.Command = 'cineDisplay';
            Control.Parameters = {'displayWindow', 1, 'frameNumber', round(cv)};
            runAcq(Control);
        else
            cv = get(cineSldr,'Value');
            set(cineValue,'String',num2str(cv,'%2.0f'));
        end
    end

% cineSave Callback
    function cineSave_Callback(src,~,freezeBtn)
        if get(freezeBtn,'Value')==0   % no action if not in freeze
            set(src,'Value',0);
            return
        end
        islinux = ~ismac && isunix;
        savingOptions = {'*.mp4';'*.avi'};
        Resource = evalin('base','Resource');

        % TODO this should be made chnageable, right now we only support
        % one displaywindow
        windowIndex = 1;
        DisplayWindow = Resource.DisplayWindow(windowIndex);
        if isfield(DisplayWindow,'Type')&&~isempty(DisplayWindow.Type)&&strcmp(DisplayWindow.Type,'Verasonics')
            Type = 'Verasonics';
        else
            Type = 'Matlab';
            fig = DisplayWindow.figureHandle;
            if islinux   % The mpeg video encoder is not available in Linux Matlab
                savingOptions = '*.avi';
            end
        end

        filename = datestr(now,'dd-mmmm-yyyy_HH-MM-SS');
        [fn,pn,indx] = uiputfile(savingOptions,'Save cineloop as',filename);
        if ~isequal(fn,0) % fn will be zero if user hits cancel
            fn = strrep(fullfile(pn,fn), '''', '''''');
        else
            disp('The cineloop is not saved.');
            return
        end

        frameRate = 7;
        n = 1; % n keeps track of frame no. from first to last.
        switch Type
            case 'Verasonics'
                import com.verasonics.viewer.tools.movie.FrameCaptureTool
                import com.verasonics.imageprocessing.formats.ffmpeg.FfmpegMovie
                import com.verasonics.imageprocessing.formats.attributes.*
                % relying on the workspace variables is the killer caused a
                % nasty bug!
                imageViewer = evalin('base','imageViewer');
                % we have to select the right image viewer here
                imageViewer = imageViewer(windowIndex);

                % number of digits is required for createMovieFromImageFiles
                frameDigit = numel(num2str(DisplayWindow.numFrames));
                folderName = 'MovieCreation';
                % If the temp folder exists, empty it.
                % Else, create it.
                if exist(folderName, 'dir')
                    delete(fullfile(folderName, '*.jpg'));
                else
                    mkdir(folderName);
                end

                captureTool = FrameCaptureTool();
                captureTool.enableCapturing(imageViewer.getImageCanvas,1,folderName,'frame',DisplayWindow.numFrames+1,1);

                while n <= DisplayWindow.numFrames
                    captureTool.beginFrameCapture();
                    Control.Command = 'cineDisplay';
                    Control.Parameters = {'displayWindow', 1, 'frameNumber', n};
                    runAcq(Control);
                    n = n+1;
                    % 2019-01-28, the software will hang if the pause duration is less than 0.04 s
                    pause(0.05)
                    captureTool.endFrameCapture();
                end

                captureTool.disableCapturing();

                FfmpegMovie.registerImageFormat()
                switch indx
                    case 2
                        movieFormat = MovieContainer.avi;
                        encoder = VideoEncoder.mpeg;
                        quality = MovieEncodingMode.maxBitrate;
                    otherwise
                        movieFormat = MovieContainer.mp4;
                        encoder = VideoEncoder.h264;
                        quality = MovieEncodingMode.constantQuality;
                end

                newMovieAttributes = MovieAttributes(movieFormat,...
                    encoder,...
                    quality,...
                    MovieQualityLevel.medium,...
                    frameRate,...
                    -1, -1,...
                    MoviePlayMode.playForward, 1);

                ffmpegMovie = FfmpegMovie();
                success = ffmpegMovie.createMovieFromImageFiles(fullfile(pwd, folderName), ...
                                                                sprintf('frame[0-9]{%d}+.jpg', ...
                                                                frameDigit), ...
                                                                newMovieAttributes, ...
                                                                fn(1:end-4), ...
                                                                false);
            case 'Matlab'
                success = 0;
                F = repmat(struct('cdata',[],'colormap',[]),1,DisplayWindow.numFrames);
                while n <= DisplayWindow.numFrames
                    Control.Command = 'cineDisplay';
                    Control.Parameters = {'displayWindow', 1, 'frameNumber', n};
                    runAcq(Control);
                    F(n) = getframe(fig);
                    n = n+1;
                end

                % Initialize v to empty
                v = [];
                switch indx
                    case 2
                        v = VideoWriter(fn);
                    otherwise
                        if islinux
                            v = VideoWriter(fn);
                        else
                            v = VideoWriter(fn,'MPEG-4');
                        end
                end

                if ~isempty(v)
                    v.FrameRate = frameRate;
                    open(v);
                    writeVideo(v,F);
                    close(v);
                    success = 1;
                end
        end

        if success
            fprintf('The cineloop has been saved at %s \n',fn);
            set(src,'Value',0);
        else
            fprintf(2, 'There was a failure in creating the cineloop movie!\n');
        end
    end

% HV timer callback for actual push capacitor voltage; used for both
% internal and external profile 5 power supply configurations.
    function hvTimer_Callback(~,~,hv2Sldr,hv2Actual)
        [Result,extCapVoltage] = getHardwareProperty('TpcExtCapVoltage');
        if ~strcmp(Result,'Success')
            % If a hardware fault has occurred, an error was already reported,
            % don't also report this because it confuses the root cause issue.
            if com.verasonics.hal.faults.Faults.getNumberOfFaults() == 0
                error('VSX: Error from getHardwareProperty call to read push capacitor Voltage.');
            end
        end
        set(hv2Actual,'String',num2str(extCapVoltage,'%.1f'));
        if 1-(extCapVoltage/get(hv2Sldr,'Value')) > 0.20 % if 20% low
            set(hv2Actual,'BackgroundColor',[0.7,0.7,1.0]);
        else
            set(hv2Actual,'BackgroundColor',[0.8,0.8,0.8]);
        end
    end


% preSet Callbacks

    function stopSequence()
        
        Control = vsv.seq.getBaseComp('Control');
        if isempty(Control)
            Control.Command = 'stopSequence';
        else
            Control(end).Command = 'stopSequence';
        end
        assignin( 'base', 'Control', Control );
        assignin( 'base', 'freeze', 1);
        
        evalin( 'base', 'runAcq(Control);' );
        evalin( 'base', 'Control = struct(''Command'', [], ''Parameters'', []);');
    end

    function savePreSet_Callback(~, ~)

        % new app control support 
        app = vsv.seq.getBaseComp('vsxApplication');
            
        if ~isempty(app)

            stopSequence();
            
            hasStorage = getStorageParameterFromUI(app);
                        
            % create a storage parser
            storageParser = vsv.vsx.storage.StorageParser();
            storageParser.setStorageParameter( hasStorage );
            
            fileName = evalin('base','displayWindowTitle');
            preFix = fullfile(vpfRoot,'MatFiles',fileName);
            preSet.preFix = preFix;
            preSet.storageParser = storageParser;
            
            % From VSX: Trans, Resource, TW, TX and Event are required, but
            % other structures must be checked for existence before saving.
            Resource = evalin('base','Resource');
            preSet.SWversion = Resource.SysConfig.SWversion;

            % Trans, all fields are required
            preSet.Trans = evalin('base','Trans');

            % TW and TX, remove TXPD to avoid a huge preset file size
            preSet.TW = evalin('base','TW');
            TX = evalin('base','TX');
            if isfield(TX,'TXPD')
                for j=1:length(TX)
                    TX(j).hasTXPD = ~isempty(TX(j).TXPD);
                end
                TX = rmfield(TX,'TXPD');
            end
            preSet.TX = TX;

            % added feature to save speed correct factor as well
            if evalin('base','exist(''speedCorrect'',''var'')')
                preSet.speedCorrect = evalin('base','speedCorrect');
            end

            % Receive
            if evalin('base','exist(''Receive'',''var'')')
                preSet.Receive = evalin('base','Receive');
            end

            % Image may not be required
            if evalin('base','exist(''PData'',''var'')')
                preSet.PData = evalin('base','PData');

                % Save colormap
                for winNum = 1:length(Resource.DisplayWindow)
                    preSet.DisplayWindow(winNum).Position = Resource.DisplayWindow(winNum).Position;
                    preSet.DisplayWindow(winNum).ReferencePt = Resource.DisplayWindow(winNum).ReferencePt;
                    preSet.DisplayWindow(winNum).Colormap = Resource.DisplayWindow(winNum).Colormap;
                end

                % Save custom gamma curve if ColorMapTool is used for
                % adjustment
                if evalin('base','exist(''customGamma'',''var'')')
                    preSet.Display.customGamma = evalin('base','customGamma');
                end

            end

            % Recon may not be required
            if evalin('base','exist(''Recon'',''var'')')
                preSet.Recon = evalin('base','Recon');
                preSet.ReconInfo = evalin('base','ReconInfo');
            end

            % Process may not be used
            if evalin('base','exist(''Process'',''var'')')
                preSet.Process = evalin('base','Process');
                % if PTool has been called, compFactorAll should be saved
                if evalin('base','exist(''compFactorAll'',''var'')')
                    preSet.compFactorAll = evalin('base','compFactorAll');
                end
            end

            % persf and perst for Doppler, if exist
            if evalin('base','exist(''persf'',''var'')')
                preSet.Doppler.persf = evalin('base','persf');
            end
            if evalin('base','exist(''persp'',''var'')')
                preSet.Doppler.persp = evalin('base','persp');
            end

            % for Personal definition
            if evalin('base','exist(''P'',''var'')')
                preSet.P = evalin('base','P');
            end

            % for backward compability
            if evalin('base','exist(''SFormat'',''var'')')
                preSet.SFormat = evalin('base','SFormat');
            end

            % TGC
            if evalin('base','exist(''TGC'')')
                preSet.TGCparam.TGC = evalin('base','TGC');
                preSet.TGCparam.nTGC = evalin('base','nTGC');
                preSet.TGCparam.TGCAll = TGCAll;
                preSet.TGCparam.SPAll = SPAll;
            end

            % TPC structure always exists (created by VSX if not by SetUp script)
            preSet.TPC = evalin('base','TPC'); 

            [fn,pn] = uiputfile('*.mat','Save preSet as',[preFix,'_preSet']);
            if ~isequal(fn,0) % fn will be zero if user hits cancel
                fn = strrep(fullfile(pn,fn), '''', '''''');
                % Avoid warning about 'creation of very large files'
    %             warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
                save(fn, 'preSet' );
    %             warning('on', 'MATLAB:Figure:FigureSavedToMATFile');
                fprintf('The preSet has been saved at %s\n',fn);
            else
                disp('The preSet is not saved.');
            end
        else
            disp('The preSet is not saved because vsxApplication is not present in VSX.');
        end

    end

    function applyPreset(preSet, Resource)
    % From the VSX: Trans, Resource, TW, TX and Event are required,
    % but other structures must be checked before loading
        
        Trans = preSet.Trans; assignin('base','Trans',Trans);

        % assign TX and TW
        TX = preSet.TX; 
        assignin('base','TX',TX);
        
        TW = preSet.TW; 
        assignin('base','TW',TW);

        if isfield(preSet, 'Receive')
            Receive = preSet.Receive; 
            assignin('base','Receive',Receive);
        end

        % Image may not be required
        if isfield(preSet, 'PData')
            PData = preSet.PData; assignin('base','PData',PData);
            for winNum = 1:length(Resource.DisplayWindow)
                Resource.DisplayWindow(winNum).Position = preSet.DisplayWindow(winNum).Position;
                Resource.DisplayWindow(winNum).ReferencePt = preSet.DisplayWindow(winNum).ReferencePt;
                Resource.DisplayWindow(winNum).Colormap = preSet.DisplayWindow(winNum).Colormap;
            end
            % Restore custom gamma curve if it's saved in preSet
            if isfield(preSet,'Display')
                if isfield(preSet.Display,'customGamma')
                    assignin('base','customGamma',preSet.Display.customGamma);
                end
            end
            assignin('base','Resource',Resource);
        end

        if isfield(preSet, 'SFormat')
            SFormat = preSet.SFormat; assignin('base','SFormat',SFormat);
        end

        % Recon may not be required
        if isfield(preSet, 'Recon')
            Recon = preSet.Recon; assignin('base','Recon',Recon);
            ReconInfo = preSet.ReconInfo; assignin('base','ReconInfo',ReconInfo);
        end

        % Process may be used
        if isfield(preSet, 'Process')
            Process = preSet.Process; assignin('base','Process',Process);
            if isfield(preSet,'compFactorAll')
                assignin('base','compFactorAll',preSet.compFactorAll);
            end
        end

        % for Personal definition
        if isfield(preSet, 'P')
            P = preSet.P; assignin('base','P',P);
        end

        % TGC and TGCAll
        if isfield(preSet, 'TGCparam')
            TGC = preSet.TGCparam.TGC; assignin('base','TGC',TGC);
            nTGC = preSet.TGCparam.nTGC; assignin('base','nTGC',nTGC);
            TGCAll = preSet.TGCparam.TGCAll;
            SPAll = preSet.TGCparam.SPAll;
            TGCValues = TGC(nTGC).CntrlPts/1023;
            SP = TGCValues/SPAll(nTGC);
            for i=1:length(TGC(nTGC).CntrlPts)
                assignin('base',sprintf('tgc%d',i),TGCValues(i));
            end
            
            % apply the value
            if ~isempty(tgcSelect) ...
                    && isa(tgcSelect, 'matlab.ui.control.UIControl') ...
                    && isvalid(tgcSelect)
                
                tgcSelect.Value = nTGC;
                % call the callback
                tgcSelect_Callback(tgcSelect,[],tgcSldrs,tgcAllSldr);
            else
                updateTGCSliders( nTGC, tgcSldrs, tgcAllSldr);
            end
            
            
        end

        TPC = preSet.TPC; assignin('base','TPC',TPC);
        % High Voltage 1
        if TPC(1).inUse
            evalin('base','tStartHvSldr = tic;'); % set time the slider was moved for error suppression.
            
            applySliderValue( hv1Sldr, TPC(1).hv );
            applyEditValue( hv1Value, TPC(1).hv );
            
            setTpcVoltage(TPC(1).hv, 1);
        end

        % High Voltage 2
        if hvps2 ~= 0
            
            applySliderValue( hv2Sldr, TPC(hvps2).hv );
            applyEditValue( hv2Value, TPC(hvps2).hv );
            
            setTpcVoltage(TPC(hvps2).hv, hvps2);
        end

        % If PTool window is open, reload it
        hPTool = findobj('tag','ProcessTool');
        if ishandle(hPTool)
            posPTool = get(hPTool,'position');
            PTool;
            set(hPTool,'position',posPTool);
        end

        % If ColorMapTool window is open, reload it
        hCMTool = findobj('tag','ColorMapTool');
        if ishandle(hCMTool)
            posGSTool = get(hCMTool,'Position');
            ColorMapTool;
            set(hCMTool,'position',posGSTool)
        end

        % persf and perst for Doppler, if exist
        if isfield(preSet, 'Doppler')
            assignin('base','persf',preSet.Doppler.persf);
            assignin('base','persp',preSet.Doppler.persp);
        end

        Control = [];
        Control.Command = 'update&Run';
        Control.Parameters = {'PData','InterBuffer','ImageBuffer',...
            'DisplayWindow','Parameters','Trans','Media','TW','TX',...
            'Receive','TGC','TPC','Recon','Process'};

        % Recompute TXPD and put in workspace, if necessary.
        if isfield(preSet.TX, 'hasTXPD')
            h = waitbar(0,'Program TX parameters, please wait!');
            n = 1;
            for i = 1:length(Recon)
                txnum = unique([ReconInfo(Recon(i).RINums(1):Recon(i).RINums(end)).txnum]);
                pdatanum = Recon(i).pdatanum;
                for j = 1:length(txnum)
                    if TX(j).hasTXPD
                        TX(txnum(j)).TXPD = computeTXPD(TX(txnum(j)),PData(pdatanum));
                        waitbar(n/size(TX,2));
                        n = n+1;
                    end
                end
            end
            close(h);
            TX = rmfield(TX, 'hasTXPD');
            assignin('base','TX',TX);
        end
        
        % now we need to reapply the changes from the UI controls
        if isfield(preSet, 'storageParser')
            
            % storageParser = vsv.vsx.storage.StorageParser
            storageParser = preSet.storageParser;
            
            % new app control support 
            app = vsv.seq.getBaseComp('vsxApplication');
            
            % get new ui controls
            hasStorage = getStorageParameterFromUI(app);
            
            % we need to assign in base to make sure that the controls will
            % take all the changes
            assignin('base','Control',Control);
            
            % reapply to current GUI components what was saved in the
            % storage parameter. The parser will find the GUI components by
            % their label definition (the label defines the ID of the
            % storage parameter, e.g., Range (mm) ) 
            storageParser.applyStorageParameter(hasStorage);
            
            % make sure we do not have multiple update and runs and make
            % sure the update and run is all in the right order. This needs
            % to be tested and might not work in any case
            vsv.vsx.command.SimpleCommandQue.mergeBaseControl();
            Control = vsv.seq.getBaseComp('Control'); 
        end

        % Add 'set&Run' to set colormap.
        for winNum = 1:length(Resource.DisplayWindow)
            n=numel(Control)+1;
            Control(n).Command = 'set&Run'; %#ok
            Control(n).Parameters = {'DisplayWindow',winNum,...
                'colormap',preSet.DisplayWindow(winNum).Colormap};%#ok
        end
        
        if isfield( preSet, 'speedCorrect' )
            speedCorrect = preSet.speedCorrect;
            % OK so this cannot be applied :( because we can only apply one
            % action. grr
            speedSldr.Value   = speedCorrect;
            speedValue.String = num2str(speedCorrect);
            
            n = numel(Control)+1;
            
            Resource = vsv.seq.getBaseComp('Resource');
            
            if ~isempty(Resource)
                Resource.Parameters.speedCorrectionFactor = speedCorrect;
            
                assignin('base','Resource',Resource);
                Control(n).Command = 'set&Run';
                Control(n).Parameters = {'Parameters',1,...
                                     'speedCorrectionFactor',speedCorrect};
            end
            
        end
        
        assignin('base','Control',Control);
        assignin('base', 'action', 'displayChange');
    end

    function loadPreSet_Callback(~, ~)

        if evalin('base','isequal(initialized,0)')
            return
        end

        stopSequence();

        fileName = evalin('base','displayWindowTitle');
        preFix = fullfile(vpfRoot,'MatFiles',fileName);

        [FileName,PathName] = uigetfile('*.mat', ...
            ['Select the preSet file for ',fileName], ...
            [preFix '_preSet']);
        
        
        if ischar( FileName ) && ~isempty(FileName)
            
            buffer = load(fullfile(PathName, FileName));
            
            % Exit if user pressed Cancel.
            % Following has been reformatted to follow single exit strategy
            % check whether the mat file contains a variable called preSet
            if isfield(buffer,'preSet')

                preSet = buffer.preSet;

                Resource = evalin('base','Resource');
                SWversion = Resource.SysConfig.SWversion;
                if preSet.SWversion(1) == SWversion(1) && ...
                        preSet.SWversion(2) == SWversion(2)

                    % Make sure that the transducer is correct
                    if strcmp(evalin('base', 'Trans.name'), preSet.Trans.name)
                        applyPreset(preSet, Resource);
                    else
                        warndlg('Incorrect Transducer!');
                    end

                else
                    warndlg({sprintf('Presets were created for release %d.%d', preSet.SWversion(1:2))...
                        sprintf('Please recreate your presets file for release %d.%d', SWversion(1:2))}, ...
                        'Compatibility Error', 'warn')

                end                
            else
                warndlg('The selected MAT file does not contain presets!');
            end
        end


    end

    % close PTool and ColorMapTool if open
    delete(findobj('tag','ProcessTool'));
    delete(findobj('tag','ColorMapTool'));

    % if saveRF is in the temp directory, delete it
    saveRFFilePath = fullfile(vsTempDir,'saveRF.m');
    if exist(saveRFFilePath,'file')
        delete(saveRFFilePath);
    end

end

function hasStorage = getStorageParameterFromUI(app)
    % get new ui controls
    uiManager = app.getUserControlManager();
    hasStorage = uiManager.UiControls(:,1);

    % VTS-2203 check whether ui control actually is a storage
    % parameter, this might not be the case for example if ui
    % component is a button.
    doRemove = false(size(hasStorage));
    for i = 1:numel(hasStorage)
        if ~isa( hasStorage{i}, 'vsv.seq.storage.HasStorageParameter' )
            doRemove(i) = true;
        end
    end
    hasStorage(doRemove) = [];
end

function hvset = setTpcVoltage(hv, ntpc)
    % Attempt to set high voltage.
    % On error, setTpcProfileHighVoltage() returns voltage range minimum.
    [result, hvset] = setTpcProfileHighVoltage(hv, ntpc);
    if ~strcmpi(result, 'Success')
        % ERROR!  Failed to set high voltage.
        error('ERROR!  Failed to set Verasonics TPC high voltage for profile %d because \"%s\".', ntpc, result);
    end
end
