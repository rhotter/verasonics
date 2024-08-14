function EventAnalysisTool(varargin)
% Copyright 2001-2018 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   tool for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%   Starting from software release 4.2, the initializeOnly feature and
%   VsUpdate will not be applied by default. The initializeOnly mode can be
%   executed by clicking the 'initializeOnly mode' push button or type
%   'EventAnalysisTool i' in the command window for more information (slower execution).
%
% File name: EventAnalysisTool.m - A tool to analysize and display important parameters in Event Sequence.
%
% 20191112 changed default behavior and added a push button for initializeOnly mode

persistent ADhndl Apodhndl Delayhndl  ...
    AX H1 H2 rcvReset txReset PersistAcqDispArray...
    imageholdReset imageholdhndl imagehold reconReset processReset...
    seqControlReset reconVal reconEvent ReconInfoTable ReconInfos...
    tableReconEvent  ReconInfoColnames InOutTable ReconInfoWarning...
    notxhndl messageFont centerPos tableFontSize;

%Set default values
rcvReset = 0;
txReset = 0;
reconReset = 0;
processReset = 0;
seqControlReset = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the master figure size
delete(findobj('tag','EATool'));

sa = figure('Visible','on',...
    'Units','normalized',...
    'Position',[.05,.1,.9,.8],...
    'NumberTitle','off',...
    'Color',[.7 .7 .7],...
    'MenuBar','none', ...
    'Name','Event Analysis Tool',...
    'tag','EATool');

set(sa,'CloseRequestFcn',{@closefunc});

% These are settings to use for any OS other than Windows
titleFont = .7; % font for titles in gui window
textFont = .65; % font for gui text labels
messageFont = .55;
tableFontSize = 12;

hori1x = .015; % reference x position for gui items: table and text output
vert1x = .95; % reference y position for first set of text parameter output
% vert2x = 730; % reference y position for second set of text parameter output
% vert3x = 710; % reference y position for third set of text parameter output

%Default positions for checkboxes
checkboxhori = .265;
checkboxvert = .9;

%Default TX/Rcv plot Positions
plot1hori = .3;
plot2hori = .4;
plot1vert = .9;
plot2vert = .425;

%Center Position
centerPos = [.1 .65 .8 .04];

% dispmessage
dispmessage = uicontrol('Parent',sa,'Style','text',...
    'FontUnits','normalized',...
    'FontSize',messageFont,...
    'Units','normalized',...
    'Position',centerPos,...
    'HorizontalAlignment','center',...
    'BackgroundColor',[.7 .7 .7],...
    'HandleVisibility','off',...
    'FontWeight','bold',...
    'Visible','off',...
    'String','No Script Variables Found in Matlab Workspace...');

if ~evalin('base','exist(''Resource'',''var'');')
    dispmessage.String = 'No Resource Variable Found in Matlab Workspace...';
    dispmessage.Visible = 'on';
    return
end

try
    %Event Table Preallocation
    %Evaluate in Event Structure
    EventStruct = evalin('base','Event');
    EventSize = size(EventStruct,2); %Horizontal Direction
catch
    %Check for vaiables in the matlab workspace
    if exist('EventStruct','var') == 0
        dispmessage.Visible = 'on';
        return
    end
end

%%Preallocate Event Cell array
EventSequence = cell(6,EventSize); %Get the total number of events in the sequence
for i = 1:EventSize
    EventSequence{1,i} = num2str(EventStruct(1,i).info);
    EventSequence{2,i} = num2str(EventStruct(1,i).tx);
    EventSequence{3,i} = num2str(EventStruct(1,i).rcv);
    EventSequence{4,i} = num2str(EventStruct(1,i).recon);
    EventSequence{5,i} = num2str(EventStruct(1,i).process);
    EventSequence{6,i} = num2str(EventStruct(1,i).seqControl);
end

%Replace white space in variables with 2 or more arguments
EventSequence = regexprep(EventSequence, '\s\s', ',  ', 'ignorecase');


%%  ==== general debug functionality ====

%% InitializeOnly
if ~isempty(varargin)  % if any argument is defined, run initialization
    needIniz = 1;
else
    needIniz = 0;
end

alreadyIniz = 0;
Resource = evalin('base','Resource');
if isfield(Resource,'VDAS') % already finished initialization?
    if isfield(Resource.VDAS,'exportDelta') && isfield(Resource.VDAS,'cgEnaDma')
        alreadyIniz = 1;
    end
end

if needIniz && ~alreadyIniz
    dispmessage.String = 'Parsing the script using initializeOnly. Please wait...';
    dispmessage.Visible = 'on';
    drawnow
    evalin('base','Resource.Parameters.verbose = 1;');
    evalin('base','Resource.Parameters.initializeOnly = 1;');
    evalin('base','Resource.Parameters.simulateMode = 1;');
    oldMcrHide = [];
    if evalin('base','exist(''Mcr_GuiHide'',''var'');')
        oldMcrHide = evalin('base','Mcr_GuiHide');
    end
    evalin('base','Mcr_GuiHide = 1;');
    if isfield(Resource,'DisplayWindow')
        evalin('base','Resource = rmfield(Resource,''DisplayWindow'');');
    end
    oldfilename = [];
    if evalin('base','exist(''filename'',''var'');')
        oldfilename = evalin('base','filename');
    end
    
    if evalin('base','isfield(TX,''TXPD'')')
        evalin('base','TX = rmfield(TX,''TXPD'');');
    end
    if evalin('base','exist(''UI'',''var'');')
        evalin('base','clear UI');
    end
    evalin('base','filename = ''tempFileForEA'';');
    evalin('base','save(filename,''-v6'');');
    
    try
        evalin('base','VSX');
    catch errMsg
        assignin('base','errMsg',errMsg);
        evalin('base','delete tempFileForEA.mat');
    end
    
    % remove Mcr_GuiHide if it exists
    if ~isempty(oldMcrHide)
        assignin('base','Mcr_GuiHide',oldMcrHide);
    else
        evalin('base','clear Mcr_GuiHide');
    end
    
    % put Resource back
    if isfield(Resource,'DisplayWindow')
        assignin('base','oldDisplayWindow',Resource.DisplayWindow);
        evalin('base','Resource.DisplayWindow = oldDisplayWindow;');
        evalin('base','clear oldDisplayWindow');
    end
    
    % put filename back
    if ~isempty(oldfilename)
        assignin('base','filename',oldfilename);
    end
    
    % set alreadyIniz and delete the temp file
    alreadyIniz = 1;
    if evalin('base','exist(''tempFileForEA.mat'',''file'');')
        evalin('base','delete tempFileForEA.mat');
    end
    
end

%% Other debug features
% Check required structures, from VSX
vars = evalin('base','whos');
RequiredStructs = {'Trans','Resource','TW','TX','Event'};
m = 0;

for i = 1:size(vars,1)
    if any(strcmp(vars(i).name, RequiredStructs)), m=m+1; end
end

if (m~=5), error('Trans, Resource, TW, TX and Event structures are required as minimum. Exiting...\n'); end

% find any empty element, runAcq will have errors if Event has empty
% elements (except for 'info' field)
EventFields = fieldnames(EventStruct);
MsgSize = 1;
warningMsg = cell(20,1); % preallocate 20 cells should be enough
errorColor = cell(20,1);

emptyEvent = cellfun('isempty',EventSequence);
if any(emptyEvent(:))
    [fieldIdx,eventIdx] = find(emptyEvent);
    for Num = 1:length(fieldIdx)
        warningMsg{MsgSize,1} = ['Error!! The ', EventFields{fieldIdx(Num)} ,' field is empty in event(s):  ',num2str(eventIdx(Num))];
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
end

% make sure that each TX, Receive, Recon, Process, and SeqControl is used
% in Event sequence and vise versa
Trans = evalin('base','Trans');
% TX is required
TX = evalin('base','TX');

% check empty feilds
txFields = {'waveform','Apod','Delay'};

for i = 1:length(txFields)
    emptyInd = find(cellfun(@isempty,{TX.(txFields{i})}));
    if ~isempty(emptyInd)
        warningMsg{MsgSize,1} = ['Error!! ',txFields{i}, ' in TX ', num2str(emptyInd), ' is empty!'];
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
end

[TXDiffNum,~] = setdiff(1:1:length(TX),[EventStruct.tx]);
if ~isempty(TXDiffNum)
    if length(TXDiffNum) < 5
        warningMsg{MsgSize,1} = ['Warning!! TX ', num2str(TXDiffNum),' is/are not used in the Event Structure.'];
    else
        warningMsg{MsgSize,1} = ['Warning!! TX ', num2str(TXDiffNum(1:5)),' ... and more are not used in the Event Structure.'];
    end
    errorColor{MsgSize} = 'k';
    MsgSize = MsgSize+1;
end

[TXDiffNum,~] = setdiff([EventStruct.tx],1:1:length(TX));
TXDiffNum(TXDiffNum==0) = [];
if ~isempty(TXDiffNum)
    if length(TXDiffNum) < 5
        warningMsg{MsgSize,1} = ['Error!! TX ', num2str(TXDiffNum),' is/are used in the Event Sequence but not predefined!'];
    else
        warningMsg{MsgSize,1} = ['Error!! TX ', num2str(TXDiffNum(1:5)),' ... and more are used in the Event Sequence but not predefined!'];
    end
    errorColor{MsgSize} = 'r';
    MsgSize = MsgSize+1;
end

if evalin('base', 'exist(''Receive'', ''var'')')
    Receive = evalin('base','Receive');
    
    % check empty feilds
    rcvFields = {'Apod','startDepth','endDepth','TGC','bufnum','framenum','acqNum','mode'};
    
    for i = 1:length(rcvFields)
        emptyInd = find(cellfun(@isempty,{Receive.(rcvFields{i})}));
        if ~isempty(emptyInd)
            warningMsg{MsgSize,1} = ['Warning! ',rcvFields{i}, ' in Receive ', num2str(emptyInd), ' is empty!'];
            errorColor{MsgSize} = 'k';
            MsgSize = MsgSize+1;
        end
    end
    
    rcvNumInEvents = [EventStruct.rcv];  % rcv number shown in Events
    rcvCount = arrayfun(@(x) sum(rcvNumInEvents==x),(1:length(Receive)));
    % rcvCount(1) shows the counts of 0, rcvCount(2) shows the counts of
    % Receive(1) and so on...
    
    repeatInd = find(rcvCount>1); % get the index of Receive used more than once
    if ~isempty(repeatInd)
        dispCount = 1;
        while dispCount < 5
            try rcvInd = repeatInd(dispCount);
                if ~isequal(rcvInd,0)
                    warningMsg{MsgSize,1} = ['Warning!! Receive ', num2str(rcvInd),' is used in multiple events ', num2str(find(rcvNumInEvents == rcvInd)) ,'.'];
                    errorColor{MsgSize} = 'k';
                    MsgSize = MsgSize+1;
                end
            catch
            end
            dispCount = dispCount + 1;
        end
        
        if size(repeatInd,2) > 4
            warningMsg{MsgSize,1} = 'Many identical Receive objects are used in multiple events, please check the Event sequence in the SetUp script';
            errorColor{MsgSize} = 'k';
            MsgSize = MsgSize+1;
        end
    end
    
    [RcvDiffNum,~] = setdiff(1:1:length(Receive),[EventStruct.rcv]);
    if ~isempty(RcvDiffNum)
        if length(RcvDiffNum) < 5
            warningMsg{MsgSize,1} = ['Warning!! Receive ', num2str(RcvDiffNum),' is/are not used in the Event Structure.'];
        else
            warningMsg{MsgSize,1} = ['Warning!! Receive ', num2str(RcvDiffNum(1:5)),' ... and more are not used in the Event Structure.'];
        end
        errorColor{MsgSize} = 'k';
        MsgSize = MsgSize+1;
    end
    
    [RcvDiffNum,~] = setdiff([EventStruct.rcv],1:1:length(Receive));
    RcvDiffNum(RcvDiffNum==0) = [];
    if ~isempty(RcvDiffNum)
        if length(RcvDiffNum) < 5
            warningMsg{MsgSize,1} = ['Error!! Receive ', num2str(RcvDiffNum),' is/are used in the Event Sequence but not predefined!'];
        else
            warningMsg{MsgSize,1} = ['Error!! Receive ', num2str(RcvDiffNum(1:5)),' ... and more are not used in the Event Sequence but not predefined!'];
        end
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
    % Each Receive should be unique, in terms of endDepth, bufnum, framenum,
    % acqNum with mode 0; otherwise the RcvData will be overwritten by later rcv event
    RcvTable = cat(2,[Receive.endDepth]',[Receive.bufnum]',[Receive.framenum]',[Receive.acqNum]',[Receive.mode]');
    [~,~,ic] = unique(RcvTable,'rows','Stable');
    icN = arrayfun(@(x) sum(ic==x),(1:length(RcvTable)));
    multiRcvNum = find(icN > 1);  % get receive number used more than twice
    if ~isempty(multiRcvNum)
        for i = 1:length(multiRcvNum)
            index = find(ismember(ic, multiRcvNum(i)));
            if isequal(RcvTable(index,5),zeros(size(index))) % error only occurs with mode 0 (no accumulation)
                if length(index) < 10
                    warningMsg{MsgSize,1} = ['Warning!!  endDepth, bufnum, framenum, acqNum within Receive ', num2str(index'),' should not be identical in mode 0'];
                else
                    warningMsg{MsgSize,1} = ['Warning!!  endDepth, bufnum, framenum, acqNum within Receive ', num2str(index(1:10)'),' ... and more should not be identical in mode 0'];
                end
                errorColor{MsgSize} = 'k';
                MsgSize = MsgSize+1;
            end
        end
    end
    
elseif nnz([EventStruct.rcv])
    RcvInd = unique(nonzeros([EventStruct.rcv]));
    warningMsg{MsgSize,1} = ['Error!! Receive ', num2str(RcvInd'),' is/are used in the Event Sequence but Receive does not exist!'];
    errorColor{MsgSize} = 'r';
    MsgSize = MsgSize+1;
end

% The Receive.framenum and Receive.endDepth across all ReconInfos of a Recon should be
% identical, otherwise the image would be incorrect
if evalin('base', 'exist(''Recon'', ''var'')')
    Recon = evalin('base','Recon');
    ReconInfo = evalin('base','ReconInfo');
    
    [ReconDiffNum,~] = setdiff(1:1:length(Recon),[EventStruct.recon]);
    if ~isempty(ReconDiffNum)
        if length(ReconDiffNum) < 5
            warningMsg{MsgSize,1} = ['Warning!! Recon ', num2str(ReconDiffNum),' is/are not used in the Event Structure.'];
        else
            warningMsg{MsgSize,1} = ['Warning!! Recon ', num2str(ReconDiffNum(1:5)),' ... and more are not used in the Event Structure.'];
        end
        errorColor{MsgSize} = 'k';
        MsgSize = MsgSize+1;
    end
    
    [ReconDiffNum,~] = setdiff([EventStruct.recon],1:1:length(Recon));
    ReconDiffNum(ReconDiffNum==0) = [];
    if ~isempty(ReconDiffNum)
        if length(ReconDiffNum) < 5
            warningMsg{MsgSize,1} = ['Error!! Recon ', num2str(ReconDiffNum),' is/are used in the Event Sequence but not predefined!'];
        else
            warningMsg{MsgSize,1} = ['Error!! Recon ', num2str(ReconDiffNum(1:5)),' ... and more are used in the Event Sequence but not predefined!'];
        end
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
    for i = 1:length(Recon)
        % Corrent the size of Recon.RINums
        if size(Recon(i).RINums,1)>size(Recon(i).RINums,2)
            Recon(i).RINums = Recon(i).RINums';
        end
        
        % Retrive all receive numbers in ReonInfo
        rcvInd = [ReconInfo(Recon(i).RINums(1):Recon(i).RINums(end)).rcvnum];
        
        RcvDepth = unique([Receive(rcvInd).endDepth]);
        if numel(RcvDepth) > 1
            if length(rcvInd) < 5
                warningMsg{MsgSize,1} = ['Error!! Receive.endDepth across Receive ', num2str(rcvInd),' must be identical for Recon ',num2str(i)];
            else
                warningMsg{MsgSize,1} = ['Error!! Receive.endDepth across Receive ', num2str(rcvInd(1:5)),'... must be identical for Recon ',num2str(i)];
            end
            errorColor{MsgSize} = 'r';
            MsgSize = MsgSize+1;
        end
        
        RcvBuf = unique([Receive(rcvInd).bufnum]);
        if numel(RcvBuf) > 1
            if length(rcvInd) < 5
                warningMsg{MsgSize,1} = ['Error!! Receive.bufnum across Receive ', num2str(rcvInd),' must be identical for Recon ',num2str(i)];
            else
                warningMsg{MsgSize,1} = ['Error!! Receive.bufnum across Receive ', num2str(rcvInd(1:5)),'... must be identical for Recon ',num2str(i)];
            end
            errorColor{MsgSize} = 'r';
            MsgSize = MsgSize+1;
        end
        
        RcvFrame = unique([Receive(rcvInd).framenum]);
        if numel(RcvFrame) > 1
            if length(rcvInd) < 5
                warningMsg{MsgSize,1} = ['Error!! Receive.framenum across Receive ', num2str(rcvInd),' must be identical for Recon ',num2str(i)];
            else
                warningMsg{MsgSize,1} = ['Error!! Receive.framenum across Receive ', num2str(rcvInd(1:5)),'... must be identical for Recon ',num2str(i)];
            end
            errorColor{MsgSize} = 'r';
            MsgSize = MsgSize+1;
        end
    end
    
    % ReconInfo needs to be unique in multiple Recons
    AllReconInfo = unique([Recon.RINums]);
    reUsedRINum = find(arrayfun(@(x) sum([Recon.RINums]==x),AllReconInfo)>1);
    if ~isempty(reUsedRINum)
        warningMsg{MsgSize,1} = ['Error!! ReconInfo ', num2str(reUsedRINum),' can not be used in more than one Recon!'];
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
    
elseif nnz([EventStruct.recon])
    ReconInd = unique(nonzeros([EventStruct.recon]));
    warningMsg{MsgSize,1} = ['Error!! Recon ', num2str(ReconInd'),' is/are used in the Event Sequence but Recon does not exist!'];
    errorColor{MsgSize} = 'r';
    MsgSize = MsgSize+1;
end

if evalin('base', 'exist(''Process'', ''var'')')
    Process = evalin('base','Process');
    
    % check empty feilds
    processFields = {'classname';'method';'Parameters'};
    
    for i = 1:length(processFields)
        try Process.(processFields{i});
        catch msg
            warningMsg{MsgSize,1} = ['Error! ',msg.message, ' The field name, ', processFields{i} , ', is case sensitive.'];
            errorColor{MsgSize} = 'r';
            MsgSize = MsgSize+1;
        end
        if ~exist('msg','var')
            emptyInd = find(cellfun(@isempty,{Process.(processFields{i})}));
            nonCellInd = find(~cellfun(@iscell,{Process.(processFields{i})})); % Process.Parameters can be an empty "Cell" but not an empty array
            if ~isempty(emptyInd)&&~isequal(i,3)
                warningMsg{MsgSize,1} = ['Error! Process(', num2str(emptyInd),').',processFields{i},' is empty!'];
                errorColor{MsgSize} = 'r';
                MsgSize = MsgSize+1;
            elseif ~isempty(emptyInd)&&isequal(i,3)&&isequal(emptyInd,nonCellInd)
                warningMsg{MsgSize,1} = ['Error! Process(', num2str(emptyInd),').',processFields{i},' needs to be a cell (empty or correct definition)!'];
                errorColor{MsgSize} = 'r';
                MsgSize = MsgSize+1;
            end
        end
    end
    
    [ProcessDiffNum,~] = setdiff(1:1:length(Process),[EventStruct.process]);
    if ~isempty(ProcessDiffNum)
        if length(ProcessDiffNum) < 5
            warningMsg{MsgSize,1} = ['Warning!! Process ', num2str(ProcessDiffNum),' is/are not used in the Event Structure.'];
        else
            warningMsg{MsgSize,1} = ['Warning!! Process ', num2str(ProcessDiffNum(1:5)),' ... and more are not used in the Event Structure.'];
        end
        errorColor{MsgSize} = 'k';
        MsgSize = MsgSize+1;
    end
    
    [ProcessDiffNum,~] = setdiff([EventStruct.process],1:1:length(Process));
    ProcessDiffNum(ProcessDiffNum==0) = [];
    if ~isempty(ProcessDiffNum)
        if length(ProcessDiffNum) < 5
            warningMsg{MsgSize,1} = ['Error!! Process ', num2str(ProcessDiffNum),' is/are used in the Event Sequence but not predefined!'];
        else
            warningMsg{MsgSize,1} = ['Error!! Process ', num2str(ProcessDiffNum(1:5)),' ... and more are used in the Event Sequence but not predefined!'];
        end
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
elseif nnz([EventStruct.process])
    ProcInd = unique(nonzeros([EventStruct.process]));
    warningMsg{MsgSize,1} = ['Error!! Process ', num2str(ProcInd'),' is/are used in the Event Sequence but Process does not exist!'];
    errorColor{MsgSize} = 'r';
    MsgSize = MsgSize+1;
end

if evalin('base', 'exist(''SeqControl'', ''var'')')
    SeqControl = evalin('base','SeqControl');
    
    % SeqControl.command can't be empty
    emptyInd = find(arrayfun(@isempty,{SeqControl.command}));
    if ~isempty(emptyInd)
        warningMsg{MsgSize,1} = ['Error! command in SeqControl ', num2str(emptyInd), ' is empty!'];
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
    [SeqDiffNum,~] = setdiff(1:1:length(SeqControl),[EventStruct.seqControl]);
    if ~isempty(SeqDiffNum)
        warningMsg{MsgSize,1} = ['Warning!! SeqControl ', num2str(SeqDiffNum),' is/are not used in the Event Structure.'];
        errorColor{MsgSize} = 'k';
        MsgSize = MsgSize+1;
    end
    
    [SeqDiffNum,~] = setdiff([EventStruct.seqControl],1:1:length(SeqControl));
    SeqDiffNum(SeqDiffNum==0) = [];
    if ~isempty(SeqDiffNum)
        warningMsg{MsgSize,1} = ['Error!! SeqControl ', num2str(SeqDiffNum),' is/are used in the Event Sequence but not predefined!'];
        errorColor{MsgSize} = 'r';
        MsgSize = MsgSize+1;
    end
    
    % returnToMatlab is required to prevent an endless loop
    % jump back to the first event will add returnToMatlab automatically
    indRTNM = find(strcmp({SeqControl.command}',{'returnToMatlab'}),1);
    indJump = find(strcmp({SeqControl.command}',{'jump'}),1); % jump
    seqNumInEvents = [EventStruct.seqControl];  % seqControl number shown in Events
    if ~ismember(indRTNM,seqNumInEvents) % if returnToMatlab is not used in the Event sequence
        if ismember(indJump,seqNumInEvents) % if jump is used
            if ~isequal(SeqControl(indJump).argument,1) % if the sequence doesn't jump to the first event
                warningMsg{MsgSize,1} = 'Error!! VSX Control will not respond. The only way to break out of the loop is to force quit Matlab.';
                errorColor{MsgSize} = 'r';
                MsgSize = MsgSize+1;
                warningMsg{MsgSize,1} = 'Please add Return to Matlab or jump to the first event in the event sequence';
                errorColor{MsgSize} = 'r';
                MsgSize = MsgSize+1;
            end
        end
    end
    
    % RcvData may only contain few frames for 'one-shot' sequence if the
    % script doesn't jump without sync, triggerIn, or waitForProcessing
    if isempty(indJump) || ~ismember(indJump,seqNumInEvents) % if jump is not used
        indSync = find(strcmp({SeqControl.command}',{'sync'}),1);
        indTriggerIn = find(strcmp({SeqControl.command}',{'triggerIn'}),1);
        waitForProcessing = false;
        if isfield(Resource.Parameters,'waitForProcessing') && ...
                Resource.Parameters.waitForProcessing
            waitForProcessing = true;
        elseif isfield(SeqControl,'condition')
            indWaitProcessing = find(strcmp({SeqControl.condition}',{'waitForProcessing'}),1);
            waitForProcessing = ~isempty(indWaitProcessing) && ismember(indWaitProcessing,seqNumInEvents);
        end
        if (isempty(indSync) || ~ismember(indSync,seqNumInEvents)) && ...
                (isempty(indTriggerIn) || ~ismember(indTriggerIn,seqNumInEvents)) && ...
                ~waitForProcessing
            warningMsg{MsgSize,1} = 'RcvData may only contain a few frames because the script is asynchronous without jump: please add a sync seqControl at the end of Event sequence';
            errorColor{MsgSize} = 'k';
            MsgSize = MsgSize+1;
        end
    end
    
    % TTNA warning estimation
    if evalin('base', 'exist(''Receive'', ''var'')')
        eventOverhead = 4; % 3-5 us, use 4 here
        indTTNA = find(strcmp({SeqControl.command}',{'timeToNextAcq'}));
        for i = 1:length(indTTNA)
            indTTNAevent = find(cellfun(@(x)ismember(indTTNA(i),x),{EventStruct.seqControl})); % get all event numbers with the same TTNA
            % get all seqControl index along with TTNA event
            seqInd = unique([EventStruct(indTTNAevent).seqControl]);
            % TTNA is useless with triggerIn  or sync
            if any(strcmp({SeqControl(seqInd).command},'triggerIn')) || any(strcmp({SeqControl(seqInd).command},'sync'))
                warningMsg{MsgSize,1} = 'A timeToNextAcq is useless when triggerIn or sync is used in the same event: please remove either one, and perhaps add another Event.';
                errorColor{MsgSize} = 'r';
                MsgSize = MsgSize+1;
                break  % once it's detected, jump to line 526
            end
            
            rcvInd = [EventStruct(indTTNAevent).rcv]; % get associated Receive number with the same TTNA
            rcvInd = rcvInd(rcvInd > 0);
            if nnz(rcvInd) > 0
                rcvDepth = unique([Receive(rcvInd).endDepth]); % How many sets of Receive.endDpeth are used?
                for ir = 1:length(rcvDepth)
                    minTTNA = round(rcvDepth(ir)/Trans.frequency)*2+eventOverhead; % in us
                    if SeqControl(indTTNA(i)).argument < minTTNA
                        warningMsg{MsgSize,1} = ['A timeToNextAcq warning might be displayed, try to reduce the receive depth or increase the TTNA value of the SeqControl ', num2str(indTTNA(i))];
                        errorColor{MsgSize} = 'k';
                        MsgSize = MsgSize+1;
                    end
                end
            end
        end
    end
    
    % transferToHost needs to be unique in the event sequences
    indDMA = find(strcmp({SeqControl.command},{'transferToHost'}));
    seqCount = arrayfun(@(x) sum([EventStruct.seqControl]==x),(1:length(SeqControl)));
    for i = 1:length(indDMA)
        if seqCount(indDMA(i))>1
            warningMsg{MsgSize,1} = ['Error!! transferToHost must be unique in the event sequence! SeqControl ', num2str(indDMA(i)),...
                ' cannot be used in multiple events!'];
            errorColor{MsgSize} = 'r';
            MsgSize = MsgSize+1;
        end
    end
    
    % get event number of each DMA event
    indDMAevent = find(cellfun(@any,(cellfun(@(x)ismember(x,indDMA),{EventStruct.seqControl},'UniformOutput',false))));
    
    if alreadyIniz
        if Trans.connType > 0
            activeCG = evalin('base','UTA.activeCG');
            totalActiveCG = sum(activeCG);
            dmaSpeed = min(1.5*sum(activeCG),6.6); % DMA speed is 1.5 MB/ms per channel group (32 channels) and the maximum is 6.6 MB/ms
        else
            totalActiveCG = 4; % simulationOnly - Assume it's a 4 boards system
            dmaSpeed = 6; % simulationOnly - Assume it's 6 GB/s for 128 channels
            % Add notification in the first warning message
            tempMsg = 'Note! This is a simulateOnly script. Use 6 GB/s as the estimated DMA speed.';
            tempColor = 'b';
            warningMsg = [tempMsg;warningMsg];
            errorColor = [tempColor;errorColor];
            MsgSize = MsgSize+1;
        end
        dmaOverheadPerMB = 32; %  per MB for estimation.
    end
    
    showDMAwarning = 0;
    warningDMAind = [];
    sizeOfEachDMA = zeros(length(indDMA),1);
    timeForEachDMA = zeros(length(indDMA),1);
    intervalToNextDMA = zeros(length(indDMA),1);
    
    for i = 1:length(indDMA)
        
        endDMAevent = indDMAevent(i);
        startDMAevent = 1;
        % Get current buffer and frame number
        for n = endDMAevent:-1:1
            if EventStruct(n).rcv > 0
                dmaBufNum = Receive(EventStruct(n).rcv).bufnum;
                dmaFrameNum = Receive(EventStruct(n).rcv).framenum;
                break;
            end
        end
        
        % Search backwards until the event of:
        % 1) another transerToHost
        % 2) different frame of the same buffer
        % 3) different buffer number
        stopByDMA = 1; % The search should stop at another DMA or the first event
        for n = endDMAevent-1:-1:1
            % check seqControl first because the last acquisition usually happends with DMA event
            if EventStruct(n).seqControl > 0
                if any(strcmp({SeqControl(EventStruct(n).seqControl).command},'transferToHost'))
                    startDMAevent = n+1;
                    break
                end
            end
            if EventStruct(n).rcv > 0
                if ~isequal(Receive(EventStruct(n).rcv).bufnum,dmaBufNum) || ...
                        ~isequal(Receive(EventStruct(n).rcv).framenum,dmaFrameNum)
                    startDMAevent = n+1;
                    stopByDMA = 0;
                    break;
                end
            end
        end
        
        dmaRcvInd = [EventStruct(startDMAevent:endDMAevent).rcv];
        dmaRcvInd = dmaRcvInd(dmaRcvInd ~= 0);
        
        % If the search is not stopped by the DMA,
        % check any dummy acquisitions
        if ~stopByDMA && ~isequal(startDMAevent,1)
            previousDMAevent = 1;
            for n = startDMAevent-1:-1:1
                if EventStruct(n).seqControl > 0
                    if any(strcmp({SeqControl(EventStruct(n).seqControl).command},'transferToHost'))
                        previousDMAevent = n;
                    end
                end
            end
            
            previousRcvInd = [EventStruct(previousDMAevent:startDMAevent).rcv];
            previousRcvInd = previousRcvInd(previousRcvInd ~= 0);
            
            if any(~ismember([Receive(previousRcvInd).bufnum],dmaBufNum))
                
            elseif any(~ismember([Receive(previousRcvInd).framenum],dmaFrameNum))
                warningMsg{MsgSize,1} = ['Warning!! Receive ', num2str(previousRcvInd(~ismember([Receive(previousRcvInd).framenum],dmaFrameNum))),...
                    ' will not be transferred by SeqControl(',num2str(indDMA(i)),').'];
                errorColor{MsgSize} = 'k';
                MsgSize = MsgSize+1;
            end
        end
        
        if alreadyIniz && isfield(Receive,'endSample')
            % Calculate data size of one DMA and esitmate the TTNA warning
            perChannelDMA = 2 * (Receive(dmaRcvInd(end)).endSample - Receive(dmaRcvInd(1)).startSample + 1)/2^20; % in MB, not byte
            perFrameDMA = perChannelDMA * 32 * totalActiveCG; % in MB
            timeForDMA = perFrameDMA/dmaSpeed + dmaOverheadPerMB*perFrameDMA/1e3; % in ms
            timeForEachDMA(i) = timeForDMA; % will be used while DMA event is selected
            sizeOfEachDMA(i) = perFrameDMA; % will be used while DMA event is selected
            
            % Calculate the interval between two DMAs because the next DMA
            % won't be initialized before the current one is completed!!
            dmaEnd = 0;
            dmaInterval = 0;
            en = indDMAevent(i);
            seqInd = EventStruct(en).seqControl;
            
            % TTNA warning only happens when TTNA and transferToHost are at
            % the same event
            if any(strcmp({SeqControl(seqInd).command},'timeToNextAcq'))
                % DMA warning estimation doesn't work with synchronized script
                if Resource.Parameters.waitForProcessing == 0
                    while isequal(dmaEnd,0)
                        seqInd = EventStruct(en).seqControl;
                        if ~isequal(seqInd,0)
                            % DMA interval calculation doesn't work with triggerIn or sync command
                            if any(strcmp({SeqControl(seqInd).command},'triggerIn')) || any(strcmp({SeqControl(seqInd).command},'sync'))
                                dmaInterval = 0;
                                break  % break while loop
                            end
                            % jump to assigned event
                            jumpSeqInd = seqInd(strcmp({SeqControl(seqInd).command},'jump'));
                            if jumpSeqInd
                                en = SeqControl(jumpSeqInd).argument;
                                seqInd = EventStruct(en).seqControl;
                                dmaInterval = dmaInterval + eventOverhead; % in us
                            end
                            % check if TTNA is used
                            if any(ismember(seqInd,indTTNA))
                                dmaInterval = dmaInterval + SeqControl(seqInd(ismember(seqInd,indTTNA))).argument; % in us
                            end
                        elseif any(EventStruct(en).tx) || any(EventStruct(en).rcv) % no SeqControl, only Hardware Sequencer needs to be taken into account
                            dmaInterval = dmaInterval + eventOverhead; % in us
                        end
                        en = en+1;
                        
                        % stop if the next DMA or the end of sequence is reached
                        % or it's one frame with jump... (rare case)
                        if en > EventSize || any(ismember(EventStruct(en).seqControl,indDMA)) || isequal(length(indDMA),1)
                            dmaEnd = 1;
                        end
                    end
                end
            end
            intervalToNextDMA(i) = dmaInterval/1e3; %ms
            
            if (dmaInterval > 0) && (dmaInterval/1e3 < timeForDMA)
                showDMAwarning = 1;
                if length(warningDMAind) < 6
                    warningDMAind = [warningDMAind indDMA(i)];
                end
            end
            
        end
    end
    
    if exist('Recon', 'var') % checked from workspace already
        if isfield(Recon,'newFrameTimeout')
            newFrameTimeout = Recon.newFrameTimeout;
        else
            newFrameTimeout = 1000; % default to 1000 ms
        end
        if intervalToNextDMA(1) > newFrameTimeout
            warningMsg{MsgSize,1} = ['A message about timeout waiting for next acquisition frame or next DMA may be displayed. Try to set Recon.newFrameTimeout greater than ',num2str(newFrameTimeout),' (ms).'];
            errorColor{MsgSize} = 'k';
            MsgSize = MsgSize+1;
        end
    end
    
    
    if showDMAwarning
        warningMsg{MsgSize,1} = ['A timeToNextAcq warning might be displayed, check SeqControl  ', num2str(warningDMAind) ,'... for more details'];
        errorColor{MsgSize} = 'k';
        MsgSize = MsgSize+1;
    end
    
elseif nnz([EventStruct.seqControl])
    SeqInd = unique(nonzeros([EventStruct.seqControl]));
    warningMsg{MsgSize,1} = ['Error!! SeqControl ', num2str(SeqInd'),' is/are used in the Event Sequence but SeqControl does not exist!'];
    errorColor{MsgSize} = 'r';
    MsgSize = MsgSize+1;
end

% display error from initialization if it exists
if evalin('base','exist(''errMsg'',''var'')')
    errMsg = evalin('base','errMsg');
    if ~isempty(errMsg)
        warningMsg{MsgSize,1} = ['Error!! ', errMsg.message,];
        errorColor{MsgSize} = 'r';
    end
end

warningMsg(cellfun('isempty',warningMsg)) = [];
errorColor(cellfun('isempty',errorColor)) = [];

%%%% ==== Done with debug ====

%%  ==== regular start ====
%Selection prompt
if isempty(warningMsg)
    showWarning = 0;
    set(dispmessage,'String','Select an event parameter in the table below...');
    dispmessage.Visible = 'on';
else
    showWarning = 1;
    dispmessage.Visible = 'off';
    % Print warning/Error message
    hor = -0.1;
    if length(warningMsg) > 7
        finalMsg = warningMsg(1:6,1);
        finalMsg(7:8,1) = {':'};
        finalMsg(9,1) = {'More possible errors are detected, please see the command window'};
        errorColor(7:9,1) = {'k'};
        disp('---- All detected errors are shown below ----');
        disp(warningMsg);
    else
        finalMsg = warningMsg;
    end
    %     dispWarningMsg = zeros(length(finalMsg),1);
    for L = 1:length(finalMsg)
        dispWarningMsg(L) = uicontrol('Parent',sa,'Style','text',...
            'FontUnits','normalized',...
            'FontSize',messageFont,...
            'Units','normalized',...
            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
            'HorizontalAlignment','center',...
            'BackgroundColor',[.7 .7 .7],...
            'HandleVisibility','off',...
            'FontWeight','bold',...
            'ForegroundColor',errorColor{L},...
            'String',finalMsg{L,1});
        hor = hor+0.05;
    end
end

%Set Column names to max number of events in the sequence
eventColnames = cell(1,EventSize);

for i = 1:EventSize
    eventColnames{1,i} =  strcat(sprintf('Event %s',num2str(i)));
end

%Set Row names for table
eventRownames = {'info','tx','rcv','recon','process','seqControl'};
columnWidth = 80;
ScrnSize = get(groot,'ScreenSize');
tableHeight = 0.22; % Normalized table height

%Reload Workspace Pushbutton
rMode = uicontrol('Style','pushbutton',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',textFont,...
    'Units','normalized',...
    'Position',[.015,.95,.1,.03],...
    'BackgroundColor',[.7 .7 .7],...
    'HandleVisibility','off',...
    'String','Regular mode',...
    'Callback',@selcReload);

if~alreadyIniz
    initStr = 'Click ''initializeOnly mode'' button or type ''EventAnalysisTool i'' in the command window for more information (slower execution)';
    %Reload Workspace Pushbutton with initializedOnly
    uicontrol('Style','pushbutton',...
        'Parent',sa,...
        'FontUnits','normalized',...
        'FontSize',textFont,...
        'Units','normalized',...
        'Position',[.120,.95,.11,.03],...
        'BackgroundColor',[.7 .7 .7],...
        'HandleVisibility','off',...
        'String','InitializeOnly mode',...
        'Callback',{@selcReload,1});
    uicontrol('Style','text',...
        'Parent',sa,...
        'FontUnits','normalized',...
        'FontSize',textFont,...
        'Units','normalized',...
        'Position',[.25,.95,.7,.03],...
        'BackgroundColor',[.7 .7 .7],...
        'HandleVisibility','on',...
        'FontAngle','italic',...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'String',initStr);
else
    rMode.String = 'Reload workspace';
end

% Text for Event info
EventInfoTxt = uicontrol('Style','text',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',titleFont,...
    'Units','normalized',...
    'Position',[.22,.24,.77,.024],...
    'HandleVisibility','off',...
    'FontAngle','italic',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'Tag','EventInfoTxt',...
    'String','EventInfo');

% Jump to
uicontrol('Style','text',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',titleFont,...
    'Units','normalized',...
    'Position',[.01,.24,.04,.024],...
    'BackgroundColor',[.7 .7 .7],...
    'HandleVisibility','off',...
    'FontWeight','bold',...
    'String','Jump to');

JumpMenu = uicontrol('Style','popupmenu',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',titleFont,...
    'Units','normalized',...
    'Position',[.05,.243,.1,.023],...
    'HandleVisibility','off',...
    'FontWeight','bold',...
    'String',{'Event';'tx';'rcv';'recon';'process';'seqControl'});

uicontrol('Style','text',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',titleFont,...
    'Units','normalized',...
    'Position',[.15,.24,.01,.024],...
    'BackgroundColor',[.7 .7 .7],...
    'HandleVisibility','off',...
    'FontWeight','bold',...
    'String','#');

JumpNum = uicontrol('Style','Edit',...
    'Parent',sa,...
    'FontUnits','normalized',...
    'FontSize',titleFont,...
    'Units','normalized',...
    'Position',[.16,.24,.05,.024],...
    'HandleVisibility','off',...
    'FontWeight','bold',...
    'String','1',...
    'Callback',@JumpToEvent);

%Create Event Table
EventTable = uitable('Parent',sa,...
    'ColumnName',eventColnames,...
    'FontSize',tableFontSize,...
    'ColumnWidth', {columnWidth},...
    'HandleVisibility','off',...
    'RowName',eventRownames,...
    'Units','normalized',...
    'Tag','EventTable',...
    'Position',[.01,.015,.98,tableHeight]);

%Set data in table
set(EventTable,'Data',EventSequence);
set(EventTable,'CellSelectionCallback',@EventTable_Callback);
rcvReset = 1; % for the first click of the event table

%% call back functions
%Function Allows Reloading of Parameters in the EventAnalysisTool
    function selcReload(varargin)
        close(sa)
        if isequal(nargin,3)
            evalin('base','EventAnalysisTool i');
        else
            evalin('base','EventAnalysisTool');
        end
    end

AppendInfoNote = [];

% Jump to an event
    function JumpToEvent(varargin)
        % import java library at the top of the function even if it's used
        % at the end of this function
        import java.awt.Robot;
        import java.awt.event.*
        
        % 1) check the max value based on the structures
        fieldInd = get(JumpMenu,'Value');
        minValue = 1;
        switch fieldInd
            case 1 % Event
                maxValue = EventSize;
            case 2 % tx
                maxValue = max([EventStruct.tx]);
            case 3 % rcv
                maxValue = max([EventStruct.rcv]);
            case 4 % recon
                maxValue = max([EventStruct.recon]);
            case 5 % process
                maxValue = max([EventStruct.process]);
            case 6 % seqControl
                maxValue = max([EventStruct.seqControl]);
                % if command is entered, get the associated id here
                if strcmpi('DMA',JumpNum.String)
                    JumpNum.String = 'transferToHost';
                end
                
                if strcmpi('TTNA',JumpNum.String)
                    JumpNum.String = 'timeToNextAcq';
                end
                
                if strcmpi('Return',JumpNum.String)
                    JumpNum.String = 'returnToMatlab';
                end
                
                if isnan(str2double(JumpNum.String))
                    idSet = find(contains({SeqControl.command},JumpNum.String,'IgnoreCase',true));
                    if ~isempty(idSet)
                        if length(idSet) > 1
                            AppendInfoNote = ['   More ', JumpNum.String, ' events can be found in SeqControl ', num2str(idSet)];
                        end
                        JumpNum.String = num2str(idSet(1));
                    else
                        EventInfoTxt.String = [JumpMenu.String{fieldInd},' ',JumpNum.String, ' is not used in the Event Sequence!!'];
                        EventInfoTxt.ForegroundColor = 'r';
                        return
                    end
                end
        end
        
        jumpToInd = min(max(minValue,str2double(JumpNum.String)),maxValue);
        set(JumpNum,'String',num2str(jumpToInd,'%5.0f'));
        
        % 2) find associated event number
        if ~isequal(fieldInd,1)
            EventNum = find(cellfun(@(x)ismember(jumpToInd,x),{EventStruct.(JumpMenu.String{fieldInd})}));
            if ~isempty(EventNum)
                EventNum = EventNum(1);
            else
                EventInfoTxt.String = [JumpMenu.String{fieldInd},' ',JumpNum.String, ' is not used in the Event Sequence!!'];
                EventInfoTxt.ForegroundColor = 'r';
                return
            end
        else
            EventNum = jumpToInd;
        end
        
        % 3) jump and simulate mouse click.
        % NOTE: the getJtable is copied from the findjobj and the copyright
        % is Copyright (c) 2017, Yair Altman
        % All rights reserved.
        
        tableJ = getJtable(EventTable);
        rowHeader = tableJ.RowHeader;
        columnHeader = tableJ.columnHeader;
        rowHeight = rowHeader.getView.getRowHeight;
        
        HorizonScrollBarJ = tableJ.getHorizontalScrollBar;
        
        % Get the pixel location of the top-left corner
        sa.Units = 'pixels'; FigPos = sa.Position; sa.Units = 'normalized';
        
        % Start location of the mouse pointer is located at the center of
        % "Event 1" cell with the unit of pixel
        % Figure Bottom-Left corner + UItable offset + header width + columnWidth/2
        StartPosX = FigPos(1)+EventTable.Position(1)*FigPos(3)+rowHeader.getWidth+columnWidth/2;
        % Screen height - Figure Bottom-Left corner - UItable offset - table height + rowHeight/2
        StartPosY = ScrnSize(1,4)-FigPos(2)-EventTable.Position(2)*FigPos(4)-tableJ.getHeight+columnHeader.getHeight-rowHeight/2+2;
        
        % Move Schroll bar
        barValue = max(0,(EventNum-9) * columnWidth); % put the desired event to the center
        StartEvent = floor(barValue/columnWidth)+1;
        
        MousePosX = floor(StartPosX+(EventNum-StartEvent)*columnWidth);
        MousePosY = floor(StartPosY+(fieldInd)*rowHeight);
        
        % Reach to the end of event
        if EventNum > EventSize - 10
            barValue = HorizonScrollBarJ.getMaximum-HorizonScrollBarJ.getWidth;
            MousePosX = FigPos(1)+15+columnHeader.getX+columnHeader.getWidth-(EventSize-EventNum)*columnWidth-columnWidth/2;
        end
        
        set(tableJ.getHorizontalScrollBar,'Value',barValue);
        pause(0.05)
        tableJ.repaint;
        
        mouse = Robot;
        mouse.mouseMove(MousePosX,MousePosY);
        mouse.mousePress(InputEvent.BUTTON1_MASK);
        mouse.mouseRelease(InputEvent.BUTTON1_MASK);
        
    end

%Set up an event table for the Event sequence
    function EventTable_Callback(~,eventdata)
        %Remove displayed message first run
        if strcmp(dispmessage.Visible,'on')
            dispmessage.Visible = 'off';
        end
        
        if showWarning
            set(dispWarningMsg,'Visible','off'); % set is useful for all handles..
        end
        
        % Set Jump Menu and Number
        eventInd = eventdata.Indices;
        set(JumpMenu,'Value',eventInd(1));
        if ~isequal(eventInd(1),1)
            set(JumpNum,'String',EventSequence{eventInd(1),eventInd(2)});
        else
            set(JumpNum,'String',eventInd(2));
        end
        
        tableMapping = {'info','tx','rcv','recon','process','seqControl'};
        EventParameter = tableMapping{1,eventInd(1,1)};
        
        EventInfo = sprintf('Event %s   -   "%s".',num2str(eventInd(1,2)),cell2mat(EventSequence(1,eventInd(1,2))));
        EventInfoTxt.String = [EventInfo,AppendInfoNote]; AppendInfoNote = [];
        EventInfoTxt.ForegroundColor = 'k';
        EventInfoTxt.Visible = 'on';
        
        switch EventParameter
            case 'info'
                
            case 'tx'
                %If open close previous TW plot
                close(findobj('type','figure','name','Simulated Transmit Channel Plot'))
                
                %Track hold plot for TX case
                if exist('ADhndl','var') == 0 || rcvReset == 1 || reconReset == 1 || processReset == 1 || seqControlReset == 1
                    clf
                    if imageholdReset == 1
                        clear PersistAcqDispArray
                        imageholdReset = 0;
                    end
                    rcvReset = 0;
                end
                
                %Clear no tx message
                if exist('notxhndl','var') == 0
                    if size(notxhndl,1) == 1
                        delete(notxhndl)
                        clear notxhndl
                    end
                end
                
                %Event TX Callback Graphing
                %Get eventdata to parse the
                tableRow = eventInd(1,1)-1;
                tableCol = eventInd(1,2);
                
                %Evaluate in Event Structure
                Event = evalin('base','Event');
                txEvent = Event(tableRow,tableCol).tx;
                assignin('base','txEvent',txEvent);
                
                %Display Message if Event == 0
                if txEvent == 0
                    %Print Trans Title Structure
                    notxhndl = uicontrol('Parent',sa,'Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',messageFont,...
                        'Units','normalized',...
                        'Position',centerPos,...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'HandleVisibility', 'on',...
                        'FontWeight','bold',...
                        'String','No Transmit Specified...');
                    txReset = 1;
                    rcvReset = 1;
                end
                
                %Check for non transmit event
                if txEvent ~= 0
                    if ~ismember(txEvent,1:1:length(TX))
                        clf
                        hor = 0;
                        dispmessage = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
                            'HorizontalAlignment','center',...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'String',['Error!! TX(', num2str(txEvent),') is used in the Event Sequence but not predefined!']);
                        txReset = 1;
                        rcvReset = 1;
                        return
                    end
                    
                    %Checkbox Parameter Structure
                    uicontrol('Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'Position',[checkboxhori-.0125 checkboxvert-.14 .105 .024],...
                        'HorizontalAlignment','center',...
                        'FontWeight','bold',...
                        'String','Programmed TW Plot');
                    
                    %Pushbutton TW waveform
                    uicontrol('Style','pushbutton',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'pos',[checkboxhori+.01 checkboxvert-.18 .06 .025],...
                        'BackgroundColor',[.7 .7 .7],...
                        'String','TW Plot',...
                        'Callback',@selcTW);
                    
                    %Checkbox Parameter Structure
                    uicontrol('Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'Position',[checkboxhori-.01 checkboxvert-.245 .1 .024],...
                        'HorizontalAlignment','center',...
                        'FontWeight','bold',...
                        'String','Simulated TX Plot');
                    
                    %Pushbutton TXPD
                    uicontrol('Style','pushbutton',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'pos',[checkboxhori+.01 checkboxvert-.28 .06 .025],...
                        'BackgroundColor',[.7 .7 .7],...
                        'String','TXPD',...
                        'Callback',@selcTXPD);
                    
                    if evalin('base','exist(''TX_Limits'',''var'');')
                        if evalin('base','~isempty(TX_Limits);') % using profile 5
                            TXEvent = evalin('base','TX_Limits.TXEvent');
                            if ~isempty(TXEvent)
                                if ismember(txEvent,[TXEvent.tx])
                                    % Show Voltage Droop for profile 5
                                    uicontrol('Style','pushbutton',...
                                        'FontUnits','normalized',...
                                        'FontSize',textFont,...
                                        'Units','normalized',...
                                        'pos',[checkboxhori-.01 checkboxvert-.35 .1 .025],...
                                        'BackgroundColor',[.7 .7 .7],...
                                        'String','Show Voltage Droop',...
                                        'Tag','Vdroop',...
                                        'Callback',@plotVoltDroop);
                                else
                                    delete(findobj('Tag','Vdroop'));
                                end
                            end
                        end
                    end
                    
                    %Display Trans structure in uitable
                    colWidth = {150,95};
                    hori1xx = .015;
                    vert1xx = .635;
                    vert1xxOffset = .025;
                    horiSize = .215;
                    vertSize = .3;
                    
                    %Evalin Trans Structure, set row names, create data for table
                    Trans = evalin('base','Trans');
                    RowNames = fieldnames(Trans);
                    StructName = Trans;
                    StructText = 'Trans';
                    EventNum = 1;
                    EventNumText = '';%sprintf('(%d)',EventNum);
                    
                    %Call Function To create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    %Create Trans Table
                    TransTable = uitable('Parent',sa,...
                        'ColumnName',{'Transducer Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx vert1xx horiSize vertSize]);
                    %Set Trans data in table
                    set(TransTable,'Data',TableData);
                    clear TableData RowNames StructName StructText EventNum
                    
                    %Evaluate in TX structures
                    TX = evalin('base','TX');
                    Trans = evalin('base','Trans');
                    
                    % some VDAS parameters will not be displayed
                    txFields = fieldnames(TX);%store all field names in cell array
                    txVDAS = strfind(txFields,'VDAS');
                    TX = rmfield(TX,txFields(cellfun(@(x)isequal(x,1),txVDAS)));
                    TX = rmfield(TX,txFields(cellfun(@(x)strcmp(x,'Numpulses'),txFields)));
                    TX = rmfield(TX,txFields(cellfun(@(x)strcmp(x,'CumOnTime'),txFields)));
                    TX = rmfield(TX,txFields(cellfun(@(x)strcmp(x,'imgProfileMaxHv'),txFields)));
                    RowNames = fieldnames(TX);
                    StructName = TX;
                    StructText = 'TX';
                    EventNum = txEvent;
                    EventNumText = sprintf('(%d)',EventNum);
                    
                    %Call Function To create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    %Create TX Event Table
                    TXTable = uitable('Parent',sa,...
                        'ColumnName',{'TX Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx (vert1xx-0.2)-vert1xxOffset horiSize vertSize-0.1]);
                    
                    %Set TX data in table
                    set(TXTable,'Data',TableData);
                    clear TableData RowNames StructName StructText EventNum
                    
                    
                    %Evalin TX Structure, set row names, create data for table
                    TW = evalin('base','TW');
                    RowNames = fieldnames(TW);
                    StructName = TW;
                    StructText = 'TW';
                    EventNum = TX(txEvent).waveform;
                    EventNumText = sprintf('(%d)',EventNum);
                    
                    %Call Function To create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    %Create TW Event Table
                    TWTable = uitable('Parent',sa,...
                        'ColumnName',{'TW Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx (vert1xx-0.2-0.1)-2*vert1xxOffset horiSize vertSize-0.2]);
                    
                    %Set TW data in table
                    set(TWTable,'Data',TableData);
                    clear TableData RowNames StructName StructText EventNum
                    
                    % Plot
                    
                    %Variables to reformat gui
                    %                     ApOffsetX = 0;
                    %                     ApOffsetY = 0;
                    ApMoveX = 0;
                    ApMoveY = 0;
                    %Plot Aperture graphically
                    if isfield(Trans,'HVMux') == 1
                        %Set TX/Delay plot size
                        %                         ApOffsetX = -.2;
                        %                         ApOffsetY = 0;
                        ApMoveX = 0;
                        ApMoveY = -.1;
                        %Display Aperture Data
                        if isfield(TX(1,txEvent),'aperture') == 1
                            Aperture = TX(1,txEvent).aperture;  %Import Aperture Data
                        else
                            Aperture = 1;
                        end
                        if isfield(Trans.HVMux,'ApertureES')
                            ActEle = Trans.HVMux.ApertureES(:,Aperture);
                        else
                            ActEle = Trans.HVMux.Aperture(:,Aperture);
                        end
                        ActEle(ActEle ~= 0) = 1;
                        ActEle(ActEle == 0) = NaN;
                        MuxPlot = axes('units','normalized','position',[plot1hori+.17 plot1vert-.045 .4 .1]);
                        plot(ActEle,'r.');
                        axis([1 Trans.numelements .5 1.5]);
                        ylabel('On/Off','fontsize',12);
                        xlabel(['Aperture 1:' num2str(Trans.numelements)],'fontsize',12);
                        set(gca, 'Yticklabel', []);
                        set(gca,'FontSize',12);
                        title('Active Aperture Position')
                        %                         if ismember(Trans.type,[0,1]),text(Aperture,1.2,(['\downarrow  ' num2str(Aperture)]),'FontSize',12); end
                    else
                        ActEle = ones(1,Trans.numelements);
                    end
                    
                    %Plot Apod
                    ActiveTX = size(TX(1,txEvent).Apod);%Can only be a 1 x #ofTX
                    %Plot Delay
                    ActiveDelay = size(TX(1,txEvent).Delay); %Can only be a 1 x #ofTX
                    
                    %Find max delay value
                    NumberOfTX = size(TX,2);
                    for j = 1:NumberOfTX
                        if j == 1
                            delayval = ceil(max(TX(1,j).Delay));
                        else
                            if delayval < ceil(max(TX(1,j).Delay))
                                delayval = ceil(max(TX(1,j).Delay));
                            end
                        end
                    end
                    
                    switch Trans.type
                        case {0,1}
                            
                            %Check for var and handle
                            if isempty(ADhndl)||~ishandle(ADhndl)
                                %Plot both on same axis
                                ADhndl = axes('Parent',sa,'units','normalized',...
                                    'NextPlot','replacechildren','Xlim',[1 ActiveTX(1,2)],...
                                    'position',[plot2hori+ApMoveX plot2vert+ApMoveY .55 .525+ApMoveY]);
                            end
                            
                            %Logic to turn on and off apod
                            x = 1:ActiveTX(1,2);
                            y1 = TX(1,txEvent).Apod(1:ActiveTX(1,2));
                            y1(y1 == 0) = NaN;
                            y2 = TX(1,txEvent).Delay(1:ActiveDelay(1,2));
                            y2(y2 == 0) = NaN;
                            [AX, H1, H2] = plotyy(ADhndl,x,y1,x,y2,'plot');
                            
                            %Fix Apod axis
                            set(AX(1),'ylim',[-1.5 1.5],'Ytick',[-1 -0.5 0 0.5 1]);
                            grid on;
                            %                 set(AX(1), 'Ytick',0:.2:1);
                            
                            if delayval ~= 0
                                %                     Fix Delay axis
                                set(AX(2), 'ylim',[0 delayval],'Ytick',[0 delayval]);
                                if (delayval/10) > 1
                                    set(AX(2), 'Ytick',0:10:delayval);
                                end
                            end
                            
                            %Set plot properties
                            set(H1,'LineStyle','-','Marker','.','MarkerSize',14)
                            set(H2,'LineStyle','-','Marker','.','MarkerSize',10)
                            get(AX(1),'Nextplot');
                            get(AX(2),'Nextplot');
                            set(AX,'xlim',[1 ActiveTX(1,2)]);
                            set(H2,'color','r');
                            set(AX(2), 'YColor', 'r')
                            set(get(AX(1),'Ylabel'),'String','TX Apod Value','fontsize',12)
                            set(get(AX(2),'Ylabel'),'String','TX Delay Value','color','r','fontsize',12)
                            xlabel(['Transmit Channels 1:',num2str(ActiveTX(1,2))],'fontsize',12)
                            title('TX Apod and TX Delay','fontsize',12)
                            set(AX,'FontSize',12);
                            txReset = 1;
                            
                        case {2,4}
                            if isempty(Apodhndl)||~ishandle(Apodhndl)
                                %Spearate Apod and Delay into two axes
                                Apodhndl = axes('Parent',sa,'units','normalized','NextPlot','replacechildren',...
                                    'position',[plot2hori+ApMoveX-0.02 plot2vert+ApMoveY .55/2 .525+ApMoveY]);
                                
                                Delayhndl = axes('Parent',sa,'units','normalized','NextPlot','replacechildren',...
                                    'position',[plot2hori+ApMoveX+0.3 plot2vert+ApMoveY .55/2 .525+ApMoveY]);
                            end
                            
                            % subplot TX.Apod
                            if isequal(length(TX(txEvent).Apod),Trans.numelements) % dynamic HVMux
                                plot3(Apodhndl,Trans.ElementPos((ActEle == 1),1), Trans.ElementPos((ActEle == 1),2), TX(txEvent).Apod(ActEle==1), 'ro');
                            else
                                plot3(Apodhndl,Trans.ElementPos((ActEle == 1),1), Trans.ElementPos((ActEle == 1),2), TX(txEvent).Apod, 'ro');
                            end
                            set(Apodhndl,'XLim',[min(Trans.ElementPos(:,1)),max(Trans.ElementPos(:,1))]);
                            set(Apodhndl,'YLim',[min(Trans.ElementPos(:,2)),max(Trans.ElementPos(:,2))]);
                            set(Apodhndl,'ZLim',[-1, 1.5]);title(Apodhndl,'Apod');grid(Apodhndl,'on');
                            %                             xlabel(Apodhndl,'Element Position (wls)'),ylabel(Apodhndl,'Element Position (wls)'),
                            zlabel(Apodhndl,'TX.Apod');view(Apodhndl,[15,60]);
                            
                            % subplot TX.Delay
                            if isequal(length(TX(txEvent).Apod),Trans.numelements)
                                plot3(Delayhndl,Trans.ElementPos((ActEle == 1),1), Trans.ElementPos((ActEle == 1),2), TX(txEvent).Delay(ActEle==1), 'b+');
                            else
                                plot3(Delayhndl,Trans.ElementPos((ActEle == 1),1), Trans.ElementPos((ActEle == 1),2), TX(txEvent).Delay, 'b+');
                            end
                            set(Delayhndl,'XLim',[min(Trans.ElementPos(:,1)),max(Trans.ElementPos(:,1))]);
                            set(Delayhndl,'YLim',[min(Trans.ElementPos(:,2)),max(Trans.ElementPos(:,2))]);
                            if ~isequal(max([TX.Delay]),0)
                                set(Delayhndl,'ZLim',[0,max([TX.Delay])]);
                            else
                                set(Delayhndl,'ZLim',[0,1]);
                            end
                            title(Delayhndl,'Delay');grid(Delayhndl,'on');
                            
                            %                             xlabel(Delayhndl,'Element Position (wls)'),ylabel(Delayhndl,'Element Position (wls)'),
                            zlabel(Delayhndl,'TX.Delay (wavelength)');view(Delayhndl,[15,60]);
                            
                            h = rotate3d;
                            h.Enable = 'on';
                            if exist('MuxPlot','var')
                                setAllowAxesRotate(h,MuxPlot,false);
                            end
                            txReset = 1;
                            
                        otherwise
                            if isempty(ADhndl)||~ishandle(ADhndl)                                
                                ADhndl = uicontrol('Parent',sa,'Style','text',...
                                    'FontUnits','normalized',...
                                    'FontSize',messageFont,...
                                    'Units','normalized',...
                                    'Position',centerPos+[hori1xx+checkboxhori,0,0,0],...
                                    'HorizontalAlignment','center',...
                                    'BackgroundColor',[.7 .7 .7],...
                                    'HandleVisibility', 'on',...
                                    'FontWeight','bold',...
                                    'String','Not supported transducer type');
                                txReset = 1;
                            end
                    end
                end
                
            case 'rcv'
                %If open close previous TW plot
                close(findobj('type','figure','name','Simulated Transmit Channel Plot'))
                
                %Get eventdata to parse the
                tableRow = eventInd(1,1)-2;
                tableCol = eventInd(1,2);
                
                %Position variables for receive
                plot1rcvhori = .365;
                plot2rcvhori = .685;
                plot1rcvvert = .825;
                plot2rcvvert = .3;
                
                %Evaluate in Event Structure
                Event = evalin('base','Event');
                rcvEvent = Event(tableRow,tableCol).rcv;
                
                %Evaluate in Resource Structure
                Resource = evalin('base','Resource');
                
                %Keep screen if event is set to 0
                if rcvEvent ~= 0
                    %Track hold plot for TX case
                    if exist('ADhndl','var') == 0 || txReset == 1 || reconReset == 1 || processReset == 1 || seqControlReset == 1 || rcvReset == 1
                        clf
                        rcvReset = 0;
                        txReset = 0;
                    end
                end
                
                %Display Message if Event == 0
                if rcvEvent == 0
                    clf
                    %Print Trans Title Structure
                    uicontrol('Parent',sa,'Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',messageFont,...
                        'Units','normalized',...
                        'Position',centerPos,...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'HandleVisibility', 'on',...
                        'FontWeight','bold',...
                        'String','No Receive Specified...');
                    rcvReset = 1;
                end
                
                if rcvEvent ~= 0
                    errorMsg = [];
                    if evalin('base', 'exist(''Receive'', ''var'')')
                        if ~ismember(rcvEvent,1:1:length(Receive))
                            errorMsg = ['Error!! Receive(', num2str(rcvEvent),') is used in the Event Sequence but not predefined!'];
                        end
                    else
                        errorMsg = ['Error!! Receive(', num2str(rcvEvent),') is used in the Event Sequence but Receive does not exist!'];
                    end
                    
                    if ~isempty(errorMsg)
                        clf
                        hor = 0;
                        dispmessage = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
                            'HorizontalAlignment','center',...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'String',errorMsg);
                        rcvReset = 1;
                        return
                    end                    
                    
                    %Evalute in Structures
                    Receive = evalin('base','Receive');
                    Trans = evalin('base','Trans');                   
                                        
                    %Display Trans structure in uitable
                    colWidth = {150,95};
                    hori1xx = .015;
                    vert1xx = .81;
                    vert1xxOffset = .025;
                    horiSize = .215;
                    vertSize = .125;
                    
                    %Evalin Trans Structure, set row names, create data for table
                    RowNames = fieldnames(Trans);
                    StructName = Trans;
                    StructText = 'Trans';
                    EventNum = 1;
                    EventNumText = '';%sprintf('(%d)',EventNum);
                    
                    %                     Call Function To create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    %Create Event Table
                    TransTable = uitable('Parent',sa,...
                        'ColumnName',{'Transducer Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx vert1xx horiSize vertSize]);
                    
                    %Set Trans data in table
                    set(TransTable,'Data',TableData);
                    clear TableData RowNames StructName StructText EventNum
                    
                    % All VDAS parameters won'tbe displayed
                    rcvFields = fieldnames(Receive);%store all field names in cell array
                    rcvVDAS = strfind(rcvFields,'VDAS');
                    Receive = rmfield(Receive,rcvFields(cellfun(@(x)isequal(x,1),rcvVDAS)));
                    RowNames = fieldnames(Receive);
                    StructName = Receive;
                    StructText = 'Receive';
                    EventNum = rcvEvent;
                    EventNumText = sprintf('(%d)',EventNum);
                    
                    %Call function to create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    vertSize = 0.5;
                    %Create Receive Table
                    ReceiveTable = uitable('Parent',sa,...
                        'ColumnName',{'Receive Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx (vert1xx-vertSize)-vert1xxOffset horiSize vertSize]);
                    
                    %Set Trans data in table
                    set(ReceiveTable,'Data',TableData);
                    clear TableData RowNames StructName StructText EventNum
                    
                    %Pushbutton Filters
                    filterButton = uicontrol('Style','pushbutton',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'Position',[hori1xx+horiSize+0.025, (vert1xx-vertSize)-vert1xxOffset+0.01 , 0.07, 0.03],...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'String','Show Filters',...
                        'Callback',{@filterTool,rcvEvent});
                    
                    %                     ar = annotation('arrow');
                    %                     ar.LineWidth = 2;
                    %                     ar.Position = [0.232 0.31 0.02 0];
                    
                    %Reset Rcv Plot
                    HoldPlotTxt = uicontrol('Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',textFont,...
                        'Units','normalized',...
                        'Position',[checkboxhori-.025 vert1x-.225 .08 .025],...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'FontWeight','bold',...
                        'String','Rcv Buffer Reset');
                    
                    %Pushbutton Reset Rcv Plot
                    %Check for var and handle
                    if isempty(imageholdhndl) == 1
                        imagehold = 1;
                        imageholdhndl = uicontrol('Style','Check','String','Hold Plot','Units','normalized',...
                            'Position',[checkboxhori-.015 vert1x-.255 .055 .02],'Value',imagehold);
                        set(imageholdhndl,'Callback',@rcvbuffplotreset);
                    else
                        imageholdhndl = uicontrol('Style','Check','String','Hold Plot','Units','normalized',...
                            'Position',[checkboxhori-.015 vert1x-.255 .055 .02],'Value',imagehold);
                        set(imageholdhndl,'Callback',@rcvbuffplotreset);
                    end
                    
                    %Variables to reformat gui
                    ApOffsetX = 0;
                    %                     ApOffsetY = 0;
                    %                     ApMoveX = 0;
                    %                     ApMoveY = 0;
                    %Plot Aperture graphically
                    if isfield(Trans,'HVMux') == 1
                        %Set TX/Delay plot size
                        ApOffsetX = -.325;
                        %                         ApOffsetY = -80;
                        %                         ApMoveX = -100;
                        %                         ApMoveY = -70;
                        %Display Aperture Data
                        if isfield(Receive(1,rcvEvent),'aperture') == 1
                            Aperture = Receive(1,rcvEvent).aperture;  %Import Aperture Data
                        else
                            Aperture = 1;
                        end
                        if isfield(Trans.HVMux,'ApertureES')
                            ActEle = Trans.HVMux.ApertureES(:,Aperture);
                        else
                            ActEle = Trans.HVMux.Aperture(:,Aperture);
                        end
                        ActEle(ActEle ~= 0) = 1;
                        ActEle(ActEle == 0) = NaN;
                        axes('units','normalized','position',[plot2rcvhori plot1rcvvert .601+ApOffsetX .125]);
                        aper = plot(ActEle,'r.');
                        axis([1 Trans.numelements 0 2]);
                        set(gca,'Xtick',[1 Trans.numelements]);
                        set(gca,'Ytick',1);
                        set(gca,'Yticklabel',{'On/Off'})
                        set(aper,'Color','red','LineWidth',2);
                        xlabh1 = get(gca,'XLabel');
                        set(xlabh1,'Position',get(xlabh1,'Position') + [0 .075 0])
                        xlabel(['Aperture 1:' num2str(Trans.numelements)],'fontsize',11)
                        title('Active Aperture Position','fontsize',11)
                        %                         if ismember(Trans.type,[0,1]),text(Aperture,1.2,(['\downarrow  ' num2str(Aperture)]),'FontSize',12); end
                    end
                    
                    %Plot apod every time
                    x2 = [Receive(1,rcvEvent).Apod];
                    x2(x2 == 0) = NaN;
                    
                    %Plot setting/properties
                    axes('units','normalized','position',[plot2rcvhori-.322 plot1rcvvert .601+ApOffsetX .125]);
                    plot(x2,'b.');
                    xlabel(['Apod 1:',num2str(length(x2))],'units','normalized','fontsize',11)
                    xlabh1 = get(gca,'XLabel');
                    set(xlabh1,'Position',get(xlabh1,'Position') + [0 .075 0])
                    xlim([1 length(x2)])
                    set(gca,'Xtick',[1 Trans.numelements]);
                    title('Apod','units','normalized','fontsize',11)
                    
                    %Display RcvSamples in Acq Event
                    
                    numFrames = Resource.RcvBuffer(Receive(rcvEvent).bufnum).numFrames;
                    RcvInd = arrayfun(@(x) isequal(x.bufnum,Receive(rcvEvent).bufnum),Receive);
                    numAcqs = max(arrayfun(@(x) x.acqNum, Receive(RcvInd)));
                    
                    AcqArray = zeros(numAcqs,numFrames);
                    AcqDispArray = zeros(numAcqs,numFrames);
                    AcqArray(Receive(1,rcvEvent).acqNum,Receive(1,rcvEvent).framenum) = Receive(1,rcvEvent).acqNum;
                    AcqDispArray(Receive(1,rcvEvent).acqNum,Receive(1,rcvEvent).framenum) = 1;
                    
                    %Set up persitence on image plot
                    if imagehold == 1
                        %Assignin AcqDispArray
                        if isempty(imageholdReset) == 1 || imageholdReset == 0 || ~isequal(size(PersistAcqDispArray),size(AcqDispArray))
                            %Assign Acq the first time
                            PersistAcqDispArray = AcqDispArray;
                            imageholdReset = 1;
                        else
                            %OR new info to array
                            PersistAcqDispArray = PersistAcqDispArray | AcqDispArray;
                            %Set imagehold bit
                            imageholdReset = 1;
                            %Assign persistent data to image plot variable
                            AcqDispArray = PersistAcqDispArray;
                        end
                    end
                    
                    %Image must be displayed two times in order to line up grid
                    %axes
                    rcvImg = axes('Parent',sa,'units','normalized','position',[plot1rcvhori plot2rcvvert .6 .45]);
                    imagesc(AcqDispArray);
                    set(rcvImg, 'layer', 'top');
                    if numFrames <= 30
                        set(rcvImg, 'Xtick',0:1:numFrames);
                    end
                    if size(AcqArray(:,1),1) <40
                        set(rcvImg, 'Ytick',0:1:size(AcqArray(:,1),1));
                    end
                    set(rcvImg, 'xcolor', 'w', 'ycolor', 'w');
                    set(rcvImg, 'YDir', 'reverse')
                    hold on
                    %Set the grid at .5 between the number ticklabel
                    h = axes('Parent',sa,'units','normalized','position',[plot1rcvhori plot2rcvvert .6 .45]);
                    imagesc(AcqDispArray)
                    if numFrames <= 30
                        set(h, 'Xtick',.5:1:numFrames);
                    end
                    if size(AcqArray(:,1),1) < 40
                        set(h, 'Ytick',.5:1:size(AcqArray(:,1),1));
                    end
                    set(h, 'yGrid', 'on');
                    set(h, 'xGrid', 'on');
                    set(h, 'Yticklabel', []);
                    set(h, 'Xticklabel', []);
                    set(h, 'xcolor', 'w', 'ycolor', 'w');
                    title(['Receive Buffer ',num2str(Receive(rcvEvent).bufnum)],'fontsize',14,'color','w')
                    %Set the position of xlabel
                    xlabel('Frame Number','fontsize',14,'units','normalized')
                    xlabh = get(gca,'XLabel');
                    set(xlabh,'Position',get(xlabh,'Position') - [0 .01 0])
                    %Set the position of y label
                    ylabel('Acquisition Number','fontsize',14,'units','normalized')
                    ylabh = get(gca,'YLabel');
                    set(ylabh,'Position',get(ylabh,'Position') - [.01 0 0])
                    hold off
                    assignin('base','rcvEvent',rcvEvent);
                    rcvReset = 1;
                    
                    clear AcqDispArray
                end
                
            case 'recon'
                
                if evalin('base','exist(''Recon'',''var'')')
                    Recon = evalin('base','Recon');
                    ReconInfo = evalin('base','ReconInfo');
                else
                    reconEvent = 0;
                end
                
                %Get eventdata to parse the
                tableRow = eventInd(1,1)-3;
                tableCol = eventInd(1,2);
                
                %Evaluate in Event Structure
                Event = evalin('base','Event');
                reconEvent = Event(tableRow,tableCol).recon;
                
                %Clear Previous Events
                if exist('ADhndl','var') == 0 || txReset == 1 || rcvReset == 1 || reconReset == 1 || processReset == 1 || seqControlReset == 1
                    clf
                    rcvReset = 0;
                    txReset = 0;
                    reconReset = 0;
                    processReset = 0;
                    seqControlReset = 0;
                end
                
                %Display Message if Event == 0
                if reconEvent == 0
                    %Print Trans Title Structure
                    uicontrol('Parent',sa,'Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',messageFont,...
                        'Units','normalized',...
                        'Position',centerPos,...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'HandleVisibility', 'on',...
                        'FontWeight','bold',...
                        'String','No Recon Specified...');
                    reconReset = 1;
                end
                
                %                 hori = 0;
                
                if reconEvent ~= 0
                    errorMsg = [];
                    if evalin('base', 'exist(''Recon'', ''var'')')
                        if ~ismember(reconEvent,1:1:length(Recon))
                            errorMsg = ['Error!! Recon(', num2str(reconEvent),') is used in the Event Sequence but not predefined!'];
                        end
                    else
                        errorMsg = ['Error!! Recon(', num2str(reconEvent),') is used in the Event Sequence but Recon does not exist!'];
                    end
                    
                    if ~isempty(errorMsg)
                        clf
                        hor = 0;
                        dispmessage = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
                            'HorizontalAlignment','center',...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'String',errorMsg);
                        rcvReset = 1;
                        return
                    end
                    
                    %%%% ==== Recon Structure Display ====
                    colWidth = {160,100};
                    hori1xx = .015;
                    vert1xx = .735;
                    vert1xxOffset = .025;
                    horiSize = .215;
                    vertSize = .2;
                    
                    %Evalin Recon Structure, set row names, create data for table
                    Recon = evalin('base','Recon');
                    RowNames = fieldnames(Recon);
                    StructName = Recon;
                    StructText = 'Recon';
                    EventNum = reconEvent;
                    
                    %Logic to address multiple entries in Event
                    for ii = 1:size(EventNum,2)
                        EventNumSelect = EventNum(ii);
                        EventNumText = sprintf('(%d)',EventNumSelect);
                        %Call function to create data in table format
                        TableData = dataTableCreate(RowNames,StructName,StructText,EventNumSelect,EventNumText);
                        
                        %Create Event Table
                        ReconTable = uitable('Parent',sa,...
                            'ColumnName',{'Reconstruction Object';'Entry/Value'},...
                            'FontSize',tableFontSize,...
                            'ColumnWidth', colWidth,...
                            'RowName',[],...
                            'Units','normalized',...
                            'Position',[hori1xx (vert1xx-(ii-1)*vertSize)-(ii-1)*vert1xxOffset horiSize vertSize]);
                        
                        %Set Trans data in table
                        set(ReconTable,'Data',TableData);
                        clear TableData EventNumSelect EventNumText
                    end
                    clear RowNames StructName StructText EventNum
                    
                    %%%% ==== ReconInfo and In/Out Structure Display ====
                    
                    %Get size of event can be more than one recon per event
                    reconEventSize = size(reconEvent,2);
                    %Check for multiple recons in one event set to the
                    %first event
                    if reconEventSize > 1
                        reconVal = 1;
                        tableReconEvent = reconEvent(1,reconVal);
                    else
                        tableReconEvent = reconEvent;
                    end
                    
                    % === 1st uitable shows Source and Output ===
                    if ~isfield(Recon(1,tableReconEvent),'rcvBufFrame')
                        Recon(1,tableReconEvent).rcvBufFrame = [];
                    end
                    
                    InOutData = {[blanks(8),'RcvBuffer'],[blanks(8),'InterBuffer'],[blanks(8),'ImageBuffer'];...
                        '[ ]','[ ]','[ ]';...
                        [blanks(15),num2str(Recon(1,tableReconEvent).rcvBufFrame)],'[ ]','[ ]'};
                    if isfield(Recon(1,tableReconEvent),'IntBufDest') && ~isempty(Recon(1,tableReconEvent).IntBufDest)
                        InOutData(2:3,2) = {[blanks(15),num2str(Recon(1,tableReconEvent).IntBufDest(1))];...
                            [blanks(15),num2str(Recon(1,tableReconEvent).IntBufDest(2))]};
                    end
                    if isfield(Recon(1,tableReconEvent),'ImgBufDest') && ~isempty(Recon(1,tableReconEvent).ImgBufDest)
                        InOutData(2:3,3) = {[blanks(15),num2str(Recon(1,tableReconEvent).ImgBufDest(1))];...
                            [blanks(15),num2str(Recon(1,tableReconEvent).ImgBufDest(2))]};
                    end
                    InOutRowName = {[],'buffer num','frame num'};
                    InOutColName = {'Source','Destination','Destination'};
                    
                    % === 2nd uitable shows ReconInfo and associated Receive Structures ===
                    
                    %Get the ReconInfo Field names, only four fields are
                    %required for eventanalysis
                    ReconInfoNew = struct('mode',{ReconInfo.mode},...
                        'txnum',{ReconInfo.txnum},...
                        'rcvnum',{ReconInfo.rcvnum},...
                        'pagenum',[],...
                        'regionnum',{ReconInfo.regionnum});
                    
                    if isfield(ReconInfo,'pagenum')
                        [ReconInfoNew(:).pagenum] = deal(ReconInfo.pagenum);
                    end
                    
                    ReconInfoColnames = {'mode','txnum','rcvnum','bufnum','framenum','acqNum','pagenum','regionnum'};
                    
                    ReconInfo = ReconInfoNew;
                    Receive = evalin('base','Receive');
                    
                    %Get the number of how many field names
                    sizeReconInfoCol = size(ReconInfoColnames,2);
                    ReconInfoSize = length(Recon(1,tableReconEvent).RINums);
                    
                    %Preallocate Event Cell array
                    ReconInfos = cell(ReconInfoSize,sizeReconInfoCol); %Get the total number of events in the sequence
                    
                    % Preallocate the space for warning message
                    if isempty(ReconInfoWarning) || ~ishandle(ReconInfoWarning)
                        %                         ReconInfoWarning = [];
                        ReconInfoWarning = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',[hori1x+.3,vert1x-0.625,.6,.04],...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'Visible','off');
                    else
                        ReconInfoWarning.Visible = 'off';
                    end
                    
                    for n = 1:sizeReconInfoCol
                        %Set field name for reconinfo
                        reconInfoFieldName = ReconInfoColnames{1,n};
                        
                        %Use switch statement to input data in the reconinfo structure
                        switch reconInfoFieldName
                            case 'mode'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).mode);
                                end
                                
                            case 'txnum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).txnum);
                                end
                                
                            case 'rcvnum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).rcvnum);
                                end
                                
                            case 'bufnum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    rcvNum = ReconInfo(RIindex).rcvnum;
                                    ReconInfos{k,n}  = num2str(Receive(rcvNum).bufnum);
                                end
                                InOutData{2,1} = [blanks(15),num2str(Receive(rcvNum).bufnum)];
                                % If the RcvBuffer among all ReconInfos is not the consistent, show Warning
                                if ReconInfoSize > 1
                                    RcvNum = ReconInfo(Recon(1,tableReconEvent).RINums(1)).rcvnum:ReconInfo(Recon(1,tableReconEvent).RINums(ReconInfoSize)).rcvnum;
                                    [BufNum,~,ind] = unique([Receive(RcvNum).bufnum]);
                                    in = find(ind~=mode(ind));
                                    if numel(BufNum) > 1
                                        for rowNum = in(1):in(end)
                                            ReconInfos(rowNum,n) = {horzcat('<html><table border=0 width=65 bgcolor="red"><TR><TD>',ReconInfos{rowNum,n},'</TD></TR> </table></html>')};
                                        end
                                        %Print warning message
                                        set(ReconInfoWarning,'String','Error!! The bufnum is not consistent!','Visible','on');
                                        InOutData(2,1) = {horzcat('<html><table border=0 width=115 bgcolor="red"><TR><TD>Inconsistent!</TD></TR> </table></html>')};
                                    end
                                end
                                
                            case 'framenum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    rcvNum = ReconInfo(RIindex).rcvnum;
                                    ReconInfos{k,n}  = num2str(Receive(rcvNum).framenum);
                                end
                                
                                RcvNum = ReconInfo(Recon(1,tableReconEvent).RINums(1)).rcvnum:ReconInfo(Recon(1,tableReconEvent).RINums(ReconInfoSize)).rcvnum;
                                
                                
                                if numel(unique([Receive(RcvNum).endDepth])) > 1
                                    depthWarning = copy(ReconInfoWarning);
                                    depthWarning.Position(2) = vert1x-0.66;
                                    set(depthWarning,'Parent',sa,'String',['Error!! Receive.endDepth across Receive ', num2str(RcvNum(1)),' to ', ...
                                        num2str(RcvNum(end)),' must be identical for Recon ',num2str(tableReconEvent)],'ForegroundColor','r','Visible','on');
                                end
                                
                                % All frame number in the ReconInfo should be consistent before beling replaced by the Recon.rcvBufFrame
                                [frameNum,~,ind] = unique([Receive(RcvNum).framenum]);
                                in = find(ind~=mode(ind));
                                if numel(frameNum) > 1
                                    for rowNum = in(1):in(end)
                                        ReconInfos(rowNum,n) = {horzcat('<html><table border=0 width=65 bgcolor="red"><TR><TD>',ReconInfos{rowNum,n},'</TD></TR> </table></html>')};
                                    end
                                    set(ReconInfoWarning,'String','Error!! The framenum is not consistent!','Visible','on');
                                    %                                     InOutData(3,1) = {horzcat('<html><table border=0 width=115 bgcolor="red"><TR><TD>Inconsistent!</TD></TR> </table></html>')};
                                else
                                    ReconFrame = Recon(1,tableReconEvent).rcvBufFrame;
                                    if ~isempty(ReconFrame)
                                        % if Recon.rcvBufFrame is not empty, applys it to the frame no. in Reconfo and display warning
                                        ReconInfos(1:ReconInfoSize,n) = {num2str(ReconFrame)};
                                        set(ReconInfoWarning,'String','Note! Recon.rcvBufFrame over-rides the frame number specified by the ReconInfo','ForegroundColor','k','Visible','on');
                                    end
                                end
                                
                            case 'acqNum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    rcvNum = ReconInfo(RIindex).rcvnum;
                                    ReconInfos{k,n}  = num2str(Receive(rcvNum).acqNum);
                                end
                                
                            case 'pagenum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).pagenum);
                                end
                                
                            case 'regionnum'
                                for k = 1:ReconInfoSize
                                    RIindex = Recon(1,tableReconEvent).RINums(k);
                                    ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).regionnum);
                                end
                        end
                        
                        %                         if ~isequal(n,1) && ReconInfoSize < 200 % don't shift the number to center with too many ReconInfos
                        %                             dat = ReconInfos(:,n);
                        %                             if cellfun(@(x) ~contains(x,'html'), dat) % center aliment only with numetric array
                        %                                 ReconInfos(:,n) = cellfun(@(x) sprintf('%*s',max(cellfun(@length,dat)+7),x),dat,'UniformOutput',false);
                        %                             end
                        %                         end
                    end
                    
                    % === generate uitables ===
                    
                    %Create InOutTable Table
                    InOutTable = uitable('Parent',sa,...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', {110,115,115,115},...
                        'Units','normalized',...
                        'Position',[hori1x+.45 vert1x-.15 .36 .05+.025*3]);
                    
                    %Set Data in RINums Table
                    set(InOutTable,'Data',InOutData);
                    set(InOutTable,'RowName',InOutRowName);
                    set(InOutTable,'ColumnName',InOutColName);
                    
                    % Input and Output for Image Reconstruction
                    uicontrol('Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',titleFont,...
                        'Units','normalized',...
                        'Position',[hori1x+.45,vert1x-0.023,.36,.02],...
                        'FontAngle','italic',...
                        'HorizontalAlignment','center',...
                        'FontWeight','bold',...
                        'String','Source and Destination for Image Reconstruction');
                    
                    %Create ReconInfoTxt Table
                    ReconInfoTxtTable = uitable('Parent',sa,...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', {140+40,65*5,65*2},...
                        'Units','normalized',...
                        'Position',[hori1x+.35 vert1x-.25 .54 .03]);
                    set(ReconInfoTxtTable,'ColumnName',{'ReconInfo Structures','Source','Destination'});
                    
                    %Create ReconInfo Table
                    ReconInfoTable = uitable('Parent',sa,...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', {140,65,65,65,65,65,65,65},...
                        'Units','normalized',...
                        'Position',[hori1x+.35 vert1x-.55 .553 .05+.025*10]);
                    
                    ReconInfoRowNames = cell(ReconInfoSize,1);
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfoRowNames{k} = ['RI ', num2str(RIindex, '%03.0f')];
                    end
                    
                    % Check if ReconInfo is used in another Recon
                    thisReconInfo = [Recon(1,tableReconEvent).RINums];
                    ReconInd = 1:length(Recon);
                    ReconInd(tableReconEvent) = [];
                    for rn = ReconInd
                        repeatRI = ismember(thisReconInfo,[Recon(rn).RINums]);
                        if any(repeatRI)
                            set(ReconInfoWarning,'String',['Error!! ReconInfo ', num2str(thisReconInfo(repeatRI)), ' is/are also used in Recon ', num2str(rn)],'ForegroundColor','r','Visible','on');
                            break
                        end
                    end
                    
                    %Set Data for Table
                    set(ReconInfoTable,'Data',ReconInfos(1:ReconInfoSize,:))
                    set(ReconInfoTable,'RowName',ReconInfoRowNames);
                    set(ReconInfoTable,'ColumnName',ReconInfoColnames);
                    
                    % === more than one Recon in a event ===
                    
                    if reconEventSize > 1
                        str = ['Recon ',num2str(reconEvent(1,1))];
                        
                        for h = 2:reconEventSize
                            str = [str,' |Recon ',num2str(reconEvent(1,h))];
                        end
                        
                        uicontrol('Style','text',...
                            'Units','normalized',...
                            'FontUnits','normalized',...
                            'FontSize',textFont,...
                            'Position',[checkboxhori-.025 checkboxvert-.1 .1 .022],...
                            'HorizontalAlignment','center',...
                            'FontWeight','bold',...
                            'String','Reconstruction Select');
                        uicontrol('Style','popupmenu',...
                            'Units','normalized',...
                            'pos',[checkboxhori-.011 checkboxvert-.22 .08 .1],...
                            'String',str,...
                            'Callback',@reconNumSelect);
                    end
                    
                    reconReset = 1;
                end
                
            case 'process'
                %Get eventdata to parse the
                tableRow = eventInd(1,1)-4;
                tableCol = eventInd(1,2);
                
                %Evaluate in Event Structure
                Event = evalin('base','Event');
                processEvent = Event(tableRow,tableCol).process;
                
                %Track hold plot for TX case
                if exist('ADhndl','var') == 0 || txReset == 1 || rcvReset == 1 || reconReset == 1 || processReset == 1 || seqControlReset == 1
                    clf
                    rcvReset = 0;
                    txReset = 0;
                    reconReset = 0;
                    processReset = 0;
                    seqControlReset = 0;
                end
                
                if processEvent == 0
                    %Print Trans Title Structure
                    uicontrol('Parent',sa,'Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',messageFont,...
                        'Units','normalized',...
                        'Position',centerPos,...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'HandleVisibility', 'on',...
                        'FontWeight','bold',...
                        'String','No Process Specified...');
                    processReset = 1;
                end
                
                if processEvent ~= 0
                    errorMsg = [];
                    if evalin('base', 'exist(''Process'', ''var'')')
                        if ~ismember(processEvent,1:1:length(Process))
                            errorMsg = ['Error!! Process(', num2str(processEvent),') is used in the Event Sequence but not predefined!'];
                        end
                    else
                        errorMsg = ['Error!! Process(', num2str(processEvent),') is used in the Event Sequence but Process does not exist!'];
                    end
                    
                    if ~isempty(errorMsg)
                        clf
                        hor = 0;
                        dispmessage = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
                            'HorizontalAlignment','center',...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'String',errorMsg);
                        rcvReset = 1;
                        return
                    end
                    
                    colWidth = {160,100};
                    hori1xx = .015;
                    vert1xx = .785;
                    vert1xxOffset = .025;
                    horiSize = .215;
                    vertSize = .15;
                    
                    %Evalin Process Structure, set row names, create data for table
                    Process = evalin('base','Process');
                    RowNames = fieldnames(Process);
                    StructName = Process;
                    StructText = 'Process';
                    EventNum = processEvent;
                    EventNumText = sprintf('(%d)',EventNum); %Set in loop below
                    
                    %Call function to create data in table format
                    TableData = dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText);
                    
                    %Table verticle dim
                    ii = 1;
                    
                    %Create Event Table
                    ProcessTable = uitable('Parent',sa,...
                        'ColumnName',{'Process Object';'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx (vert1xx-(ii-1)*vertSize)-(ii-1)*vert1xxOffset horiSize vertSize]);
                    %Set Trans data in table
                    set(ProcessTable,'Data',TableData);
                    clear TableData RowNames StructName StructText
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    %NOTE WEIRD DATA FORMAT IN PROCESS.PARAMETERS SUBJECT TO CHANGE
                    %Create Table Data for Process Parameters
                    ProcessParameters = cell(size(Process(EventNum).Parameters,2)/2,2); %Preallocate array
                    for jj = 1:2:size(Process(EventNum).Parameters,2)
                        ProcessParameters{(jj+1)/2,1} = Process(EventNum).Parameters{jj};
                        ProcessParameters{(jj+1)/2,2} = num2str(Process(EventNum).Parameters{jj+1}); %Create Structure dynamically
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %Call function to create data in table format
                    TableData = ProcessParameters;
                    
                    %Create Event Table
                    ProcessParametersTable = uitable('Parent',sa,...
                        'ColumnName',{['Process' EventNumText '.Parameters'];'Entry/Value'},...
                        'FontSize',tableFontSize,...
                        'ColumnWidth', colWidth,...
                        'RowName',[],...
                        'Units','normalized',...
                        'Position',[hori1xx (vert1xx-3*vertSize)-vert1xxOffset horiSize 3*vertSize]);
                    
                    %Set Trans data in table
                    set(ProcessParametersTable,'Data',TableData);
                    clear RowNames StructName StructText TableData EventNum EventNumText ii
                    processReset = 1;
                end
                
            case 'seqControl'
                %Get eventdata to parse the
                tableRow = eventInd(1,1)-5;
                tableCol = eventInd(1,2);
                
                %Evaluate in Event Structure
                Event = evalin('base','Event');
                seqControlEvent = Event(tableRow,tableCol).seqControl;
                
                %Track hold plot for TX case
                if exist('ADhndl','var') == 0 || txReset == 1 || rcvReset == 1 || reconReset == 1 ||processReset == 1 || seqControlReset == 1
                    clf %Clears the figure window of all visible handles
                    rcvReset = 0;
                    txReset = 0;
                    reconReset = 0;
                    processReset = 0;
                    seqControlReset = 0;
                end
                
                if seqControlEvent == 0
                    %Print Trans Title Structure
                    uicontrol('Parent',sa,'Style','text',...
                        'FontUnits','normalized',...
                        'FontSize',messageFont,...
                        'Units','normalized',...
                        'Position',centerPos,...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[.7 .7 .7],...
                        'HandleVisibility', 'on',...
                        'FontWeight','bold',...
                        'String','No Sequency Control Specified...');
                    seqControlReset = 1;
                end
                
                hori = 0;
                
                if seqControlEvent ~= 0
                    
                    errorMsg = [];
                    if evalin('base', 'exist(''SeqControl'', ''var'')')
                        if ~ismember(seqControlEvent,1:1:length(SeqControl))
                            errorMsg = ['Error!! SeqControl(', num2str(seqControlEvent),') is used in the Event Sequence but not predefined!'];
                        end
                    else
                        errorMsg = ['Error!! SeqControl(', num2str(seqControlEvent),') is used in the Event Sequence but SeqControl does not exist!'];
                    end
                    
                    if ~isempty(errorMsg)
                        clf
                        hor = 0;
                        dispmessage = uicontrol('Parent',sa,'Style','text',...
                            'FontUnits','normalized',...
                            'FontSize',messageFont,...
                            'Units','normalized',...
                            'Position',centerPos+[-0.3 -0.05-hor 0.6 0],...
                            'HorizontalAlignment','center',...
                            'BackgroundColor',[.7 .7 .7],...
                            'HandleVisibility', 'on',...
                            'FontWeight','bold',...
                            'ForegroundColor','r',...
                            'String',errorMsg);
                        rcvReset = 1;
                        return
                    end
                    
                    %Display SeqControl structure in uitable
                    colWidth = {160,100};
                    hori1xx = .015;
                    vert1xx = .785;
                    vert1xxOffset = .025;
                    horiSize = .215;
                    vertSize = .15;
                    
                    %Evalin SeqControl Structure, set row names, create data for table
                    SeqControl = evalin('base','SeqControl');
                    RowNames = fieldnames(SeqControl);
                    StructName = SeqControl;
                    StructText = 'SeqControl';
                    EventNum = seqControlEvent;
                    
                    
                    %Logic to address multiple entries in Event
                    for ii = 1:size(EventNum,2)
                        EventNumSelect = EventNum(ii);
                        EventNumText = sprintf('(%d)',EventNumSelect);
                        %Call function to create data in table format
                        TableData = dataTableCreate(RowNames,StructName,StructText,EventNumSelect,EventNumText);
                        
                        %Create Event Table
                        SeqControlTable = uitable('Parent',sa,...
                            'ColumnName',{'Sequence Cntl Object';'Entry/Value'},...
                            'FontSize',tableFontSize,...
                            'ColumnWidth', colWidth,...
                            'RowName',[],...
                            'Units','normalized',...
                            'Position',[hori1xx (vert1xx-(ii-1)*vertSize)-(ii-1)*vert1xxOffset horiSize vertSize]);
                        
                        %Set Trans data in table
                        set(SeqControlTable,'Data',TableData);
                        clear TableData EventNumSelect EventNumText
                        
                    end
                    clear RowNames StructName StructText EventNum
                    seqControlReset = 1;
                    
                    %Get size of event can be more than one seqControl per event
                    seqControlEventSize = size(seqControlEvent,2);
                    seqControlFields = fieldnames(SeqControl);%store all field names in cell array
                    seqControlFieldsize = size(seqControlFields,1);%find size of the cell array created
                    
                    hori = hori+.025;
                    
                    for ns = 1:seqControlEventSize
                        
                        %Prepare string
                        for k = 1:seqControlFieldsize
                            %Format string
                            
                            if strcmp(SeqControl(1,seqControlEvent(1,ns)).(seqControlFields{k,1}),'setRcvProfile') && evalin('base','exist(''RcvProfile'',''var'')')
                                RcvProfileAll = evalin('base','RcvProfile');
                                argument = SeqControl(1,seqControlEvent(1,ns)).(seqControlFields{k+1,1});
                                RcvProfile = RcvProfileAll(:,argument);
                                RcvProfileFields = fieldnames(RcvProfile);
                                RcvProfilesize = length(RcvProfileFields);
                                %Prepare string
                                
                                for j = 1:RcvProfilesize(1)
                                    if isempty(RcvProfile.(RcvProfileFields{j,1}))
                                        RcvDisplayArray{j,1} = sprintf('RcvProfile(%g).%s = [ ] \n',argument,RcvProfileFields{j,1});
                                        %Check for array results
                                    else
                                        RcvCase = RcvProfileFields{j,1};
                                        switch RcvCase
                                            case 'VDASRcvProfile'
                                                RcvDisplayArray{j,1} = sprintf('RcvProfile(%g).%s =  < %s x %s >\n',argument,RcvProfileFields{j,1},num2str(size(RcvProfile.(RcvProfileFields{j,1}),1)),num2str(size(RcvProfile.(RcvProfileFields{j,1}),2)));
                                            case 'VDASRcvProfileSize'
                                                RcvDisplayArray{j,1} = sprintf('RcvProfile(%g).%s = [%s,%s] \n',argument,RcvProfileFields{j,1},num2str(RcvProfile.(RcvProfileFields{j,1})(1)),num2str(RcvProfile.(RcvProfileFields{j,1})(2)));
                                            otherwise
                                                RcvDisplayArray{j,1} = sprintf('RcvProfile(%g).%s =  %s \n',argument,RcvProfileFields{j,1},num2str(RcvProfile.(RcvProfileFields{j,1})));
                                        end
                                    end
                                end
                                
                                %Format RcvProfile string
                                for j = 1:RcvProfilesize(1)
                                    RcvProfileParam{j,1} = sprintf('%s \n',RcvDisplayArray{j,1});
                                end
                                
                                horiRcv = 0.04;
                                %Print RcvProfile structure
                                uicontrol('Style','text',...
                                    'FontUnits','normalized',...
                                    'FontSize',textFont,...
                                    'Units','normalized',...
                                    'Position',[hori1x+0.4 vert1x-horiRcv .25 .025],...
                                    'HorizontalAlignment','center',...
                                    'FontWeight','bold',...
                                    'String','RcvProfile Parameters');
                                horiRcv = horiRcv+.025;
                                
                                for j = 1:RcvProfilesize(1)
                                    uicontrol('Style','text',...
                                        'FontUnits','normalized',...
                                        'FontSize',textFont,...
                                        'Units','normalized',...
                                        'Position',[hori1x+0.4 vert1x-horiRcv .25 .025],...
                                        'HorizontalAlignment','left',...
                                        'String',RcvProfileParam{j,1});
                                    horiRcv = horiRcv+.025;
                                end
                            end
                            
                            % show Round Trip Time and timeToNextAcq comparision
                            if strcmp(SeqControl(1,seqControlEvent(1,ns)).command,'timeToNextAcq')
                                if any(strcmp({SeqControl(seqControlEvent).command},'triggerIn')) || any(strcmp({SeqControl(seqControlEvent).command},'sync'))
                                    % TTNA is useless with triggerIn
                                    EventInfoTxt.String = 'A timeToNextAcq is useless when triggerIn or sync is used in the same event, please remove either one.';
                                    EventInfoTxt.ForegroundColor = 'r';
                                elseif Event(tableCol).rcv > 0  % display rount trip info if the rcv is non zero in the same event
                                    minTTNA = round(Receive(Event(tableCol).rcv).endDepth/Trans.frequency)*2+eventOverhead; % in us
                                    if SeqControl(1,seqControlEvent(1,ns)).argument < minTTNA
                                        EventInfoTxt.String = ['Estimated round trip time + event overhead: ',num2str(minTTNA), ' us. A timeToNextAcq warning might be displayed'];
                                        EventInfoTxt.ForegroundColor = 'r';
                                    else
                                        EventInfoTxt.String = ['Estimated round trip time + event overhead: ',num2str(minTTNA), ' us.'];
                                    end
                                end
                            end
                            
                            
                            %Plot DMA Frame
                            if strcmp(SeqControl(1,seqControlEvent(1,ns)).command,'transferToHost')
                                endDMAevent = tableCol;
                                startDMAevent = 1;
                                
                                % Get current buffer and frame number
                                for n = endDMAevent:-1:1
                                    if EventStruct(n).rcv > 0
                                        dmaBufNum = Receive(EventStruct(n).rcv).bufnum;
                                        dmaFrameNum = Receive(EventStruct(n).rcv).framenum;
                                        break;
                                    end
                                end
                                
                                % Search backwards until the event of:
                                % 1) another transerToHost
                                % 2) different frame of the same buffer
                                % 3) different buffer number
                                stopByDMA = 1; % The search should stop at another DMA or the first event
                                for n = endDMAevent-1:-1:1
                                    if EventStruct(n).seqControl > 0
                                        if any(strcmp({SeqControl(EventStruct(n).seqControl).command},'transferToHost'))
                                            startDMAevent = n+1;
                                            break
                                        end
                                    end
                                    if EventStruct(n).rcv > 0
                                        if ~isequal(Receive(EventStruct(n).rcv).bufnum,dmaBufNum) || ...
                                                ~isequal(Receive(EventStruct(n).rcv).framenum,dmaFrameNum)
                                            startDMAevent = n+1;
                                            stopByDMA = 0;
                                            break;
                                        end
                                    end
                                end
                                
                                dmaRcvInd = [EventStruct(startDMAevent:endDMAevent).rcv];
                                dmaRcvInd = dmaRcvInd(dmaRcvInd ~= 0);
                                
                                % acqNum must be consecutive in one DMA
                                minAcq = min([Receive(dmaRcvInd).acqNum]);
                                maxAcq = max([Receive(dmaRcvInd).acqNum]);
                                
                                %Plot the DMA frame
                                Resource = evalin('base','Resource');
                                
                                %Position variables for receive
                                plot1rcvhori = .365;
                                %                                 plot2rcvhori = .685;
                                %                                 plot1rcvvert = .825;
                                plot2rcvvert = .3;
                                
                                numFrames = Resource.RcvBuffer(Receive(dmaRcvInd(1)).bufnum).numFrames;
                                RcvInd = arrayfun(@(x) isequal(x.bufnum,Receive(dmaRcvInd(1)).bufnum),Receive);
                                numAcqs = max(arrayfun(@(x) x.acqNum, Receive(RcvInd)));
                                
                                AcqArray = zeros(numAcqs,numFrames);
                                AcqDispArray = zeros(numAcqs,numFrames);
                                AcqArray(Receive(1,dmaRcvInd(1)).acqNum,Receive(1,dmaRcvInd(1)).framenum) = Receive(1,dmaRcvInd(1)).acqNum;
                                AcqDispArray(minAcq:maxAcq,Receive(1,dmaRcvInd(1)).framenum) = 1;
                                
                                %Image must be displayed two times in order to line up grid
                                %axes
                                rcvImg = axes('Parent',sa,'units','normalized','position',[plot1rcvhori plot2rcvvert .6 .45]);
                                imagesc(AcqDispArray);
                                set(rcvImg, 'layer', 'top');
                                if numFrames <= 30
                                    set(rcvImg, 'Xtick',0:1:numFrames);
                                end
                                if size(AcqArray(:,1),1) <40
                                    set(rcvImg, 'Ytick',0:1:size(AcqArray(:,1),1));
                                end
                                set(rcvImg, 'xcolor', 'w', 'ycolor', 'w');
                                set(rcvImg, 'YDir', 'reverse')
                                hold on
                                %Set the grid at .5 between the number ticklabel
                                h = axes('Parent',sa,'units','normalized','position',[plot1rcvhori plot2rcvvert .6 .45]);
                                imagesc(AcqDispArray)
                                if numFrames <= 30
                                    set(h, 'Xtick',.5:1:numFrames);
                                end
                                if size(AcqArray(:,1),1) < 40
                                    set(h, 'Ytick',.5:1:size(AcqArray(:,1),1));
                                end
                                set(h, 'yGrid', 'on');
                                set(h, 'xGrid', 'on');
                                set(h, 'Yticklabel', []);
                                set(h, 'Xticklabel', []);
                                set(h, 'xcolor', 'w', 'ycolor', 'w');
                                title(['The acquired data below will be transferred to the host PC: Receive Buffer ',num2str(Receive(dmaRcvInd(1)).bufnum)],'fontsize',14,'color','w')
                                %Set the position of xlabel
                                xlabel('Frame Number','fontsize',14,'units','normalized')
                                xlabh = get(gca,'XLabel');
                                set(xlabh,'Position',get(xlabh,'Position') - [0 .01 0])
                                %Set the position of y label
                                ylabel('Acquisition Number','fontsize',14,'units','normalized')
                                ylabh = get(gca,'YLabel');
                                set(ylabh,'Position',get(ylabh,'Position') - [.01 0 0])
                                hold off
                                
                                % display DMA transfer rate and possible warning message
                                dmaSizeTxt = uicontrol('Parent',sa,'Style','text',...
                                    'FontUnits','normalized',...
                                    'FontSize',messageFont,...
                                    'Units','normalized',...
                                    'Position',centerPos + [0.275 0.25 0 0],...
                                    'HorizontalAlignment','left',...
                                    'BackgroundColor',[.7 .7 .7],...
                                    'HandleVisibility', 'on',...
                                    'FontWeight','bold',...
                                    'String','Frame size information is not available without initialization.');
                                
                                if alreadyIniz
                                    seqNum = seqControlEvent(1,ns);
                                    
                                    dmaSizeTxt.String = ['Frame size: ', num2str(sizeOfEachDMA(seqNum == indDMA),'%.2f'), ' MB'];
                                    dmaSpeedTxt = copy(dmaSizeTxt);
                                    set(dmaSpeedTxt,'Parent',sa, 'Position',dmaSizeTxt.Position + [0.2 0 0 0],...
                                        'String',['DMA speed: ', num2str(dmaSpeed,'%.2f'), ' GB/s']);
                                    
                                    %                                 dmaOverheadTxt = copy(dmaSizeTxt);
                                    %                                 set(dmaOverheadTxt,'Parent',sa, 'Position',dmaSizeTxt.Position + [0 -0.05 0 0],...
                                    %                                     'String',['Estimated overhead for DMA initialization: ', num2str(dmaOverhead), ' ms']);
                                    
                                    dmaTimeTxt = copy(dmaSizeTxt);
                                    set(dmaTimeTxt,'Parent',sa, 'Position',dmaSizeTxt.Position + [0 -0.05 0.2 0],...
                                        'String',['Estimated duration for this DMA (including initialization time): ', num2str(timeForEachDMA(seqNum == indDMA),'%.2f'), ' ms']);
                                    
                                    dmaIntervalTxt = copy(dmaSizeTxt);
                                    set(dmaIntervalTxt,'Parent',sa, 'Position',dmaTimeTxt.Position + [0 -0.05 0 0],...
                                        'String',['Interval to next DMA: ', num2str(intervalToNextDMA(seqNum == indDMA),'%.2f'), ' ms']);
                                    
                                    currentDMAtime = timeForEachDMA(seqNum == indDMA);
                                    currentDMAinterval = intervalToNextDMA(seqNum == indDMA);
                                    
                                    
                                    if isequal(currentDMAinterval, 0)
                                        dmaIntervalTxt.Visible = 'off';
                                    elseif  currentDMAtime > currentDMAinterval
                                        dispValue = ceil(10*(currentDMAtime - currentDMAinterval));
                                        dmaIntervalTxt.String = ['Interval to next DMA: ', num2str(currentDMAinterval,'%.2f'),...
                                            ' ms. --> timeToNextAcq duration too short warning, around ', num2str(dispValue+2,'%.0f'),'00 +- 250 us'];
                                        EventInfoTxt.String = 'Warning!! A timeToNextAcq warning might be displayed because next DMA will not be initialized before the completion of current DMA';
                                        EventInfoTxt.ForegroundColor = 'r';
                                    end
                                end
                                
                            end
                            
                        end
                        
                    end
                    seqControlReset = 1;
                end
        end
    end

%Data Table Creation Function
    function TableData =  dataTableCreate(RowNames,StructName,StructText,EventNum,EventNumText)
        %%%% Note Function is used to create/organize data to be used in tables: Trans,TX,TW,Receive,Process,and SeqControl
        %Data for Table
        TableData = cell(size(RowNames,1),2); %Preallocate array
        for fn = 1:size(RowNames,1) %Loop for # of entries in structure
            TableData{fn,1} = sprintf('%s%s.%s',StructText,EventNumText,RowNames{fn,1});
            %If data type is double and dimensisons are greater than [1,1] show dimensions and type
            if ~ischar(StructName(EventNum).(RowNames{fn})) && ~isequal(size(StructName(EventNum).(RowNames{fn})),[1,1])
                if size(StructName(EventNum).(RowNames{fn}),1) == 1 && size(StructName(EventNum).(RowNames{fn}),2) < 5 && isnumeric(StructName(EventNum).(RowNames{fn}))
                    if size(StructName(EventNum).(RowNames{fn}),2) == 2
                        TableData{fn,2} = sprintf('[%2.4g, %2.4g]',StructName(EventNum).(RowNames{fn})(1),StructName(EventNum).(RowNames{fn})(2));
                    elseif size(StructName(EventNum).(RowNames{fn}),2) == 3
                        TableData{fn,2} = sprintf('[%2.4g, %2.4g, %2.4g]',StructName(EventNum).(RowNames{fn})(1),StructName(EventNum).(RowNames{fn})(2),StructName(EventNum).(RowNames{fn})(3));
                    elseif size(StructName(EventNum).(RowNames{fn}),2) == 4
                        TableData{fn,2} = sprintf('[%2.4g, %2.4g, %2.4g, %2.4g]',StructName(EventNum).(RowNames{fn})(1),StructName(EventNum).(RowNames{fn})(2),StructName(EventNum).(RowNames{fn})(3),StructName(EventNum).(RowNames{fn})(4));
                    end
                else
                    if isreal(StructName(EventNum).(RowNames{fn})) || iscell(StructName(EventNum).(RowNames{fn})) %Class type
                        dataType = class(StructName(EventNum).(RowNames{fn}));
                    else
                        dataType = 'complex double';
                    end
                    %Set the output in table to the size of the array
                    if isempty(StructName(EventNum).(RowNames{fn}))
                        TableData{fn,2} = '[ ]';
                    else
                        TableData{fn,2} = sprintf('%dx%d %s',size(StructName(EventNum).(RowNames{fn}),1),size(StructName(EventNum).(RowNames{fn}),2),dataType);
                    end
                end
            elseif strcmpi(class(StructName(EventNum).(RowNames{fn})),'double')
                TableData{fn,2} = num2str(StructName(EventNum).(RowNames{fn}));
            elseif isstruct(StructName(EventNum).(RowNames{fn}))
                %Set the output in table to the size of the structure
                TableData{fn,2} = sprintf('%dx%d %s',size(StructName(EventNum).(RowNames{fn}),1),size(StructName(EventNum).(RowNames{fn}),2),'struct');
            else
                TableData{fn,2} = StructName(EventNum).(RowNames{fn});
            end
        end
    end

%TW Plot Callback
    function selcTW(~,~)
        txEvent = evalin('base','txEvent'); %Current TX highlighted
        TW = evalin('base','TW');
        TX = evalin('base','TX');
        TWnum = TX(1,txEvent).waveform; %Use current TX to choose TW numbe
        if ~isfield(TW,'Wvfm1Wy')
            [~,~,~,~,TW] = computeTWWaveform(TW);
        end
        delete(findobj('Name','Transmit Waveform Plot'));
        figure('Name','Transmit Waveform Plot');
        %Plot simulated waveform response
        TWPlotSim = subplot(2,1,1);
        plot(TWPlotSim,TW(TWnum).Wvfm1Wy);
        %Plot tristate response
        TWPlotTriState = subplot(2,1,2);
        plot(TWPlotTriState,TW(TWnum).TriLvlWvfm_Sim);
    end

%TXPD Plot Callback
    function selcTXPD(~,~)
        txEvent = evalin('base','txEvent');
        exist('Region','var');
        showTXPD(txEvent);
    end

%RcvBuffPlotReset
    function rcvbuffplotreset(hObject,~)
        if (get(hObject,'Value') == get(hObject,'Max'))
            imagehold  = 1;
            imageholdReset = 1;
        else
            imagehold = 0;
            imageholdReset = 0;
        end
    end

%reconNumSelect
    function reconNumSelect(hObject,~)
        reconVal =  get(hObject,'Value');
        
        %Evalin appropriate variables
        Recon = evalin('base','Recon');
        ReconInfo = evalin('base','ReconInfo');
        
        %Set the tableReconEvent
        tableReconEvent = reconEvent(1,reconVal);
        
        % === 1st uitable shows Source and Output ===
        if ~isfield(Recon(1,tableReconEvent),'rcvBufFrame')
            Recon(1,tableReconEvent).rcvBufFrame = [];
        end
        
        InOutData = {[blanks(8),'RcvBuffer'],[blanks(8),'InterBuffer'],[blanks(8),'ImageBuffer'];...
            '[ ]','[ ]','[ ]';...
            [blanks(15),num2str(Recon(1,tableReconEvent).rcvBufFrame)],'[ ]','[ ]'};
        if isfield(Recon(1,tableReconEvent),'IntBufDest') && ~isempty(Recon(1,tableReconEvent).IntBufDest)
            InOutData(2:3,2) = {[blanks(15),num2str(Recon(1,tableReconEvent).IntBufDest(1))];...
                [blanks(15),num2str(Recon(1,tableReconEvent).IntBufDest(2))]};
        end
        if isfield(Recon(1,tableReconEvent),'ImgBufDest') && ~isempty(Recon(1,tableReconEvent).ImgBufDest)
            InOutData(2:3,3) = {[blanks(15),num2str(Recon(1,tableReconEvent).ImgBufDest(1))];...
                [blanks(15),num2str(Recon(1,tableReconEvent).ImgBufDest(2))]};
        end
        
        % === 2nd uitable shows ReconInfo and associated Receive Structures ===
        
        %Get the ReconInfo Field names, only four fields are
        %required for eventanalysis
        ReconInfoNew = struct('mode',{ReconInfo.mode},...
            'txnum',{ReconInfo.txnum},...
            'rcvnum',{ReconInfo.rcvnum},...
            'pagenum',[],...
            'regionnum',{ReconInfo.regionnum});
        
        if isfield(ReconInfo,'pagenum')
            [ReconInfoNew(:).pagenum] = deal(ReconInfo.pagenum);
        end
        
        ReconInfoColnames = {'mode','txnum','rcvnum','bufnum','framenum','acqNum','pagenum','regionnum'};
        
        ReconInfo = ReconInfoNew;
        Receive = evalin('base','Receive');
        
        %Get the number of how many field names
        sizeReconInfoCol = size(ReconInfoColnames,2);
        ReconInfoSize = length(Recon(1,tableReconEvent).RINums);
        
        %Preallocate Event Cell array
        ReconInfos = cell(ReconInfoSize,sizeReconInfoCol); %Get the total number of events in the sequence
        
        for n = 1:sizeReconInfoCol
            %Set field name for reconinfo
            reconInfoFieldName = ReconInfoColnames{1,n};
            %Use switch statement to input data in the reconinfo structure
            switch reconInfoFieldName
                case 'mode'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).mode);
                    end
                    
                case 'txnum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).txnum);
                    end
                    
                case 'rcvnum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).rcvnum);
                    end
                    
                case 'bufnum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        rcvNum = ReconInfo(RIindex).rcvnum;
                        ReconInfos{k,n}  = num2str(Receive(rcvNum).bufnum);
                    end
                    InOutData{2,1} = [blanks(15),num2str(Receive(rcvNum).bufnum)];
                    % If the RcvBuffer among all ReconInfos is not the consistent, show Warning
                    if ReconInfoSize > 1
                        RcvNum = ReconInfo(Recon(1,tableReconEvent).RINums(1)).rcvnum:ReconInfo(Recon(1,tableReconEvent).RINums(ReconInfoSize)).rcvnum;
                        [BufNum,~,ind] = unique([Receive(RcvNum).bufnum]);
                        in = find(ind~=mode(ind));
                        if numel(BufNum) > 1
                            for rowNum = in(1):in(end)
                                ReconInfos(rowNum,n) = {horzcat('<html><table border=0 width=65 bgcolor="red"><TR><TD>',ReconInfos{rowNum,n},'</TD></TR> </table></html>')};
                            end
                            %Print warning message
                            set(ReconInfoWarning,'String','Error!! The bufnum is not consistent!','ForegroundColor','r','Visible','on');
                            InOutData(2,1) = {horzcat('<html><table border=0 width=115 bgcolor="red"><TR><TD>Inconsistent!</TD></TR> </table></html>')};
                        end
                    end
                    
                case 'framenum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        rcvNum = ReconInfo(RIindex).rcvnum;
                        ReconInfos{k,n}  = num2str(Receive(rcvNum).framenum);
                    end
                    
                    RcvNum = ReconInfo(Recon(1,tableReconEvent).RINums(1)).rcvnum:ReconInfo(Recon(1,tableReconEvent).RINums(ReconInfoSize)).rcvnum;
                    
                    if numel(unique([Receive(RcvNum).endDepth])) > 1
                        depthWarning = copy(ReconInfoWarning);
                        depthWarning.Position(2) = vert1x-0.66;
                        set(depthWarning,'Parent',sa,'String',['Error!! Receive.endDepth across Receive ', num2str(RcvNum(1)),' to ', ...
                            num2str(RcvNum(end)),' must be identical for Recon ',num2str(tableReconEvent)],'ForegroundColor','r','Visible','on');
                    end
                    
                    % All frame number in the ReconInfo should be consistent before beling replaced by the Recon.rcvBufFrame
                    [frameNum,~,ind] = unique([Receive(RcvNum).framenum]);
                    in = find(ind~=mode(ind));
                    if numel(frameNum) > 1
                        for rowNum = in(1):in(end)
                            ReconInfos(rowNum,n) = {horzcat('<html><table border=0 width=65 bgcolor="red"><TR><TD>',ReconInfos{rowNum,n},'</TD></TR> </table></html>')};
                        end
                        set(ReconInfoWarning,'String','Error!! The framenum is not consistent!','Visible','on');
                        %                         InOutData(3,1) = {horzcat('<html><table border=0 width=115 bgcolor="red"><TR><TD>Inconsistent!</TD></TR> </table></html>')};
                    else
                        ReconFrame = Recon(1,tableReconEvent).rcvBufFrame;
                        if ~isempty(ReconFrame)
                            % if Recon.rcvBufFrame is not empty, applys it to the frame no. in Reconfo and display warning
                            ReconInfos(1:ReconInfoSize,n) = {num2str(ReconFrame)};
                            set(ReconInfoWarning,'String','Note! Recon.rcvBufFrame over-rides the frame number specified by the ReconInfo','ForegroundColor','k','Visible','on');
                        end
                    end
                    
                case 'acqNum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        rcvNum = ReconInfo(RIindex).rcvnum;
                        ReconInfos{k,n}  = num2str(Receive(rcvNum).acqNum);
                    end
                    
                case 'pagenum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).pagenum);
                    end
                    
                case 'regionnum'
                    for k = 1:ReconInfoSize
                        RIindex = Recon(1,tableReconEvent).RINums(k);
                        ReconInfos{k,n} = num2str(ReconInfo(1,RIindex).regionnum);
                    end
            end
        end
        
        for k = 1:ReconInfoSize
            RIindex = Recon(1,tableReconEvent).RINums(k);
            ReconInfoRowNames{:,k} = ['RI ', num2str(RIindex, '%03.0f')];
        end
        
        % Check if ReconInfo is used in another Recon
        thisReconInfo = [Recon(1,tableReconEvent).RINums];
        ReconInd = 1:length(Recon);
        ReconInd(tableReconEvent) = [];
        for rn = ReconInd
            repeatRI = ismember(thisReconInfo,[Recon(rn).RINums]);
            if any(repeatRI)
                set(ReconInfoWarning,'String',['Error!! ReconInfo ', num2str(thisReconInfo(repeatRI)), ' is/are also used in Recon ', num2str(rn)],'ForegroundColor','r','Visible','on');
                break
            end
        end
        
        %Set InOutTable
        set(InOutTable,'Data',InOutData);
        
        %Set ReconInfo table
        set(ReconInfoTable,'Data',ReconInfos)
        set(ReconInfoTable,'RowName',ReconInfoRowNames);
    end

% Voltage droop estimation for profile 5
    function plotVoltDroop(varargin)
        
        figFont = 0.028;
        % recalculate pushcapVdroop for steady status
        TW = evalin('base','TW');
        TPC = evalin('base','TPC');
        TXEvent = evalin('base','TX_Limits.TXEvent');
        
        % Set Voltage to HVLim
        % voltSldr.Value = TXEvent(1).HVLim + 0.8;
        
        auxHVPsImax = 0.5; % maximum output current in Amps for aux HV supply from TPC
        Cpush = 15e3; % push capacitor is 15 mF, or 15,000 uF (uF units are used here since time is in usec)
        
        prevVdroop = 0;
        maxVdroop = 0;
        oldmaxVdroop = 0;
        rising = 1; % controls when to exit the while loop
        numTXEvent = length(TXEvent);
        
        while rising
            for TEnum = 1:numTXEvent
                Bdur(TEnum) = TW(TXEvent(TEnum).tw).Bdur;
                idleT(TEnum) = TXEvent(TEnum).cumPRI - Bdur(TEnum);
                HVLim = TXEvent(TEnum).HVLim;
                HV = HVLim - prevVdroop;
                TXEvent(TEnum).pushcapVdroop = prevVdroop + Bdur(TEnum) * max(0, TXEvent(TEnum).TotalHVPowerIn*HV/HVLim^2 - auxHVPsImax) / Cpush;
                maxVdroop = max(maxVdroop, TXEvent(TEnum).pushcapVdroop);
                prevVdroop = max(0, TXEvent(TEnum).pushcapVdroop - auxHVPsImax*idleT(TEnum)/Cpush);
            end
            if floor(maxVdroop*1e4)<=floor(oldmaxVdroop*1e4)
                rising = 0;
            end
            oldmaxVdroop = maxVdroop;
        end
        droopPercent = min(100, 100 * TXEvent(1).pushcapVdroop/TXEvent(1).HVLim);
        
        % Plot
        
        ScrnSize = get(groot,'ScreenSize');
        fWidth = 700; fHighth = 400;
        
        % Close any previously opened GUI windows.
        hf = findobj('tag','figVdroop');
        if ishandle(hf)
            pos = get(hf,'Position');
            delete(hf);
        else
            pos = [(ScrnSize(1,3)-fWidth)/2,(ScrnSize(1,4)-fHighth)/2,fWidth,fHighth];
        end
        
        figVdroop = figure('Visible','on',...  %'Units','normalized',...
            'Position',pos,... %'Position',[0.7,0.25,0.25,0.50],...
            'Name','Voltage droop visualization',...
            'NumberTitle','off',...
            'MenuBar','none', ...
            'Toolbar','figure',...
            'Resize','on', ...
            'tag','figVdroop');
        
        axesVdroop = axes(...
            'Parent',figVdroop,...
            'Units','normalized',...
            'Position',[0.06,0.1,0.9,0.8],...
            'NextPlot','add',...
            'FontUnits','normalized',...
            'FontSize',figFont,...
            'Tag','droopAxes');
        
        for TEnum = 1:numTXEvent
            if isequal(TEnum,1)
                plotV(TEnum).vStart = TXEvent(TEnum).HVLim;
                plotV(TEnum).xStart = 0;
            else
                plotV(TEnum).vStart = plotV(TEnum-1).vEnd+plotV(TEnum-1).vCharge;
                plotV(TEnum).xStart = plotV(TEnum-1).xStart + TXEvent(TEnum).cumPRI/1e3;
            end
            plotV(TEnum).xEnd = plotV(TEnum).xStart + Bdur(TEnum)/1e3;
            plotV(TEnum).vEnd = plotV(TEnum).vStart - TXEvent(TEnum).pushcapVdroop;
            if plotV(TEnum).vEnd <= 0
                plotV(TEnum).vEnd = 0;
            end
            
            plotV(TEnum).vCharge = auxHVPsImax*idleT(TEnum)/Cpush;
            plotV(TEnum).tCharge = idleT(TEnum)/1e3;
            plotV(TEnum).noCharge = 0;
            if plotV(TEnum).vEnd + plotV(TEnum).vCharge>= plotV(TEnum).vStart
                plotV(TEnum).vCharge = plotV(TEnum).vStart-plotV(TEnum).vEnd;
                plotV(TEnum).tCharge = plotV(TEnum).vCharge*Cpush/auxHVPsImax/1e3;
                plotV(TEnum).noCharge = TXEvent(TEnum).cumPRI/1e3-Bdur(TEnum)/1e3-plotV(TEnum).tCharge;
            end
            
            line(axesVdroop,[plotV(TEnum).xStart plotV(TEnum).xEnd],[plotV(TEnum).vStart plotV(TEnum).vEnd],'LineWidth',2);
            if isequal(plotV(TEnum).noCharge,0)
                line(axesVdroop,[plotV(TEnum).xEnd TEnum*TXEvent(TEnum).cumPRI/1e3],[plotV(TEnum).vEnd plotV(TEnum).vEnd+plotV(TEnum).vCharge],...
                    'LineWidth',2,'Color','red','LineStyle','--');
                legend('Droop','Charge')
            else
                line(axesVdroop,[plotV(TEnum).xEnd plotV(TEnum).xEnd+plotV(TEnum).tCharge],[plotV(TEnum).vEnd plotV(TEnum).vEnd+plotV(TEnum).vCharge],...
                    'LineWidth',2,'Color','red','LineStyle','--');
                line(axesVdroop,[plotV(TEnum).xEnd+plotV(TEnum).tCharge plotV(TEnum).xEnd+plotV(TEnum).tCharge+plotV(TEnum).noCharge],...
                    [plotV(TEnum).vStart plotV(TEnum).vStart],'LineWidth',2,'Color','g','LineStyle','--');
                legend('Droop','Charge','Idle')
            end
            grid on
            axis([0 numTXEvent*TXEvent(TEnum).cumPRI/1e3 0 TPC(5).maxHighVoltage+20]);
            xlabel('Time (ms)','FontUnits','Normalized','FontSize',figFont);
            ylabel('Voltage (V)','FontUnits','Normalized','FontSize',figFont);
            title({'Transmit Voltage vs. Time   (internal extended burst supply)', ['Droop during burst = ' num2str(droopPercent,'%3.1f') ' %' ]}, 'fontsize', 14)
        end
        
    end

% NOTE: the getJtable is copied from the findjobj and the copyright
% is Copyright (c) 2018, Yair Altman
% All rights reserved.

    function jhandle = getJtable(hControl) % input is the Table handle and the outout is its java properties
        warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
        jf = get(sa,'JavaFrame');
        jContainer = jf.getFigurePanelContainer.getComponent(0);
        
        jhandle = [];
        counter = 20;
        
        specialTooltipStr = '!@#$%^&*';
        try  % Fix for R2018b suggested by Eddie (FEX comment 2018-09-19)
            tooltipPropName = 'TooltipString';
            oldTooltip = get(hControl,tooltipPropName);
            set(hControl,tooltipPropName,specialTooltipStr);
        catch
            tooltipPropName = 'Tooltip';
            oldTooltip = get(hControl,tooltipPropName);
            set(hControl,tooltipPropName,specialTooltipStr);
        end
        
        while isempty(jhandle) && counter>0
            counter = counter - 1;
            pause(0.005);
            jControl = findTooltipIn(jContainer,specialTooltipStr);
        end
        set(hControl,tooltipPropName,oldTooltip);
        try jControl.setToolTipText(oldTooltip); catch, end
        try jhandle = jControl.getParent.getView.getParent.getParent; catch, end
        jhandle = handle(jhandle,'callbackproperties');
    end

    function jControl = findTooltipIn(jContainer, specialTooltipStr)
        try
            jControl = [];  % Fix suggested by H. Koch 11/4/2017
            tooltipStr = jContainer.getToolTipText;
            %if strcmp(char(tooltipStr),specialTooltipStr)
            if ~isempty(tooltipStr) && tooltipStr.startsWith(specialTooltipStr)  % a bit faster
                jControl = jContainer;
            else
                for idx = 1 : jContainer.getComponentCount
                    jControl = findTooltipIn(jContainer.getComponent(idx-1), specialTooltipStr);
                    if ~isempty(jControl), return; end
                end
            end
        catch
            % ignore
        end
    end

    function closefunc(~,~)
        if evalin('base','exist(''tempFileForEA.mat'',''file'');')
            evalin('base','delete tempFileForEA.mat');
        end
        delete(sa);
    end

end
