function deployVSX(varargin)
% deployVSX(varargin)
%
% 'scriptname',     path and name of script
% -c,               close all windows
% -s#,              skip multicheck and use # of secondaries
% -T(x or g),       which terminal program to use 'xterm' or 'gnome'
% default is gnome
% -L(w or i)        enable logger (warning or info) default is not enabled
%
% Example:
% deployVSX './Example_Scripts/Volume_Imaging_Package/4system/SetUpVermon1024_Flash_IQ.m' -c -s3 -Tg -Lw
%
warning('deployVSX:deprecated', 'deployVSX has been deprecated and should be replaced with vsv.multi.deploy');

% --- Set Default Flags if there are none command line --%
%defaultFlags = '-c -Tg ';  % User should adjust this line to suit their needs for a typical run
defaultFlags = '-c -Tg -Li';  % User should adjust this line to suit their needs for a typical run

% Only start if using Linux
if ~(isunix && ~ismac)
    error('deployVSX: Linux Operating System is required for use of Volume Imaging Package.');
end

% Parse the input parameters into 'scriptname' and 'flags'.
if nargin == 0
    scriptname = '';
    flags = defaultFlags;
elseif nargin == 1
    if(startsWith(varargin,'-'))
        scriptname = '';
        flags = varargin{1};
    else
        target = varargin{1};
        flags = defaultFlags;
    end
else
    if(startsWith(varargin,'-'))
        scriptname = '';
        flags = varargin;
    else
        scriptname = varargin{1};
        flags = varargin(2:length(varargin));
    end
end

% Convert the flags array into a single string (for easier parsing).
if(iscell(flags))
    flags = strjoin(flags, ' ');
end


%-- Close windows?
if(contains(flags, '-c'))
    %--- close al windows ---%
    set(findobj('Name','VSX Control'),'CloseRequestFcn','')
    close all force
    VsClose
end

%-- MultiCheck?
if(contains(flags, '-s'))
    %--- skip multicheck ---%
    s=regexpi(flags,'-s(\d+)','tokens');numSecondaries=str2double(s{1});
    if (numSecondaries~=3)&&(numSecondaries~=4)
        error('Incorrect number of secondary systems entered. must be 3 or 4.')
    else
        fprintf('Skipping multiCheck & using %i secondary systems.\n', numSecondaries);
    end
else
    %--- perform multicheck ---%
    numSecondaries = multiCheck(0);  %check secondaries and reset if needed result is # of secondaries (3 or 4) or negative for error code
    if (numSecondaries==0)
        errror('Not all secondary systems are online. Please run "multiCheck(1)" for complete diagnostics');
    else
        setenv('NUMSLAVES',num2str(numSecondaries));
    end
end

%--- Get scriptname from user ---%
if isempty(scriptname)
    scriptname = input('Name of MultiSystem .m file to process (or enter for GUI):','s');
    if isempty(scriptname)
        [filename, pathname] = uigetfile('*.m', 'Pick a MultiSystem Script:');
        if isempty(filename)
            error('no script selected.')
        else
            scriptname = [pathname filename];
        end
    end
end


%% -- Setup environment --%
setenv('LD_LIBRARY_PATH', '/usr/lib/x86_64-linux-gnu/libexpat.so.1') %this is needed for starting gnome-terminal from Matlab. If xterm is used, this is not needed
SWDIRECTORY = getenv('VERASONICS_VPF_ROOT');

%--- Close Dead Terminals from previous runs ---%
unix('tmux kill-session -t multisys'); %kill existing terminal if exists
DIRNAME=pwd;

%-- Start new terminal --%
if contains(flags, '-Tx')
    unix('xterm -maximize -e "tmux new-session -s multisys -n ''MultiSys'' " &')  %to run on Xterm instead of gnome-terminal
elseif contains(flags, '-Tg')||~contains(flags,'-T')
    unix('gnome-terminal --window --maximize -- bash -c "tmux new-session -s multisys -n ''MultiSys'';exec bash"');
end

%-- Setup Tmux Window --%
switch numSecondaries
    case 4
        unix('tmux split-window -v ; tmux split-window -h ; tmux selectp -t 0 ; tmux split-window -h'); %4 panes
    case 3
        %unix('tmux split-window -h -p 66; tmux split-window -h -p 50'); %3 vertical panes
        unix('tmux split-window -v -p 66; tmux split-window -v -p 50'); %3 horizontal panes
    case 2
        unix('tmux split-window -v -p 50'); %3 horizontal panes        
    case 1
        %placeholder. no need to split window
end


unix(['tmux source-file ' SWDIRECTORY '/System/+vsv/+multi/tmux.conf']);


% 3. connect to each secondary
IPlist={'10.10.0.2','10.10.1.2','10.10.2.2','10.10.3.2'};
IPlist = IPlist(5-numSecondaries:4);
for a = 1:numSecondaries
    cmd = sprintf('tmux send-keys -t multisys:0.%i ''ssh %s'' Enter',(a),IPlist{a});
    unix(cmd);
end

% 4. on all secondaries
unix('tmux setw synchronize-panes on'); %multicast on

% cd to workspace
unix(sprintf('tmux send-keys -t multisys:0 "cd %s" Enter',pwd));

% clear any existing matlab processes
unix('tmux send-keys -t multisys:0 ''killall -9 MATLAB'' Enter');

% Start Matlab
[status, result] = unix('tmux send-keys -t multisys:0 ''matlab -r "run activate.m"'' Enter'); %normal
pause(7); %ugly fix. the line above needs to really block until matlab has started.

unix('tmux setw synchronize-panes off'); %multicast off


%--- Execute script remotely/Locally ---%
fprintf('Execute scripts remotely\n')
for a = 1:numSecondaries
    cmd = sprintf('tmux send-keys -t multisys:0.%i "run(''%s'')" Enter',(a),scriptname);
    unix(cmd);
end
fprintf('Execute script locally\n')
evalin('base',['run(''' scriptname ''')']);

%--- Start logger remotely & locally---%
loggerFlag = regexp(flags,'-L(\w)','tokens');
if ~isempty(loggerFlag)
    fprintf('Starting Logger in %s mode\n',char(loggerFlag{1}));
    lineCmd = sprintf('vsv.common.util.Logger.start(''%s'',''all'',''_LogFile_N%i'')',char(loggerFlag{1}), 0);
    evalin('base',lineCmd);
    for a = 1:numSecondaries
        lineCmd = sprintf('vsv.common.util.Logger.start(''%s'',''all'',''_LogFile_N%i'')',char(loggerFlag{1}), a);
        cmd = sprintf('tmux send-keys -t multisys:0.%i "%s" Enter',(a),lineCmd);
        unix(cmd);
    end
end

%--- Start VSX remotely/locally ---%
fprintf('Start VSX remotely\n')
lineCmd = 'VSX';
for a = 1:numSecondaries
    cmd = sprintf('tmux send-keys -t multisys:0.%i "%s" Enter',(a),lineCmd);
    unix(cmd);
end
fprintf('Start VSX Locally\n')
evalin('base','VSX');

%--- Stop logger ---%
if ~isempty(loggerFlag)
    fprintf('Stop logger remotely/locally\n')
    lineCmd = 'vsv.common.util.Logger.stop';
    for a = 1:numSecondaries
        cmd = sprintf('tmux send-keys -t multisys:0.%i "%s" Enter',(a),lineCmd);
        unix(cmd);
    end
    evalin('base',lineCmd);
end