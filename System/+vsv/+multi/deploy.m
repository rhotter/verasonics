function deploy(varargin)
% vsv.multi.deploy(varargin)
%
% 'scriptname',     path and name of script
% -c,               close all windows before launching
% -s<N>,            use <N> number of secondary systems
% -T(x or g),       which terminal program to use 'xterm' or 'gnome' default is gnome
% -L(w or i)        enable logger (warning or info) default is to log everything
%
% Example:
% vsv.multi.deploy './Example_Scripts/Volume_Imaging_Package/4system/SetUpVermon1024_Flash_IQ.m' -c -s4 -Tg -Lw

% Copyright 2001-2023 Verasonics, Inc.  All world-wide rights and remedies
% under all intellectual property laws and industrial property laws are
% reserved.  Verasonics Registered U.S. Patent and Trademark Office.

% --- Set Default Flags if there are none command line --%
defaultFlags = '-c -s3 -Tg -Li';  % User can adjust this line to suit their needs for a typical run

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
        scriptname = varargin{1};
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

%-- Get # of secondaries
if(contains(flags, '-s')) %if there is flag use this value else, autodetect
    s=regexpi(flags,'-s(\d+)','tokens');
    numSecondaries=str2double(s{1});
    if (numSecondaries<1)||(numSecondaries>4)
        error('Incorrect number of secondary systems entered. must be 1, 2, 3 or 4.')
    else
        fprintf('Skipping multiCheck & using %i secondary systems.\n', numSecondaries);
    end
else
    numSecondaries = vsv.multi.getNumSecondaries;  %check secondaries and reset if needed result is # of secondaries (3 or 4) or negative for error code
    if (numSecondaries==0)
        error('Not all secondary systems are online. Please run "multiCheck(1)" for complete diagnostics');
    end
end
switch numSecondaries
    case 4
        IPlist = {'10.10.0.2', '10.10.1.2', '10.10.2.2', '10.10.3.2'};  % 5 system configuration
    case 3
        IPlist = {'10.10.1.2', '10.10.2.2', '10.10.3.2'}; % 4 system configuration
    case 2
        IPlist = {'10.10.1.2', '10.10.2.2'}; % 3 system configuration
    case 1
        IPlist = {'10.10.1.2'}; % 2 system configuration
    otherwise
        error('The system network IP addresses are in a non-standard configuration.  Please contact Verasonics support.')
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

%--- Initialize Tmux Panels ---%
tmuxPID = vsv.multi.util.initTmux(flags, IPlist);
if (tmuxPID == 0)
    error('Could NOT initalize Tmux windows ')
end

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

% close opened terminal
%unix(sprintf('kill -9 %i', tmuxPID));