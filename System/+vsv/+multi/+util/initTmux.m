function tmuxPID = initTmux(flags, IPlist)

%vsv.multi.util.initTmux
%  Utility for initializing the multipanel terminal (TMUX) for communicating
%  with all secondary systems for a multisystem configuration
%
%  returns the tmux process ID (process) for killing window at the end
%  returns a zero if a PID could not be obtained.
% Usage:
%  typical usage:
%  tmuxProcessID = vsv.multi.util.initTmux('-Tg', {'10.10.1.2', '10.10.2.2', '10.10.3.2'});
%  for testing/debugging:
%  tmuxProcessID = vsv.multi.util.initTmux('-Tg', {'localhost','localhost','localhost'});
%
% Copyright (C) 2001-2023, Verasonics, Inc.  All worldwide rights and
% remedies under all intellectual property laws and industrial property
% laws are reserved.

arguments
    flags string
    IPlist string
end

numSecondaries = length(IPlist);
%--- Close Dead Terminals from previous runs ---%
[result, resultStr] = unix('tmux kill-session -t multisys'); %kill existing terminal if exists
if result > 0 && (~contains(resultStr, '/tmp/tmux-1000/default (No such file or directory)')) && (~contains(resultStr,'no server running on /tmp/tmux-1000/default'))
    error('initTmux:closeSession','Could not close previous tmux session.\nError: %s\n Try closing all open windows and restarting MATLAB.', resultStr);
end

%-- Start new terminal --%
% the option is available to decide between xterm and gnome terminal. but
% it defaults to gnome terminal.  may eventually deprecate xterm
if contains(flags, '-Tx')
    [status, resultStr] = unix('xterm -maximize -e "tmux new-session -s multisys -n ''MultiSys'' " &');  %to run on Xterm instead of gnome-terminal
    if (status>0) && strfnd('duplicate',result)
        unix('xterm attach -t multisys')
    end
elseif contains(flags, '-Tg')||~contains(flags,'-T')
    [status, resultStr] = unix('gnome-terminal --window --maximize -- bash -c "tmux new-session -s multisys -n ''MultiSys'';exec bash"');
end
if (status > 0 )
    error('initTmux:initializeTerminal','Could not initialize a new terminal.\nError: %s\n Try closing all open windows and restarting MATLAB.', resultStr);
end
%wait for tmux session
for a = 1:10
    if (isfile('/tmp/tmux-1000/default'))
        break
    end
    if a == 10
        error('initTmux:initializeTerminal','Could not initialize a new terminal.')
    end
    pause(0.1)
end

%-- Setup Tmux Window --%
switch numSecondaries
    case 4
        [status, resultStr] = unix('tmux split-window -v \; split-window -h \; selectp -t 0 \; split-window -h'); %4 panes
    case 3
        [status, resultStr] = unix('tmux split-window -v -p 66 \; split-window -v -p 50'); %3 horizontal panes
    case 2
        [status, resultStr] = unix('tmux split-window -v -p 50'); %2 horizontal panes
    case 1
        %place holder since the tmux window does not need to be split
    otherwise
        error('initTmux:incorrectNumberOfSecondaries','%i secondaries are not supported by the Verasonics MultiSystem software.', numSecondaries);
end
if (status >0 )
    error('initTmux:generateTmuxWindows','Could not generate Tmux windows to communicate with secondary systems.\nError:%s\nTry Closing all open windows and restarting MATLAB.', resultStr);
end

% Get PID for later;
tmuxPID = getTmuxPID();

%-- Load in tmux settings --
[status, resultStr] =unix(['tmux source-file ' vsv.file.getVSXDir() '/System/+vsv/+multi/tmux.conf']);
if (status > 0)
    error('initTmux:loadTmuxdefaults','Could not load Tmux settings.\nError:%s\nTry Closing all open windows and restarting MATLAB.', resultStr);
end

% 3. connect to each secondary
for a = 1:numSecondaries
    cmd = sprintf('tmux send-keys -t multisys:0.%i ''ssh %s'' Enter', a, IPlist{a});
    [status, resultStr] = unix(cmd);
    % test result @todo cannot actually see if a no route to host is
    % returned
    if (status > 0)||contains(resultStr,'No route to host')
        error('initTmux:connectToRemote','Could not connect to remote system through Tmux.\nError:%s\nTry Closing all open windows and restarting MATLAB.', resultStr);
    end
end

% 4. on all secondaries
[status, resultStr] = unix('tmux setw synchronize-panes on'); %multicast on
if (status > 0)
    error('initTmux:enableMulticast','Could not enable Tmux multicast.\nError:%s\nTry Closing all open windows and restarting MATLAB.', resultStr);
end

% 5. cd to workspace
[status, resultStr] = unix(sprintf('tmux send-keys -t multisys:0 "cd %s" Enter', vsv.file.getVSXDir()));
if (status > 0)
    error('initTmux:changeDirectory','Could not change to distribution directory.\nError:%s\nPlease run vsv.multi.configCheck to verify proper installation', resultStr);
end

% Clear any existing matlab processes on secondary systems
if (vsv.multi.util.isSecondarySystem())
    unix('tmux send-keys -t multisys:0 ''killall -9 MATLAB'' Enter');
end

% Start Matlab
[status, ~] = unix('tmux send-keys -t multisys:0 ''matlab -nosplash -nodesktop -r "run activate.m"'' Enter'); %normal
if (status>0)
    error('initTmux:startTmuxSession','Could not start tmux session.  Try closing all previous tmux terminals and retry')
end
pause(7); %ugly fix. the line above needs to really block until matlab has started.

unix('tmux setw synchronize-panes off'); %multicast off
end

function tmuxPID = getTmuxPID()
% Get the tmux PID but return with a zero if the function fails

[~, returnStr] = unix('ps -ef |grep tmux');
% Regular expression pattern to match the line containing "bash -c tmux new-session -s multisys -n 'MultiSys';exec bash"
pattern = ".*bash -c tmux new-session -s multisys -n 'MultiSys';exec bash.*";

% Match the pattern in the example string using regexp
match = regexp(returnStr, pattern, 'match');

% Extract the second integer value from the matched line using regexp
pid_pattern = "\d+";
[pid_match, ~] = regexp(match, pid_pattern, 'match', 'once', 'tokenExtents');
if ~isempty(pid_match)
    % Convert the extracted string to a numeric value
    tmuxPID = str2double(pid_match);
else
    warning('initTmux:tmuxPID','Could not obtain Tmux PID');  
    tmuxPID = 0;
end
end