function saveRF(varargin)
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   tool for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: saveRF.m - A tool to save RF
%

if ~isempty(findobj('tag','UI')) % running VSX
    if evalin('base','freeze')==0   % no action if not in freeze
        msgbox('Please freeze VSX');
        return
    else
        Control.Command = 'copyBuffers';
        runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
    end
else % not running VSX
    if evalin('base','exist(''RcvData'',''var'');')
        RcvData = evalin('base','RcvData');
    else
        disp('RcvData does not exist!');
        return
    end
end

RFfilename = ['RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];

RcvLastFrame = size(RcvData,3);
if (~evalin('base','simButton'))
    RcvLastFrame = Resource.RcvBuffer(1).lastFrame;
end

[fn,pn,~] = uiputfile('*.mat','Save RF data as',RFfilename);
if ~isequal(fn,0) % fn will be zero if user hits cancel
    fn = strrep(fullfile(pn,fn), '''', '''''');
    save(fn,'RcvData','RcvLastFrame','-v7.3');
    fprintf('The RF data has been saved at %s \n',fn);
else
    disp('The RF data is not saved.');
end

end
