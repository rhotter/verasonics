%MEXAUDIOSTREAM  PortAudio functions, interfaced through MatLab.
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.
% Help for mexaudiostream.c
%
%    >> mexaudiostream(0)         % open audio stream object
%    >> mexaudiostream(1000,x)    % load fifo with interleaved stereo (DOUBLE)
%                                 % samples in row vector x
%    >> mexaudiostream(1)         % start playing
%    >> mexaudiostream(2)         % stop playing
%    >> isopen = mexaudiostream(3) % stream is open
%    >> mexaudiostream(-1)        % destroy stream object
%    >> x = mexaudiostream(901,N) % FIFO read - test/diagnostics only - Do Not Use while
%                                 % player object is running.
%    >> mexaudiostream(902, x)    % Write row vector (DOUBLE class) to FIFO,
%                                 % format is interleaved samples Left/Right.
%    >> % note: in the following, one L/R pair counts as two data elements (two samples).
%    >> Nocc = mexaudiostream(903) % FIFO Occupancy in samples.
%    >> Nvac = mexaudiostream(904) % FIFO Vacancy in samples.
%    >> Ncap = mexaudiostream(905) % FIFO Capacity in samples.
%    >> framecount = mexaudiostream(906) % frames consumed since stream opened.
%    >> isRunning = mexaudiostream(907) % indicates if stream is running.
%    >> FsDAC = mexaudiostream(908) % DAC Sample Rate in Hz
%    >> wasEmpty = mexaudiostream(909) % Fifo ran out of data (but may have gotten more later).
%    >> minRdThrshFlt = mexaudiostream(910) % Fifo min read amount threshold - must have this many floats on fifo to read any.
%    >> DACFrameSize = mexaudiostream(920) % size in stereo samples of the
%    DAC frame, i.e. the read event quanta.

%   MEX-File function.
