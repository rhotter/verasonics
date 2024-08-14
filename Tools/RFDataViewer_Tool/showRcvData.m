function showRcvData (varargin)
% showRcvData: Stand-alone RF data viewer; run showRcvData command after quitting VSX.
%
% showRcvData is a stand-alone utility to view the RF data in the RcvData buffer using a grayscale image and line plots.
%   When no arguments are provided, the data displayed uses default settings, listed below.
%   Each of the following variables can be set in the function call using the Matlab convention:
%   showRcvData( 'variableName', numericalValue, ...)
%   with:
%       variableName    [default]   Description
%       ------------    ----------- -----------------------------------------------------------------------------------------
%       bufferNum       [1]       	Receive Buffer index
%       frameNum        [1]        	Frame number (can only plot one frame at a time)
%       acqNum          [0]         Acquisition number: =0 All  =N Nth acquisition in frame 'frameNum'
%       chNums          [32,96] 	channel numbers to plot as lines (=0 skips all line plotting) (values for 128 channels)
%       interpF         [4]        	interpolation upsampling factor for the line plots (=1 for no interpolation)
%       plotSamplePoints [0]       	=1 will plot '+' at each true data point, and 'o' at each interpolated sample, w/o line
%       cMax            [500]       absolute value of data above which color will be saturated black or white in the 2D image
%                                   and scale line plots to [-2*cMax 2*cMax]
%       ordering        ['channel'] 'channel':  plot data in channel order, exactly as stored in RcvData structure
%                                   'element':  plot data in element order, using Trans.Connector vector to unscramble the RF data
%       displayType     ['RF']      'magnitude': plots the envelope of the RF data on black background
%                                   'RF' (or anything other than 'magnitude') will plot the usual signed data on grayscale
%
%   By default, the program plots an entire frame, and marks the boundaries between acquisition events in the frame.
%   The zoom feature is turned on by default, so one can use the figure's rubber band box zoom control
%   to examine data details for a smaller region within a frame.
%   If the sampleRange is specified explicitly, the frameNum is ignored, and the plotted range corresponds to the samples in the
%   receive buffer. In this case, the acquisition boundaries are not marked.
%   If an acquisition number is provided, that acquisition (for frame=frameNum) is plotted. The data is also interpolated by
%   interpF to make the image smoother.
%   If ordering is by element, then chNums is used to represent element number(s) to plot.
%
%   This version of the function is intended ot be used from the Matlab interface after running VSX. A version suitable for
%   use as an external function, which requires that the only input argument be a local RF data buffer variable name, can be
%   developed from this function.
%
%   Known Bugs:
%       - Does not work with samplesPerWave = 1, 2, 4/3, or 8/3
%       - Does not work with interleaved sampling
%       - no separate documentation file (these help comments are the only available documentation)

% Testing: Tested with software release 3.2.0 on Vantage 256
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Update History:
%  6-JUN-2015   - added showRcvData to Tools folder
% 10-JUL-2015   - added use of Trans.connector to unscramble RcvData according to elment to channel mapping
%  4-Nov-2016   - added the option to plot a single acquisition by providing the acqNum
% 20-Mar-2017   - added support for other sampling modes; removed the ability to specify a sample range...can plot one frame
%                   or one acquisition, and zoom as desired.
%  8-JUN-2018 - added automatic scaling of cMax; revised figure sizing to larger defaults
% 13-JUN-2018 - added the displayType 'magnitude'  to be able to plot the envelope instead of the signed RF field.


%% Read in the required structures
if evalin('base', 'exist(''RcvData'', ''var'')')
    RcvData = evalin( 'base', 'RcvData' );
else
    error( [mfilename ': RcvData structure is required and is not currently in the workspace'])
end

if evalin('base', 'exist(''Resource'', ''var'')')
    Resource = evalin( 'base', 'Resource' );
else
    error( [mfilename ': Resource structure is required and is not currently in the workspace'])
end

if evalin('base', 'exist(''Receive'', ''var'')')
    Receive = evalin( 'base', 'Receive' );
else
    error( [mfilename ': Receive structure is required and is not currently in the workspace'])
end

% if evalin('base', 'exist(''Trans'', ''var'')')
% Trans = evalin( 'base', 'Trans' );
% else
%     error( [mfilename ': Trans structure is required and is not currently in the workspace'])
% end


%% Read in any input variable pairs in the 'matlab' convention, as a name/value pair.
% 'name' must be a string, and value is a self-defined type. E.g., showRcvData( 'bufferNum', 2 )
len = length(varargin);
for ii=1:2:len
    if strcmp(varargin{ii},'ordering')==1 || strcmp(varargin{ii},'displayType')==1
        str = [ varargin{ii} '= ''' varargin{ii+1} ''';' ];
    else
        str= ([varargin{ii} ' = [' num2str(varargin{ii+1}) '];']);
    end
    eval (str);
end

%% Validate inputs and set default RcvData parameters

if exist( 'bufferNum', 'var')
    numBufs = max(size(RcvData));
    if bufferNum>numBufs
        error([mfilename ': *** Variable bufferNum = ' num2str(bufferNum) ' exceeds the number of Receive buffers = ' num2str(numBufs)])
        return
    end
else
    bufferNum = 1;
end

if exist( 'frameNum', 'var')
    if frameNum > Resource.RcvBuffer.numFrames
        error([mfilename ': *** Variable frameNum = ' num2str(frameNum) ' exceeds the number of Receive frames = ' num2str(Resource.RcvBuffer.numFrames)])
        return
    end
else
    frameNum = 1;
end

acqsPerFrame = size(Receive, 2)/Resource.RcvBuffer(bufferNum).numFrames;    % Note: Receive is a 1xN structure,
                                                                            %   where N is the total number of acquisitions.
acq0 = (frameNum-1)*acqsPerFrame;   % zeroth acquisition
sampleOffset = Receive(acq0+1).startSample - 1;

if exist( 'acqNum', 'var')
    if acqNum > acqsPerFrame
        error([mfilename ': *** Variable acqNum = ' num2str(acqNum) ' exceeds the number of Receive acquisitions = ' num2str(acqsPerFrame)])
        return
    elseif acqNum < 0
        error([mfilename ': *** Variable acqNum = ' num2str(acqNum) ' cannot be less than zero.'])
        return
    else
        startSample = Receive(acqNum +acq0).startSample;
        endSample = Receive(acqNum +acq0).endSample;
    end
else
    acqNum = 0;
    startSample = Receive(acq0+1).startSample;
    endSample = Receive(acq0+acqsPerFrame).endSample;
end
sampleRange = startSample:endSample;
numSamples = length(sampleRange);

if ~exist( 'ordering', 'var')
    ordering='channel';   % set default to plot RF data as stored in RcvData array
end

% if exist( 'sampleRange', 'var')
%     if max(sampleRange) > Resource.RcvBuffer(bufferNum).rowsPerFrame
%         error([mfilename ': *** Variable sampleRange = ' num2str(sampleRange(end)) ' exceeds Resource.RcvBuffer(' num2str(bufferNum) ').rowsPerFrame  = ' num2str(Resource.RcvBuffer(bufferNum).rowsPerFrame)])
%         return
%     elseif min(sampleRange) < 0
%         error([mfilename ': *** Variable sampleRange = ' num2str(sampleRange(end)) ' cannot be less than zero.'])
%         return
%     end
%     numSamples = length(sampleRange);
%
% else
%     % find the number of samples in a single frame from the Receive structure
%     for rcv = 1:size(Receive, 2)
%         frame = Receive(rcv).framenum;
%         if frame == frameNum, break, end
%     end
%
%     if acqNum == 0
%         rcv = (frameNum-1)*acqsPerFrame +1;     % index of first Receive structure in the desired frame
%         endSample = Receive(rcv+acqsPerFrame-1).endSample;
%     else
%         rcv = (frameNum-1)*acqsPerFrame + acqNum;     % index of acqNum Receive structure in the desired frame
%         endSample = Receive(rcv).endSample;
%     end
% end

% determine the sample indices for all of the acquisition boundaries
acqBndry = [];
for acq = 1:acqsPerFrame
    acqBndry = [ acqBndry, Receive(acq0+acq).endSample ];
end

%% Verify inputs for plot parameters and set defaults
if exist( 'chNums', 'var')
    if max(chNums) > Resource.RcvBuffer(bufferNum).colsPerFrame
        error([mfilename ': *** One or more chNums = [' num2str(chNums) '] exceeds the number of Receive columns = ' num2str(Resource.RcvBuffer(bufferNum).colsPerFrame)])
        return
    end
else
    chNums = round([.25 .75]*Resource.RcvBuffer(bufferNum).colsPerFrame);
end

if exist( 'interpF', 'var')
    if interpF<1
        error([mfilename ': *** Parameter interpF = ' num2str(interpF) ' must be an integer greater than or equal to 1'])
        return
    end
else
    interpF = 4;            % interpolation upsampling factor for the line plots
end

if ~exist( 'plotSamplePoints', 'var')
    plotSamplePoints = 0; 	% =1 will plot '+' at each true data point
end

if ~exist( 'cMax', 'var')
    cMax = 0;               % absolute RcvData count above which color will be saturated black or white
end

if ~exist( 'displayType', 'var')
    displayType = 'RF';        % choices are 'magnitude' or anything else
end


%% Extract the RF data from the buffer
clear Rdata Ridata
channelRange = 1:size(RcvData{bufferNum}, 2);
Rdata1 = double(RcvData{bufferNum}(:, channelRange, frameNum)); % extract entire frame of data from the buffer
Rdata = Rdata1;
if strcmp (ordering, 'element')==1      % --- Unscramble the data from I/O channel ordering into element ordering
    % read in Trans for Trans.Connector
    if evalin('base', 'exist(''Trans'', ''var'')')
        Trans = evalin( 'base', 'Trans' );
    else
        error( [mfilename ': Trans structure is required and is not currently in the workspace'])
    end

    % unscramble the data
    for chan = channelRange
        Rdata(:,chan) = Rdata1(:, Trans.Connector(chan));
    end
end

%% Convert all sample modes to 4x sampling for plotting
% Because only the 4x sampling mode produces equally spaced samples in time, the other modes need to be processed
%   before plotting, or the data will look distorted. Users can use these approaches to reconstruct the uniformly sampled
%   RF data for other applications requiring signal processing prior to image reconstruction.
%
% The original RF data in the RcvData buffer is described in the Receive structure for the Receive(i).sampleMode used in each
%   acquisition. To simplify plotting of the RF data in the buffer, we convert all data in the desired frame to 4x sampling.
%   This means that the points in the plot will now not correspond to the sample numbers in the RcvData buffer, but the
%   y-scale (y-axis tick marks) and acquisition boundary lines should still correctly refer to the RF samples
%   in the receive buffer. This may be confusing when different sampling scehems are used in one buffer, such as for
%   interleaved B-mode and Doppler acquisitions. The new 4x buffer variables will be denoted with the 4x in the variable
%   name.

% step through all acquisitions and create a new RF buffer for the 4x data
Rdata4x = [];
if acqNum==0
    acqs = 1:acqsPerFrame;
else
    acqs = acqNum;
end

for acq = acqs
    rcv = acq + acq0;
    sS = Receive(rcv).startSample - sampleOffset;  % samples relative to that frame only
    eS = Receive(rcv).endSample - sampleOffset;
    Srange = sS:eS;
    switch Receive(rcv).sampleMode
        case 'NS200BW'      % Nyquist sampling (200% bandwidth) of demodFrequency.
            samplesPerWave = 4;
            % Rdata is already at 4x so just keep it.
            Rdata4x = [Rdata4x; Rdata(Srange,:)]; %#ok<*AGROW>

        case 'NS200BWI'     % 2:1 interleaved sampling (requires 2 acquisitions)
            % the data is transferred un-interleaved; interleaving is done as part of recon operation
            error ( [mfilename ' *** Receive(',num2str(rcv),').sampleMode = ''NS200BWI'' is not supported yet. '])
            return

        case 'BS100BW'      % 100% bandwidth sampling of demodFrequency.
            samplesPerWave = 4;
            for chan = 1:size(Rdata, 2)
                I(:,chan) = Rdata(sS  :2:eS-1, chan);   % pick out the odd samples
                Q(:,chan) = Rdata(sS+1:2:eS,   chan);   % even samples
                I2x(:,chan) = interp(I(:,chan), 2);  	% interpolate by a factor of 2
                Q2x(:,chan) = interp(Q(:,chan), 2);    	% ditto
                RF4x(:,chan) = upsample(I2x(:,chan),2)  + upsample(Q2x(:,chan),2, 1); % insert zeros and sum
            end
            Rdata4x = [Rdata4x; RF4x]; %#ok<*AGROW>

        case 'BS67BW'       % 67% bandwidth sampling of demodFrequency.
             error ( [mfilename ' *** Receive(',num2str(rcv),').sampleMode = ''BS67BW'' is not supported yet. '])
            return

        case 'BS50BW'       % 50% bandwidth sampling of demodFrequency.
            error ( [mfilename ' *** Receive(',num2str(rcv),').sampleMode = ''BS50BW'' is not supported yet. '])
            return

        case 'custom'       % sample rate set by decimSampleRate
            error ( [mfilename ' *** Receive(',num2str(rcv),').sampleMode = ''custom'' is not supported yet. '])
            return

        otherwise
            error ([mfilename, ' *** Unrecognized sampleMode: ' Receive(rcv).sampleMode ' for index ' num2str(rcv)])
            return
    end

end
%switch
%     if ~strcmp(Receive(rcv).sampleMode,'NS200BW')
%         disp( [' *** Currently, only Receive.sampleMode = ''NS200BW'' is supported by ' mfilename '.' ])
%         disp( '     Please see Vantage Programming manual for a description of Receive.sampleMode to understand how the RcvData samples are stored.')
%     else
%         samplesPerWave = 4;
%     end

%% Plot the RF data as an image
hR =figure(20); % Image of RcvData portion
    GR = groot; hR.Position = [ 1, 1, floor(GR.ScreenSize(3)/3), floor(.8*GR.ScreenSize(4))];
    if cMax==0 % autoscale if no input
        imMax = max(max(Rdata4x(:, channelRange)));
        cMax = max(cMax, imMax);
        cMax = min( cMax, 8000 ); % maximum permitted value of cMax
    end

    if strcmp(displayType, 'magnitude')
        Idata4x = envelope (Rdata4x);
        imagesc(channelRange, sampleRange, Idata4x); caxis([  0   cMax]); colormap gray(32); zoom on
    else
        imagesc(channelRange, sampleRange, Rdata4x); caxis([-cMax cMax]); colormap gray(64); zoom on
    end
    if acqNum==0
        title ( ['Rcv Buffer: ' num2str(bufferNum) '     Frame Num: ' num2str(frameNum) '     Acq Num: All     cMax: ' num2str(cMax)] )
    else
        title ( ['Rcv Buffer: ' num2str(bufferNum) '     Frame Num: ' num2str(frameNum) '     Acq Num: ' num2str(acqNum) '     cMax: ' num2str(cMax)] )
    end
    hold on
    if exist('acqBndry', 'var') && acqNum==0
        for n=1:acqsPerFrame % plot the acquisition event boundaries
            plot ([channelRange(1), channelRange(end)], acqBndry(n)*[1,1], 'k--')
            text( channelRange(10), acqBndry(n)-50, ['Acq ' num2str(n)])
        end
    end
    if chNums~=0
        Xchannels = [ chNums(1) chNums(1) ];    % array of line coordinates to draw onto image
        Ychannels = [sampleRange(1) sampleRange(end)];
        plot(Xchannels', Ychannels', '--', 'LineWidth', 2)
        Ylim = get(gca, 'Ylim');
        text( chNums(1)-1, Ylim(2)+(Ylim(2)-Ylim(1))/100, num2str(chNums(1)), 'FontSize', 10 )
        Lstring =  sprintf('% 4i', chNums(1));                 % legend strings
        if length(chNums)>1
            for ii=2:length(chNums)
                Xchannels = [chNums(ii) chNums(ii)];
                Ychannels = [sampleRange(1) sampleRange(end)];
                plot(Xchannels', Ychannels', '--', 'LineWidth', 2)
                Ylim = get(gca, 'Ylim');
                text( chNums(ii)-1, Ylim(2)+(Ylim(2)-Ylim(1))/100, num2str(chNums(ii)), 'FontSize', 10 )
                Lstring = [Lstring; sprintf('% 4i', chNums(ii))];                                  %#ok<AGROW>
            end
        end
    end
    hold off
    if strcmp(ordering,'element')==1, xlabel ('Element Number'), else xlabel ('Channel Number'), end
    ylabel ('Sample Number')

%% Plot individual channel data (interpolated by interpF)
% Note that if ordering is by element, chNums now represents element numbers to plot
if chNums~=0
    hChan = figure(21);
    hChan.Position = [ floor(GR.ScreenSize(3)/3) floor(GR.ScreenSize(4)/3) floor(GR.ScreenSize(3)*3/5) floor(GR.ScreenSize(4)/3)]; % scale to screen size
    for ii=1:length(chNums)
        Ridata(:,ii) = interp( Rdata4x(:,chNums(ii)), interpF ); %#ok<AGROW>
    end
    X = ((1/interpF:1/interpF:numSamples) + (interpF-1)/interpF + sampleRange(1)-1)';
    wavePerSample = 2*samplesPerWave;
    Xwv = X/wavePerSample;
    %     [size(X) size(Ridata)]
    if plotSamplePoints==1
        plot(sampleRange, Rdata4x(:,chNums), '+'), hold on, plot(X, Ridata(:,:), 'o'); hold off,
    else
        if strcmp(displayType, 'magnitude')
            plot(X, Ridata(:,:)), hold on, ax = gca; ax.ColorOrderIndex=1; plot(X, envelope(Ridata(:,:)), 'linewidth', 2), hold off
        else
            plot(X, Ridata(:,:));
        end
        save ('plotData', 'Xwv', 'Ridata')
    end
    hold on
    if exist('acqBndry', 'var')
        if acqNum==0
            for n=1:acqsPerFrame % plot the acquisition event boundaries
            plot ( acqBndry(n)*[1,1], cMax*[-2,2], 'k--')
            text( acqBndry(n)-50, -1.8*cMax, ['Acq ' num2str(n)], 'Rotation', 90)
            end
        else
            Xlim = get(gca, 'Xlim');
            text( Xlim(2)-numSamples/20, -1.8*cMax, ['Acq ' num2str(acqNum)], 'Rotation', 90)
        end
    end
    hold off
     if strcmp(ordering,'element')==1
         title (['RF Data Interpolated ' num2str(interpF) 'x  -- Elements ' num2str(chNums)])
     else
         title (['RF Data Interpolated ' num2str(interpF) 'x  -- Channels ' num2str(chNums)])
     end
    legend( Lstring ), zoom on
    ylabel ('Receive Data (interpolated)'), xlabel ('Samples')
    axis([sampleRange(1) sampleRange(end) -2*cMax 2*cMax]); % comment out this line to autoscale each plot
end

return
