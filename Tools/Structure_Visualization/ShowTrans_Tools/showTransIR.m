function showTransIR( figHandle )
    % showTransIR (figHandle)   plots the transducer element's Impulse Response (IR) and passband response (spectrum)
    % 	for a particular transducer defined in the current Trans structure
    %
    %  	See also: showGeometry, showPData, showTrans, showMediaPts, showCoordAxes, showTransImpedance

    verbose = evalin('base', 'Resource.Parameters.verbose');
    if verbose>2
        disp([mfilename ' Version 3-Apr-2018'])
    end

    %% Read in the Trans structure and open the figure
    if evalin('base', '~exist(''Trans'',''var'')')
        error ('Structure "Trans" does not exist. Please define the transducer and run again.')
    else
        Trans = evalin('base', 'Trans');    % load the Trans structure
    end

    if nargin    % if no figure handle provided, create a new figure
        figure(figHandle);
    else
        fh = figure (50); fh.Position=[20 20 800 600];
    end


    %% Compute and plot the spectrum from the IR if specified, or determine from the bandwidth of the transducer

    Fs = 250;   % the IR waveforms are sampled at the master clock freqeuncy of 250 MHz

    if ~isfield(Trans,'Bandwidth') % default to 75% bandwidth
        BW = 0.75; Trans.Bandwidth = Trans.frequency*[1-BW/2, 1+BW/2];
        assignin('base','Trans',Trans);
    end

    Fmin = Trans.Bandwidth(1);
    Fmax = Trans.Bandwidth(2);
    BW = Fmax-Fmin;

    % --- plot the transducer's spectral passband reponse
    if ~isfield( Trans, 'IR1wy' ) % compute the Spectrum using a 2nd order Butterworth filter and the transducer bandwidth
        [bflt,aflt]=buttervs(2, Trans.Bandwidth*2/Fs);
        [H,F] = freqz(bflt, aflt,1024,Fs);
    else
        [H,F] = freqz(Trans.IR1wy,1, 1024,Fs);
    end

    subplot (2,1,1)
        Ha = abs(H); Ha = Ha/max(Ha);
        plot(F, Ha, 'linewidth', 3)
        hold on
        plot(F, Ha.^2, 'linewidth', 3)

        A = axis; A = [max(Fmin-BW,A(1))  min(Fmax+BW,A(2)),  A(3:4)]; axis(A)
        plot( Trans.frequency*[1 1], A(3:4), 'k' )
        plot( Fmin*[1 1], A(3:4), 'color', [0.8 0.8 0.8], 'linewidth', 2 )
        plot( Fmax*[1 1], A(3:4), 'color', [0.8 0.8 0.8], 'linewidth', 2 )
        hold off

        legend ('1-way', '2-way', 'Trans.frequency', 'Fmin', 'Fmax', 'Location', 'northeast')
        xlabel( 'Frequency (MHz)')
        ylabel( 'Spectrum')
        title ([ 'Spectrum for ' Trans.name ], 'interpreter', 'none')

    % --- plot the impulse responses as a function of time, on our 1/(250 MHz) grid.
    if ~isfield( Trans, 'IR1wy' ) % compute the IR using a 2nd order Butterworth filter and the transducer bandwidth
        disp(['  *** ' mfilename ': WARNING: please run VSX first to populate the Impulse Response field in Trans ***'])
        disp( ' ' )
        return
    else
        IR1wy = Trans.IR1wy;
        IR2wy = Trans.IR2wy;
    end

    subplot (2,1,2)
        t = 0:length(Trans.IR1wy)-1; t = t/Fs;
        plot(t, IR1wy, 'linewidth', 3)
        hold on
        plot(t, IR2wy, '--', 'linewidth', 3)
        hold off
        A=axis; if A(3)~=A(4), A(4)=max(abs(A(3:4))); A(3)=-A(4); axis( A ); end % center the IR plot
        legend ('1-Way', '2-Way')
        xlabel( 'Time (microsec)')
        ylabel( 'Amplitude')
        title ([ 'Impulse Responses for ' Trans.name ], 'interpreter', 'none')

    %     A = axis; A = A+expandAxis; axis(A);


end
