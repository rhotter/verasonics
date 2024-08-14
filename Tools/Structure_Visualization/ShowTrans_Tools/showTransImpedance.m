function showTransImpedance( figHandle )
    % showTransImpedance plots the complex impedance vs. freqeucy for a particular transducer as provided in the Trans structure
    %
    %  	See also: showGeometry, showPData, showTrans, showMediaPts, showCoordAxes, showTransIR


    verbose = evalin('base', 'Resource.Parameters.verbose');
    if verbose>2
        disp([mfilename ' Version 3-Apr-2018'])
    end

    if evalin('base', '~exist(''Trans'',''var'')')
        error ('Structure "Trans" does not exist. Please define the transducer and run again.')
     else
         Trans = evalin('base', 'Trans');    % load the Trans structure
     end

   if nargin    % if no figure handle provided, create a new figure
        figure(figHandle);
    else
        figure
    end


    % --- consider the case where impedance is only specified as a scalar
    expandAxis = [0 0 0 0];
    if size(Trans.impedance) == [1 1]
        Trans.impedance = [ Trans.Bandwidth(1), Trans.impedance ;  Trans.Bandwidth(2), Trans.impedance];
        expandAxis = [-1 1 0 50];
    end

    % --- plot the complex impedance as a function of frequency (note the frequencies are also complex)
    plot(abs(Trans.impedance(:,1)), abs(Trans.impedance(:,2)), 'linewidth', 3)
    hold on
    plot(abs(Trans.impedance(:,1)), real(Trans.impedance(:,2)), '--')
    plot(abs(Trans.impedance(:,1)), imag(Trans.impedance(:,2)), '--')

    A = axis; A = A+expandAxis; axis(A);
    plot( Trans.frequency*[1 1], A(3:4), 'k' )
    plot( Trans.Bandwidth(1)*[1 1], A(3:4), 'color', [0.8 0.8 0.8], 'linewidth', 2 )
    plot( Trans.Bandwidth(2)*[1 1], A(3:4), 'color', [0.8 0.8 0.8], 'linewidth', 2 )
    hold off

    legend ('Magnitude', 'Real', 'Imaginary', 'Trans.frequency', 'Fmin', 'Fmax', 'Location', 'southeast')
    xlabel( 'Frequency (MHz)')
    ylabel( 'Impedance (Ohms)')
    title ([ 'Element Impedance for ' Trans.name ], 'interpreter', 'none')

end

