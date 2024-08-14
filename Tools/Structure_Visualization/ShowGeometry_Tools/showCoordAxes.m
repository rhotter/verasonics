function showCoordAxes (figHandle)
    % showCoordAxes (figHandle)  plots the X,Y,Z coordinate system into an existing figure or into figure 50.
    %
    % 	The function provides a graphical representation of the Verasonics coordinate system
    %  	by plotting lines through the origin in 3D
    %
    %   See also: showGeometry, showPData, showTrans, showMediaPts

    % PK -- beta version 11-JUN-2018
    % YT -- refined on 29-MAR-2019

    persistent CoorPlot
    verbose = evalin('base', 'Resource.Parameters.verbose');
    if verbose>2
        disp([mfilename ' Version 3-Apr-2018'])
    end

    %% Plot will go into figure handle argument. If none provided, use #50.
    if nargin==0
        VsClose
        figHandle = figure(50);
        set(gcf,'Position', [500 100 800 800]);
    end

    %% Draw the axes with labels on postive ends
    figure(figHandle)
    rotate3d on
    grid on

    view(-35, 15);
    DkBl = [0 0.4 0.6]; % Dark Blue / Teal color

    A = axis;
    if ishandle(CoorPlot)
        delete(CoorPlot); % Delete old CoorPlot for multitple Pdata structures
    end

    hold on
    % Draw X axis
    CoorPlot(1) = plot3(A(1:2), [0 0], [0 0], 'linewidth', 3, 'color', DkBl );
    CoorPlot(2) = text( A(2), 0, 0, ' +X', 'fontsize', 24, 'fontweight', 'bold','color', DkBl);
    % Draw Y
    CoorPlot(3) = plot3([0 0], A(3:4), [0 0], 'linewidth', 3, 'color', DkBl );
    CoorPlot(4) = text( 0, A(4), 0, ' +Y', 'fontsize', 24, 'fontweight', 'bold', 'color', DkBl);
    % Draw Z
    CoorPlot(5) = plot3([0 0], [0 0], [0 A(6)], 'linewidth', 3, 'color', DkBl );
    CoorPlot(6) = text( 0, 0, A(6), ' +Z', 'fontsize', 24, 'fontweight', 'bold', 'color', DkBl);

    set(CoorPlot,'HandleVisibility','off'); % if necessary, cla in won't clear the axis
    if nargin==0
        set (gca, 'ZDir', 'reverse','Ydir', 'reverse')
    end
    hold off

end
