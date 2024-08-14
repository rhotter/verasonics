function showGeometry (figHandle)
    % showGeometry (figHandle) plots three structures in the Verasonics coordinate system,
    % 	including MediaPts, PData, and Trans in three dimensions.
    %  	If no argument is provided, a new figure will be created.
    %  	The user is encouraged to use the MATLAB window menu
    %  	to rotate and zoom the plot as desired.
    %
    %   See also: showPData, showTrans, showMediaPts, showCoordAxes

    % PK -- beta version 11-JUN-2018
    % YT -- refined on 29-MAR-2019

    %% Plot will go into figure handle argument. If none provided, create a new one.
    if nargin==0
        VsClose
        figHandle = figure;
        GR = groot; height = floor(.75* GR.ScreenSize(4)); ystart = GR.ScreenSize(4)-height;
        set(gcf,'Position', [200 ystart height height]);
    end

    %% plot all three strutures identifying the physical geometry
    showTrans (figHandle);        hold on
    showPData (figHandle);        hold on
    showMediaPts (figHandle);     hold on

    %% modify figure properties
    showCoordAxes (figHandle);

    axisname = 'wavelengths';
    if evalin('base','exist(''Resource'',''var'')')
        Resource = evalin ('base', 'Resource');
        if isfield(Resource.DisplayWindow,'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
            if strcmpi(Resource.DisplayWindow(1).AxesUnits,'mm')
                axisname = 'mm';
            end
        end
    end

    xlabel (['X (' axisname ')']), ylabel (['Y (' axisname ')']), zlabel (['Z (' axisname ')']),
    name = 'undefined'; % if user did not provide a name set to 'undefined' by default  

    if evalin('base', 'exist(''displayWindowTitle'')')
        name = evalin('base', 'displayWindowTitle');
    else
        if evalin('base', 'exist(''Trans'')')
            Trans = evalin('base','Trans');
            if isfield(Trans, 'name') && ~isempty(Trans.name)
                name = evalin('base', 'Trans.name');
            end
        end
    end
    title (['Geometry for ' name], 'interpreter', 'none')

    set (gca, 'ZDir', 'reverse')
    set (gca, 'Ydir', 'reverse')

    hold off

end
