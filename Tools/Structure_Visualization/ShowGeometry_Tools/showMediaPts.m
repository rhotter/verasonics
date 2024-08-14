function showMediaPts (figHandle)
% showMediaPts (fighandle)  plots the point scatterers defined in the Media structure in three dimensions
%
%   The function provides a graphical representation of the scatterer location and magnitude by plotting the points in the
%   Verasonics coordinate system.
%
%   When run without a figure handle argument, the function uses figure 50, and also calls the Media.function
%   specified in the Media structure (for 100 moves) and animates the points to show how they will be displaced
%   during the normal execution of the sequence.
%
%   See also: showGeometry, showTrans, showPData, showCoordAxes

% PK -- first version 04-JUN-2018
% YT -- Fix axis revers issue on 2-APR-2019
% PK -- Made compatible with speckle points with low scattering values (on the order of 0.02) 22-APR-2019
% PK -- Fixed bug with speckle points with low scattering vlaues; normalization was incorrect. 6-APR-2021

%% read in the Media structure
if ~evalin('base','exist(''Media'',''var'')')
    error([mfilename ': No Media structure found in base workspace.']);
end

Media = evalin('base','Media');
if ~isfield('Media', 'numPoints')
    Media.numPoints = size(Media.MP, 1);
end
Trans = evalin('base','Trans');

% Define scale, use the same definition as VSX
scale = 1.0;  % default is no scaling (wavelengths)
axisname = 'wavelengths';
if evalin('base','exist(''Resource'',''var'')')
    Resource = evalin ('base', 'Resource');
    if isfield(Resource.DisplayWindow,'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmpi(Resource.DisplayWindow(1).AxesUnits,'mm')
            axisname = 'mm';
            if ~isfield(Resource.Parameters,'speedOfSound')||isempty(Resource.Parameters.speedOfSound)
                scale = 1.54/Trans.frequency; % use default speed of sound.
            else
                scale = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
            end
        end
    end
else
    warning([mfilename ': No Resource structure found in base workspace. Set the speed of sound to 1540 m/s']);
    scale = 1.54/Trans.frequency;
end

%% Plot will go into figure handle argument. If none provided, use #50.
if nargin==0
    VsClose
    figHandle = figure(50);
    set(gcf,'Position', [500 100 800 800]);
else
    figure(figHandle);
end

% Draw in the points, with circle size based on strength
dnorm = max( max(Media.MP(:,4)), 1);        % normalize to 1, unless there are some very bright reflectors
handleMP = gobjects(length(Media.MP),1);    % preallocate one graphics object for each reflector
for nP = 1:Media.numPoints
    if nP==2, hold on, end
    d =  ( (log10( abs(Media.MP(nP,4)/dnorm) ))+2 )/2; % compress the brightness using the log, assuming that the dimmest point is 1/100 and brightest is 1.
    if d<0, d=0; end            % do not permit negatives for points dimmer than 0.01
    % if d>1, d=1; end            % do not permit values greater than 1 (should be guaranteed by normalization)
    dcolor = (1 - d*[1 1 1]);   % invert color so that dim points appear light colored
    handleMP(nP) = plot3( scale*Media.MP(nP,1), scale*Media.MP(nP,2), scale*Media.MP(nP,3), 'ok', 'MarkerSize', 12*d, 'MarkerFaceColor', dcolor, 'MarkerEdgeColor', 'none'  );
    % plot3( scale*Media.MP(nP,1), scale*Media.MP(nP,2), scale*Media.MP(nP,3), 'ok', 'MarkerSize', 12*d, 'MarkerFaceColor', dcolor, 'MarkerEdgeColor', 'none'  );
end
hold off


%% If there is no figure handle provided, display Media points with motion animation
if nargin == 0
    axis equal
    % --- set axis limits greater than the media points span by rounding up to the next significant digit
    A = axis;
    A(1) = sign(A(1))*ceil(abs(A(1)/10))*10;
    A(2) = sign(A(2))*ceil(abs(A(2)/10))*10;
    A(3) = sign(A(3))*ceil(abs(A(3)/10))*10;
    A(4) = sign(A(4))*ceil(abs(A(4)/10))*10;
    A(5) = 0;
    A(6) = sign(A(6))*ceil(abs(A(6)/10))*10;
    axis(A)

    view(-30, 20)
    xlabel (['X (' axisname ')']), ylabel (['Y (' axisname ')']), zlabel (['Z (' axisname ')']), rotate3d on
    % flip 2 sets of axes to properly represent the Right Hand coordinate system with Z downward
    set (gca, 'ZDir', 'reverse')
    set (gca, 'Ydir', 'reverse')

    if isvalid(figHandle)
        showCoordAxes(figHandle); % provide a fixed reference for the points
    else
        return  % exit if figure has been closed
    end

    if isfield(Media,'function')
        numberOfMoves = 100;
        maxPts = 300;
        figure(figHandle)
        % titleName = title ({Trans.name, [Media.function ' iteration # 1 of ' num2str(numberOfMoves)]});
        title ({Trans.name, [Media.function ' iteration # 1 of ' num2str(numberOfMoves)]});
        pause(3)

        for num = 1:numberOfMoves
            evalin ('base', Media.function)
            Media = evalin ('base', 'Media');
            if num==1 && length(Media.MP)>maxPts % remove plot of most random points by setting them whilte and overlaying on one location
                for nP = maxPts+1:length(Media.MP)
                    if isvalid(handleMP(nP))
                        set(handleMP(nP),'XData',scale*Media.MP(maxPts+1,1),'YData',scale*Media.MP(maxPts+1,2),'ZData',scale*Media.MP(maxPts+1,3), 'MarkerSize', 1, 'MarkerFaceColor', [1 1 1]);
                    else
                        return  % exit if figure has been closed
                    end
                end
            end
            for nP = 1:min(length(Media.MP), maxPts) % move the first 100 points
                if isvalid(handleMP(nP))
                    set(handleMP(nP),'XData',scale*Media.MP(nP,1),'YData',scale*Media.MP(nP,2),'ZData',scale*Media.MP(nP,3));
                else
                    return  % exit if figure has been closed
                end
            end
            % titleName.String = ({Trans.name,[Media.function ' iteration # ' num2str(num) ' of ' num2str(numberOfMoves)]});
            title ({Trans.name,[Media.function ' iteration # ' num2str(num) ' of ' num2str(numberOfMoves)]})
            pause(.05)
        end
    end
end

end
