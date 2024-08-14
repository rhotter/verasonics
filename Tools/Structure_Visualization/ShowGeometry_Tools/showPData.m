function showPData (figHandle)
% showPData (figHandle) plots the PData structures as shaded patches in three dimensions into an existing figure
%
%   The function provides a graphical representation of the pixel grid defined by any number of Pdata structures in the
%   Verasonics coordinate system.
%   The function also permits viewing the Regions defined for each PData, using a slider control.
%
%   See also: showGeometry, showTrans, showMediaPts, showCoordAxes

% PK -- beta version 11-JUN-2018
% YT -- 29-MAR-2019 - support multiple PData structures, regions and polar coordinate

persistent numCols numRows numSecs xOrg yOrg zOrg ...
    Xind Yind Zind Rind Tind nregtxt nregSldr nregValue ...
    numRegions PDataPlot RegionPlot

%% read in PData structures
if ~evalin('base','exist(''PData'',''var'')')
    error([mfilename ': No PData structure found in base workspace.']);
end

PData = evalin ('base', 'PData');
Trans = evalin('base', 'Trans');

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

polarCoord = 0;

% Plot will go into figure handle argument. If none provided, use #50.
if nargin==0
    VsClose
    figHandle = figure(50);
    set(gcf,'Position', [500 100 800 800]);
else % called by showGeometry
    figure(figHandle);
end

% Additional UI Controls
bkgrnd = [0.8 0.8 0.8];
SG = struct('TO',[0.0,0.07],...    % title offset
    'TS',[0.15,0.03],...    % title size
    'TF',0.5,...            % title font size
    'SO',[0.0,0.05],...     % slider offset
    'SS',[0.15,0.025],...   % slider size
    'EO',[0.035,0.025],...  % edit box offset
    'ES',[0.08,0.025]);     % edit box size

% Multiple PData structures - DropDown menu
sizeP = size(PData,2);
if sizeP > 1
    Pos = [0.005 0.88];
    PString = cell(sizeP,1);
    for j = 1:sizeP
        PString{j} = ['   ',num2str(j,1)];
    end

    uicontrol(figHandle,'Style','text',...
        'String','PData No.',...
        'Units','normalized',...
        'Position',[Pos+SG.TO,SG.TS],...
        'FontUnits','normalized',...
        'FontSize',SG.TF,...
        'FontWeight','bold');
    nPDataMenu = uicontrol('Style','popupmenu',...
        'String',PString,...
        'Units','normalized',...
        'Position',[Pos+SG.SO+[0.025,0],0.1,0.025],...
        'FontUnits','normalized',...
        'FontSize',0.7,...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'tag','viewButton',...
        'Callback',{@numPData_Callback});
end

% Define variables for nested functions
nPD = 1;

% Plot
numPData_Callback;

%% nested function - PData Callback
    function numPData_Callback(varargin)

        if ishandle(PDataPlot), delete(PDataPlot); end
        if ishandle(RegionPlot), delete(RegionPlot); end

        if size(PData,2) > 1
            nPD = nPDataMenu.Value;
        end

        % use the same definition in computeRegions
        numRows = PData(nPD).Size(1);
        numCols = PData(nPD).Size(2);
        if size(PData(nPD).Size,2)==3
            numSecs = PData(nPD).Size(3);
        else
            numSecs = 1;
        end
        xOrg = PData(nPD).Origin(1);
        yOrg = PData(nPD).Origin(2);
        zOrg = PData(nPD).Origin(3);

        % simplified version. Only support polar and rectangular (default)
        % on 4/1/2019

        pdeltaX = PData(nPD).PDelta(1);
        pdeltaY = PData(nPD).PDelta(2);
        pdeltaZ = PData(nPD).PDelta(3);

        if isfield (PData(nPD),'Coord')
            if strcmp(PData(nPD).Coord, 'polar')
                polarCoord = 1;
                pdeltaT = PData(nPD).PDelta(1);
                pdeltaR = PData(nPD).PDelta(2);
            end
        end

        % axis and patch parameters
        C = [ .8 .9 .8 ];
        Alpha = .6;
        Axlims0 = zeros(1,6);

        figure(figHandle), hold on

        if min (PData(nPD).Size)==1 % 2D planes only
            if numSecs==1 % X,Z plane
                % range, min and max coordinates (wvlns)
                if polarCoord % Currently, polar coordinate only supports 2D XZ plane
                    % need to define multiple points for curved PData
                    Tind = -pdeltaT*(numCols-1)/2:pdeltaT:pdeltaT*(numCols-1)/2;
                    Rind = (0:pdeltaR:((numRows-1)*pdeltaR))'*scale;

                    % coordinates for the patch
                    r = pdeltaR*numRows*scale;
                    X = zeros(1,length(Tind)+1);
                    Y = zeros(1,length(Tind)+1);
                    Z = zeros(1,length(Tind)+1);

                    for t = 1:length(Tind)
                        [Z(t),X(t)] = pol2cart(Tind(t),r);
                    end

                    X(t+1) = 0; X = X + xOrg*scale; % need to shift by Origin
                    Z(t+1) = 0; Z = Z + zOrg*scale;

                    % extend the plot region beyond PData itself
                    X1 = min(X); X2 = max(X); dX = X2 - X1;
                    Z1 = min(Z); Z2 = max(Z); dZ = Z2 - Z1;

                    % Create index map for linear address of Region mapping later.
                    Xind = repmat(Tind,numRows,1);
                    Zind = repmat(Rind,1,numCols);
                    Yind = zeros(PData(nPD).Size(:,:,1));

                else % Cartesian coord
                    Xind = (xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX))*scale;
                    Zind = (zOrg:pdeltaZ:(zOrg+(numRows-1)*pdeltaZ))'*scale;

                    % coordinates for the 4 vertices in the patch
                    X1 = min(Xind); X2 = max(Xind); dX = X2 - X1;
                    Z1 = min(Zind); Z2 = max(Zind); dZ = Z2 - Z1;

                    X = [X1 X2 X2 X1];
                    Z = [Z1 Z1 Z2 Z2];
                    Y = zeros(1,4);

                    % Create index map for linear address of Region mapping later.
                    Xind = repmat(Xind,numRows,1);
                    Zind = repmat(Zind,1,numCols);
                    Yind = zeros(PData(nPD).Size(:,:,1));
                end

                Axlims = [X1-dX/10 X2+dX/10 -dX/10 +dX/10 Z1 Z2+dZ/10];

            else % XY plane
                Xind = (xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX))*scale;
                Yind = (yOrg:pdeltaY:(yOrg+(numRows-1)*pdeltaY))'*scale;

                % coordinates for the 4 vertices in the patch
                X1 = Xind(1); X2 = Xind(end); dX = diff([X1 X2]);
                Y1 = Yind(1); Y2 = Yind(end); dY = diff([Y1 Y2]);

                X = [X1 X2 X2 X1];
                Y = [Y1 Y1 Y2 Y2];
                Z = zeros(1,4);

                Axlims = [X1-dX/10 X2+dX/10 Y1-dY/10 Y2+dY/10 zOrg-dX/10 zOrg+dX/10];

                % Create index map for linear address of Region mapping later.
                Xind = repmat(Xind,numRows,1);
                Yind = repmat(Yind,1,numCols);
                Zind = zeros(PData(nPD).Size(:,:,1));
            end

            PDataPlot = patch(X,Y,Z, C, 'FaceAlpha', Alpha, 'linewidth', 1);

        else
            % The Origin of a 3-D PData space is located at [-X, +Y, 0]
            Xind = (xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX))*scale;
            Yind = (yOrg:-pdeltaY:(yOrg-(numRows-1)*pdeltaY))'*scale;
            Zind = (zOrg:pdeltaZ:(zOrg+(numSecs-1)*pdeltaZ))'*scale;

            X1 = Xind(1); X2 = Xind(end); dX = diff([X1 X2]);
            Y1 = -Yind(1); Y2 = -Yind(end); dY = diff([Y1 Y2]);
            Z1 = Zind(1); Z2 = Zind(end); dZ = diff([Z1 Z2]);

            % Below shows how to patch 6 surfaces
            X = [X1 X1 X2 X2 X1 X1 X2 X2];
            Y = [Y1 Y2 Y2 Y1 Y1 Y2 Y2 Y1];
            Z = [Z1 Z1 Z1 Z1 Z2 Z2 Z2 Z2];
            my_vertices = [X ;Y; Z]';
            my_faces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            PDataPlot = patch('Vertices', my_vertices, 'Faces', my_faces, 'FaceColor',C, 'FaceAlpha', Alpha,'edgecolor', 'k');

            % extend the plot region beyond PData itself
            Axlims = [X1-dX/10 X2+dX/10 Y1-dY/10 Y2+dY/10 Z1-dZ/10 Z2+dZ/10];

            % Create index map for linear address of Region mapping later.
            Xind = repmat(Xind,[numCols,1,numSecs]);
            Yind = repmat(Yind,[1,numRows,numSecs]);
            Zind = permute(repmat(Zind,[1,numRows,numCols]),[3 2 1]);
        end

        axis equal;
        A = [ min(Axlims(1), Axlims0(1)) max(Axlims(2), Axlims0(2)) min(Axlims(3), Axlims0(3)) max(Axlims(4), Axlims0(4)) min(Axlims(5), Axlims0(5)) max(Axlims(6), Axlims0(6))];
        axis ( A );

        % If Region field doesn't exist, call computeRegions to create one.
        if ~isfield(PData(nPD),'Region') || isempty(PData(nPD).Region)
            PData(nPD).Region = computeRegions(PData(nPD));
        end

        % delete previous handles for different PData
        delete([nregtxt,nregSldr,nregValue]);

        % Create a slider for Region display
        Pos = [0.82 0.01];
        numRegions = size(PData(nPD).Region,2);
        nregtxt = uicontrol(figHandle,'Style','text',...
            'String','Region No.',...
            'Units','normalized',...
            'Position',[Pos+SG.TO,SG.TS],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');
        nregSldr = uicontrol(figHandle,'Style','slider',...
            'Max',numRegions+1,'Min',0,'Value',0,...
            'SliderStep',[ 1/(numRegions+1), 5/(numRegions+1)],...
            'Units','normalized',...
            'Position',[Pos+SG.SO,SG.SS],...
            'BackgroundColor',bkgrnd-0.05,...
            'Interruptible','off',...
            'BusyAction','cancel',...
            'Callback',{@numRegion_Callback});
        if isequal(numRegions,1), set(nregSldr,'Max',1,'SliderStep',[1,1]); end
        nregValue = uicontrol(figHandle,'Style','edit',...
            'String','None',...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd+0.1,...
            'Callback',{@numRegion_Callback});

        % If there is no figure handle provided, label axes with correct units and
        % flip 2 sets of axes to properly represent the Right Hand coordinate system with Z downward
        showCoordAxes (figHandle);
        if nargin==0
            set (gca, 'ZDir', 'reverse','Ydir', 'reverse')
            xlabel (['X (' axisname ')']), ylabel (['Y (' axisname ')']), zlabel (['Z (' axisname ')']), rotate3d on
            title ('show PData 3D', 'interpreter', 'none')
        end
        hold off

    end

%% nested function - Region Callback
    function numRegion_Callback(hObject,~)

        Cntrl = get(hObject,'Style');
        if strcmp(Cntrl,'slider')
            nR = round(get(hObject,'Value'));
        else
            Str = hObject.String;
            if isnan(str2double(Str))
                switch Str
                    case 'All'
                        nR = numRegions+1;
                    otherwise
                        nR = 0;
                end
            else
                nR = round(str2double(Str));
                nR = max(min(numRegions+1,nR),1);
            end
        end

        % Set slider value, edit box, and update axes
        if ishandle(RegionPlot)
            delete(RegionPlot);
        end

        nregSldr.Value = nR;
        switch nR
            case 0
                nregValue.String = 'None';
            case numRegions+1
                nregValue.String = 'All';
                hold(gca,'on')
                for n = 1:numRegions
                    X = Xind(PData(nPD).Region(n).PixelsLA+1);
                    Z = Zind(PData(nPD).Region(n).PixelsLA+1);

                    if min(PData(nPD).Size)==1 % 2D plane is [Z,X,Y]
                        Y = zeros(length(X),1);
                    else
                        Y = Yind(PData(nPD).Region(n).PixelsLA+1);
                    end

                    % If polar coord, [X,Y,Z] represents [T,R,Y], and needs to be converted to cart
                    if polarCoord
                        [cZ,cX,cY] = pol2cart(X,Z,Y);
                        RegionPlot(n) = plot3(cX,cY,(cZ+zOrg*scale),'.','Color',rand(1,3));
                    else
                        RegionPlot(n) = plot3(X,Y,Z,'.','Color',rand(1,3));
                    end
                end
                hold(gca,'off');
            otherwise
                nregValue.String = num2str(nR);
                X = Xind(PData(nPD).Region(nR).PixelsLA+1);
                Z = Zind(PData(nPD).Region(nR).PixelsLA+1);

                if min(PData(nPD).Size)==1 % 2D plane is [Z,X,Y]
                    Y = zeros(length(X),1);
                else
                    Y = Yind(PData(nPD).Region(nR).PixelsLA+1);
                end

                hold(gca,'on')
                % If polar coord, [Z,X,Y] represents [R,T,Y], and needs to be converted to cart
                if polarCoord
                    [cZ,cX,cY] = pol2cart(X,Z,Y);
                    RegionPlot = plot3(cX,cY,(cZ+zOrg*scale),'.','Color','b');
                else
                    RegionPlot = plot3(X,Y,Z,'.','Color','b');
                end
                hold(gca,'off');
        end

    end
end
