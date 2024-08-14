function  hROI = drawRegionOutline(varargin)
% drawRegionOutline creates an outline of the region defined in
% PData(PDataNum).Region(RegionNum). If ROI handle has been created, it can also
% be used to set outline attributes, like visibility, and new location
%
% drawRegionOutline(hROI,'off') makes the ROI invisible ('on' is visible)
%
% drawRegionOutline(WinNum,PDataNum,RegionNum) creates an outline from
% PData(PDataNum).Region(RegionNum) on displayWindow(WinNum) with default color - white.
%
% drawRegionOutline(WinNum,PDataNum,RegionNum,LineColor) creates an outline with a custom color
%
% drawRegionOutline(hROI,WinNum,PDataNum,RegionNum) replot the outline based on the
% PData(PDataNum).Region(RegionNum)
%
% drawRegionOutline(hROI,WinNum,PDataNum,RegionNum,LineColor) replot the outline based on the
% PData(PDataNum).Region(RegionNum) with different color

% import java path for Verasonics Viewer
import java.awt.Color;
import com.verasonics.viewer.ui.overlays.*
import com.verasonics.viewer.ui.panels.*

Trans = evalin('base','Trans');
PData = evalin('base','PData');
Resource = evalin('base','Resource');

[hROI,WinNum,PDataNum,RegionNum,LineColor] = errorCheck(varargin);
% pass the errorCheck
if isequal(length(varargin),2)
    status = lower(varargin{2});
    switch status
        case 'on'
            if ishghandle(hROI)
                hROI.Visible = status;
            else
                hROI.isEnabled = true;
            end
        case 'off'
            if ishghandle(hROI)
                hROI.Visible = status;
            else
                hROI.isEnabled = false;
            end
    end
    return % turn on/off the outline only
end

% Draw a new outline or set the position
Region = PData(PDataNum).Region(RegionNum);
VsType = Resource.DisplayWindow(WinNum).Type;

if isfield(Region,'Shape')&&(~isempty(Region.Shape))
    Shape = Region.Shape;
else
    fprintf(2,'showRegionOutline: no Shape attribute.\n');
    return
end

switch Shape.Name
    case {'Rectangle','Parallelogram'}

        if ~isfield(Shape,'angle')||isempty(Shape.angle)
            Shape.angle = 0;  % default is rectangle shape
        end

        X(1) = Shape.Position(1)-Shape.width/2;
        Y(1) = Shape.Position(3);
        X(2) = X(1) + Shape.height*tan(Shape.angle);
        X(3) = X(2) + Shape.width;
        Y(2) = Y(1) + Shape.height;
        Y(3) = Y(2);
        X(4) = X(1) + Shape.width;
        Y(4) = Y(1);

        if Shape.angle > 0
            rightEnd = PData(PDataNum).Origin(1)+PData(PDataNum).Size(2)*PData(PDataNum).PDelta(1);
            if X(2) > rightEnd
               X(2) = rightEnd;
               X(3) = X(2);
               Y(2) = (X(2)-X(1))/tan(Shape.angle);
               Y(3) = (rightEnd - X(4))/tan(Shape.angle);
            elseif X(3) > rightEnd
                X(5) = X(4); Y(5) = Y(4);
                X(3) = rightEnd;
                Y(4) = (rightEnd - X(4))/tan(Shape.angle);
                X(4) = X(3);
            end
        else
            if X(2) < PData(PDataNum).Origin(1)
                Y(2) = Y(1) + (PData(PDataNum).Origin(1)-X(1))/tan(Shape.angle);
                if (X(4)-Shape.height*abs(tan(Shape.angle)))< PData(PDataNum).Origin(1)
                    X(2) = PData(PDataNum).Origin(1);
                    X(3) = X(2);
                    Y(3) = (X(4) - PData(PDataNum).Origin(1))/abs(tan(Shape.angle));
                else
                    X(5) = X(4); Y(5) = Y(4);
                    X(4) = X(3); Y(4) = Y(3);
                    X(3) = PData(PDataNum).Origin(1);
                    X(2) = X(3);
                end
            end
        end

        ROIX = [X X(1)];
        ROIY = [Y Y(1)];

        % scale to mm if the axes unit of the displaywindow is mm
        scaleToMm = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
        if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(WinNum).AxesUnits)
            if strcmp(Resource.DisplayWindow(WinNum).AxesUnits,'mm')
                ROIX = ROIX * scaleToMm; ROIY = ROIY * scaleToMm;
            end
        end

        if strcmp(VsType,'Matlab')
            bmodeFigHandle = Resource.DisplayWindow(WinNum).figureHandle;
            if ishandle(bmodeFigHandle)
                if isempty(hROI) || ~isgraphics(hROI)
                    bmodeFigHandle.CurrentAxes.NextPlot = 'add';
                    hROI = plot(bmodeFigHandle.CurrentAxes,ROIX,ROIY,'HandleVisibility','off',...
                        'Color',LineColor,'LineWidth',1.5,'LineStyle','--');
                    %@Daniel I need access to this to turn it off and on
                    setappdata( bmodeFigHandle.CurrentAxes, 'VsRegionOutline', hROI ) ;
                else
                    uistack(hROI,'top')
                    set(hROI,'XData',ROIX,'YData',ROIY,'Color',LineColor);
                end
            end
        else
            % ImageViewerPanel is the class in com.verasonics.viewer.ui.panels.*
            % getViewerFromId is the method of the ImageViewerPanel class
            uid = Resource.DisplayWindow(WinNum).imageHandle;
            imageViewer = ImageViewerPanel.getViewerFromId(uid);
            imageCanvas = imageViewer.getImageCanvas();
            srcImage = imageCanvas.getSrcImage();

            if isempty(hROI)
                hROI = ColorBoxOverlayRenderer(imageCanvas);
                hROI.isEnabled = true;
                hROI.shape = Shape.Name;
                imageViewer.addOverlayRenderer(hROI);
            end

            hROI.lineColor = LineColor;
            hROI.linearPosX = (ROIX - imageCanvas.getXPixelOffset())/srcImage.xPixelSize;
            hROI.linearPosY = (ROIY - imageCanvas.getYPixelOffset())/srcImage.yPixelSize;

        end

    case {'Sector','Circle','Annulus','SectorFT'}
        if strcmp(Shape.Name,'Annulus')&&(~isfield(Shape,'angle')||isempty(Shape.angle))
            Shape.angle = 2*pi;
        end
        if strcmp(Shape.Name,'Circle')
            if (~isfield(Shape,'r')||isempty(Shape.r))
                error('computeRegions: Missing r attribute for ''Circle'' shape in PData.Region(%d).',nr);
            else
                Shape.r1 = 0;
                Shape.r2 = Shape.r;
                Shape.angle = 2*pi;
            end
        end
        if ~isfield(Shape,'steer')||isempty(Shape.steer)
            Shape.steer = 0;
        end

        if strcmp(Shape.Name,'SectorFT')
            sectorFT = 1;
            Shape.r1 = Shape.z-Shape.Position(3);
            Shape.r2 = Shape.r;
        else
            sectorFT = 0;
        end

        r1 = Shape.r1;
        r2 = Shape.r2;
        angle = Shape.angle; % in radian
        steerAngle = Shape.steer; % in radian
        theta = pi/2-angle/2-steerAngle; % start angle
        circlePts = 180; % points for whole circle
        n = abs(ceil(circlePts*angle/(2*pi)));
        if ~iseven(n), n = n+1; end

        curve_theta = theta + (angle*(0:n)'/n);
        curve1_R = r1*ones(n+1,1);
        curve2_R = r2*ones(n+1,1);

        [curve1_X,curve1_Y] = pol2cart(curve_theta,curve1_R);
        [curve2_X,curve2_Y] = pol2cart(curve_theta,curve2_R);

        % correct ROI points if curve2 is outside of the boundry
        xmin = Resource.DisplayWindow(WinNum).ReferencePt(1);
        xmax = xmin + Resource.DisplayWindow(WinNum).Position(3)*Resource.DisplayWindow(WinNum).pdelta;
        if min(curve2_X) < xmin
            X1 = curve1_X(end); Y1 = curve1_Y(end);
            X2 = curve2_X(end); Y2 = curve2_Y(end);
            X = xmin; Y = Y1 + (X-X1)*(Y2-Y1)/(X2-X1);
            curve2_X(curve2_X<xmin) = xmin;
            curve2_X(end+1) = X;
            curve2_Y(end+1) = Y;
        elseif max(curve2_X) > xmax
            X1 = curve1_X(1); Y1 = curve1_Y(1);
            X2 = curve2_X(1); Y2 = curve2_Y(1);
            X = xmax; Y = Y1 + (X-X1)*(Y2-Y1)/(X2-X1);
            curve2_X(curve2_X>xmax) = xmax;
            curve2_X = [X; curve2_X];
            curve2_Y = [Y; curve2_Y];
        end

        curve1_X = curve1_X + Shape.Position(1); % shift X

        if sectorFT
            curve1_Y = r1*ones(size(curve1_X)) + Shape.Position(3);
        else
            curve1_Y = curve1_Y + Shape.Position(3); % shift Y
        end

        curve2_X = curve2_X + Shape.Position(1);
        curve2_Y = curve2_Y + Shape.Position(3);

        % scale to mm if the axes unit of the displaywindow is mm
        scaleToMm = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
        if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
            if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
                curve1_X = curve1_X * scaleToMm; curve2_X = curve2_X * scaleToMm;
                curve1_Y = curve1_Y * scaleToMm; curve2_Y = curve2_Y * scaleToMm;
            end
        end

        ROIX = [curve1_X; flipud(curve2_X); curve1_X(1)];
        ROIY = [curve1_Y; flipud(curve2_Y); curve1_Y(1)];

        if strcmp(VsType,'Matlab')
            bmodeFigHandle = Resource.DisplayWindow(WinNum).figureHandle;
            if ishandle(bmodeFigHandle)
                if  isempty(hROI) || ~isgraphics(hROI)
                    bmodeFigHandle.CurrentAxes.NextPlot = 'add';
                    hROI = plot(bmodeFigHandle.CurrentAxes,ROIX,ROIY,'HandleVisibility','off',...
                        'Color',LineColor,'LineWidth',1.5,'LineStyle','--');
                    %@Daniel I need access to this to turn it off and on
                    setappdata( bmodeFigHandle.CurrentAxes, 'VsRegionOutline', hROI ) ;
                else
                    uistack(hROI,'top')
                    set(hROI,'XData',ROIX,'YData',ROIY,'Color',LineColor);
                end
            end
        else
            % ImageViewerPanel is the class in com.verasonics.viewer.ui.panels.*
            % getViewerFromId is the method of the ImageViewerPanel class
            uid = Resource.DisplayWindow(WinNum).imageHandle;
            imageViewer = ImageViewerPanel.getViewerFromId(uid);
            imageCanvas = imageViewer.getImageCanvas();
            srcImage = imageCanvas.getSrcImage();

            if isempty(hROI)
                hROI = ColorBoxOverlayRenderer(imageCanvas);
                hROI.isEnabled = true;
                hROI.shape = Shape.Name;
                imageViewer.addOverlayRenderer(hROI);
            end

            % control point of the quadCurve2D is calculated by
            % B(t) = (1-t)^2*P0+2*t*(1-t)*P1+t^2*P2, 0<t<1
            % When t = 0.5, B(0.5) is the middle of curve
            % so ctrl point can be derived be the equation below
            % refer to https://en.wikipedia.org/wiki/Bezier_curve
            xOff = imageCanvas.getXPixelOffset(); xPS = srcImage.xPixelSize;
            yOff = imageCanvas.getYPixelOffset(); yPS = srcImage.yPixelSize;

            c1start = [(curve1_X(1)-xOff)/xPS,(curve1_Y(1)-yOff)/yPS];
            c1end = [(curve1_X(end)-xOff)/xPS,(curve1_Y(end)-yOff)/yPS];
            c1middle = [(curve1_X(n/2+1)-xOff)/xPS,(curve1_Y(n/2+1)-yOff)/yPS];
            ctrl1 = 2*c1middle-c1start/2-c1end/2;

            c2start = [(curve2_X(1)-xOff)/xPS,(curve2_Y(1)-yOff)/yPS];
            c2end = [(curve2_X(end)-xOff)/xPS,(curve2_Y(end)-yOff)/yPS];
            c2middle = [(curve2_X(n/2+1)-xOff)/xPS,(curve2_Y(n/2+1)-yOff)/yPS];
            ctrl2 = 2*c2middle-c2start/2-c2end/2;

            hROI.curve1Pos = floor([c1start,ctrl1,c1end]);
            hROI.curve2Pos = floor([c2start,ctrl2,c2end]);
            hROI.lineColor = LineColor;

        end

    otherwise
        fprintf(2,'showRegionOutline does not support %s shape.\n',Shape.Name);
end

    function [hROI,WinNum,PDataNum,RegionNum,LineColor] = errorCheck(args)
        import java.awt.Color;
        hROI = [];
        WinNum = [];
        PDataNum = [];
        RegionNum = [];
        LineColor = [];

        % The first argument is either region handle or WinNum
        if (isgraphics(args{1}) && strcmp(get(args{1},'Type'),'line'))||isjava(args{1})
            hROI = args{1};
            % If the first argument is region handle, the second one is
            % either on/off or WinNum
            if ~strcmpi(args{2},'on') && ~strcmpi(args{2},'off')
                WinNum = args{2};
            end
        else
            WinNum = args{1};
        end

        if ~isempty(WinNum)
            if WinNum > length(Resource.DisplayWindow) || ~isnumeric(WinNum)
                error('showRegionOutline: incorrect displayWindow number.');
            end

            switch length(args)
                case 3 % drawRegionOutline(WinNum,PDataNum,RegionNum)
                    PDataNum = args{2};
                    if PDataNum > length(PData) || ~isnumeric(PDataNum)
                        error('showRegionOutline: incorrect PData number.');
                    end
                    RegionNum = args{3};
                    if RegionNum > length(PData(PDataNum).Region) || ~isnumeric(RegionNum)
                        error('showRegionOutline: incorrect Region number.');
                    end
                    % set default line color
                    if isempty(LineColor) || isnumeric(LineColor)
                        if ~strcmp(Resource.DisplayWindow(WinNum).Type,'Matlab')
                            LineColor = Color.white;
                        else
                            LineColor = 'w';
                        end
                    end
                case 4 % drawRegionOutline(WinNum,PDataNum,RegionNum,LineColor) or drawRegionOutline(ROI,WinNum,PDataNum,RegionNum)
                    LineColor = args{4};
                    if ischar(LineColor)
                        if ~strcmp(Resource.DisplayWindow(WinNum).Type,'Matlab')
                            if ~ismember(LineColor,{'b','g','r','c','m','y','k','w'})
                                error('showRegionOutline: not supported color.');
                            else
                                switch LineColor
                                    case 'b'
                                        LineColor = Color.blue;
                                    case 'g'
                                        LineColor = Color.green;
                                    case 'r'
                                        LineColor = Color.red;
                                    case 'c'
                                        LineColor = Color.cyan;
                                    case 'm'
                                        LineColor = Color.magenta;
                                    case 'y'
                                        LineColor = Color.yellow;
                                    case 'k'
                                        LineColor = Color.black;
                                    case 'w'
                                        LineColor = Color.white;
                                end
                            end
                        end
                        [~,~,PDataNum,RegionNum,~] = errorCheck(args(1:3));
                    else % Use the same color
                        if isjava(hROI)
                            LineColor = hROI.lineColor;
                        else
                            LineColor = hROI.Color;
                        end
                        [~,~,PDataNum,RegionNum,~] = errorCheck(args(2:4));
                    end
                case 5 %drawRegionOutline(ROI,WinNum,PDataNum,RegionNum,LineColor)
                    [~,~,PDataNum,RegionNum,LineColor] = errorCheck(args(2:5));
                otherwise
                    error('showRegionOutline: incorrect input arguments.');
            end
        end
    end

end

function ise = iseven( val )
    ise = vsv.math.EquMath.iseven(val);
end
