function showTrans(figHandle)
% showTrans (figHandle)  plots the transducer elements into an existing figure or into figure 50.
%
%   The function provides a graphical representation of the element sizes and locations by drawing rectangular patches
%   in the Verasonics coordinate system. The array is placed at Z=0 but draws a rectangle at an effective depth
%   resulting from the lens correction.
%   Units in the plot are defined by the Resource.Display.AxesUnits. Units specifying the transducer coordinates in the Trans
%   structure depend on the Trans.units argument. To simplify the conversion below, the factor scale is used to convert
%   wavelengths to mm. If units are wavelengths, scale=1, and if display units are mm, scale is 1 wavelength in mm.
%
%   See also: showGeometry, showPData, showMediaPts, showCoordAxes, showTransIR, showTransImpedance

% Future additions:
%   - add elevation focal beam outline (use pairs of hyperbolic lines, or surf or surf2patch)
%
% PK -- beta version 11-JUN-2018
% YT -- refined on 29-MAR-2019 - Fully supports types 0-2

%% read in Trans structure
if ~evalin('base','exist(''Trans'',''var'')')
    error([mfilename ': No Trans structure found in base workspace.']);
end

Trans   = evalin('base', 'Trans');
ElementPos  = Trans.ElementPos;
elementWidth2   = Trans.elementWidth/2;
if ~isfield(Trans,'lensCorrection')
    lensCorrection = 0;
else
    lensCorrection  = Trans.lensCorrection;
end
if ~isfield(Trans,'name')
    Trans.name = 'custom transducer';
end
if ~isfield(Trans,'elementLength')
    elementLength2 = 0;
else
    elementLength2 = Trans.elementLength/2;
end

% Define scale, use the same definition as VSX
scale = 1.0;  % default is no scaling (wavelengths)
axisname = 'wavelengths';
if evalin('base','exist(''Resource'',''var'')')
    Resource = evalin ('base', 'Resource');
    
    verbose = evalin('base', 'Resource.Parameters.verbose');
    if verbose>2
        disp([mfilename ' Version 3-Apr-2018'])
    end
    
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

if isfield (Trans, 'elevationApertureMm') && ~isempty(Trans.elevationApertureMm) && ~isequal(Trans.elevationApertureMm,0)
    elev = Trans.elevationApertureMm/scale; %in wavelength
else
    if strcmp(Trans.units, 'wavelengths')
        elev = 5 * Trans.elementWidth;
    else
        elev = 5 * Trans.elementWidth / scale;
    end
end

%% Plot will go into figure defined by the handle argument. If none provided, use #50. Assume only 1 Display window defined in script.
if nargin==0
    VsClose
    figHandle = figure(50);
    set(gcf,'Position', [400 70 1200 950]);
else
    figure(figHandle);
end

if strcmp(Trans.units,'mm') && strcmp(axisname, 'wavelengths')
    ElementPos(:,1:3) = ElementPos(:,1:3)/scale;
    lensCorrection = lensCorrection/scale;
    elementWidth2 = elementWidth2/scale;
    elementLength2 = elementLength2/scale;
elseif strcmp(Trans.units,'wavelengths') && strcmp(axisname, 'mm')
    ElementPos(:,1:3) = ElementPos(:,1:3)*scale;
    lensCorrection = lensCorrection*scale;
    elementWidth2 = elementWidth2*scale;
    elementLength2 = elementLength2*scale;
    elev = elev*scale;
elseif strcmp(Trans.units,'mm') && strcmp(axisname, 'mm')
    elev = elev*scale;
end

%% plot all transducer elements in 3D  --- USE PATCH to make the elements oriented in the right directions
switch Trans.type
    case 0      % ---------- Linear  --------------------------------------
        for n=1:Trans.numelements   % plot each element as a red patch the size of the active area at depth given by Lens Correction
            Xc = ElementPos(n,1);
            X = Xc+ elementWidth2*[-1 -1 1 1];
            Y = elev*[-1 1 1 -1]/2;
            Z = [0 0 0 0];
            
            if n==1
                patch (X,Y,Z, [.4 .4 .4], 'edgecolor', 'none'),  view(-215,20); % element #1 is gray for easy identification
                hold on
                text (Xc(n), Y(1), Z(1), 'el 1','FontSize',14)
            else
                patch (X,Y,Z, [.8 .4 .4], 'edgecolor', 'none'),  view(-215,20);
            end
        end
        
        % lensCorrection
        plot3([ElementPos(1,1),ElementPos(1,1),ElementPos(end,1),ElementPos(end,1),ElementPos(1,1),],...
            [-1 1 1 -1 -1]*elev/2,[1 1 1 1 1]*lensCorrection,'r','LineWidth',2);
        text (ElementPos(end,1), elev/2, lensCorrection, 'lensCorrection', 'FontSize', 14)
        
    case 1      % ---------- Curvilinear ------------------------------------
        
        for n=1:Trans.numelements
            Xc = ElementPos(n,1);
            Zc = ElementPos(n,3);
            ang = ElementPos(n,4); %-(theta + (Trans.spacing*(n-0.5)/radius));
            
            % Rotate element
            El = [Xc-elementWidth2*cos(ang),Zc+elementWidth2*sin(ang)];
            Er = [Xc+elementWidth2*cos(ang),Zc-elementWidth2*sin(ang)];
            
            X = [El(1) El(1) Er(1) Er(1)];
            Y = elev*[ -1 1 1 -1 ]/2;
            Z = [El(2) El(2) Er(2) Er(2)];
            
            if n==1
                patch (X,Y,Z, [.4 .4 .4], 'edgecolor', 'none'),  view(-215,20); % element #1 is gray for easy identification
                hold on
                text (Xc(n), Y(1), Z(1), 'el 1')
            else
                patch (X,Y,Z, [.8 .4 .4], 'edgecolor', 'none'),  view(-215,20);
            end
        end
        
        numEle = Trans.numelements;
        
        % lensCorrection
        plot3([ElementPos(:,1);flip(ElementPos(:,1));ElementPos(1,1)],...
            [-1*ones(1,numEle),ones(1,numEle),-1]*elev/2,([ElementPos(:,3);flip(ElementPos(:,3));ElementPos(1,3)]+lensCorrection),'r','LineWidth',2);
        text (ElementPos(end,1), elev/2, (ElementPos(end,3)+lensCorrection+2), 'lensCorrection', 'FontSize', 14)
        
    case 2     % ---------- 2D array (elements on a spherical or planar surface) -------------------------------
        
        for n=1:Trans.numelements   % plot each element as a red patch the size of the active area at depth given by Lens Correction
            Xc = ElementPos(n,1);
            X = Xc + elementWidth2*[-1 -1 1 1];
            Yc = ElementPos(n,2);
            Y = Yc + elementWidth2*[-1 1 1 -1];
            Zc = ElementPos(n,3);
            Z = [0 0 0 0];
            
            if ismember(Trans.name,{'H-301','H-302','H-313'}) % FUS Elite 3000 3D
                az = ElementPos(n,4);
                el = ElementPos(n,5);
                % The unit of Trans.radius is always wavelengths
                R = Trans.radius;
                if strcmp(axisname, 'mm')
                    R = Trans.radius*scale;
                end
                
                r = elementWidth2;    % Radius of circle
                R1 = sqrt(R^2+r^2);
                theta=0:0.1:2*pi;
                
                % the az and el of each point on the circle needs to be
                % adjusted
                az1 = atan(r*cos(theta)/R);
                el1 = asin(r*sin(theta)/R1);
                
                az2 = az+az1;
                el2 = el+el1;
                
                Y = R1*sin(el2);
                X = R1*cos(el2).*sin(az2);
                Z = R1-R1*cos(el2).*cos(az2);
            end
            
            patch (X,Y,Z, [.8 .4 .4], 'edgecolor', 'none'),  view(-215,20);
            text(Xc,Yc,Zc, num2str(n), 'fontsize', 10, 'fontweight', 'normal')
            
        end
        
    case 4     % ---------- Row-Column array (elements on a planar surface) -------------------------------
        
        % plot row
        rInd = find(ElementPos(:,1));
        for n = rInd(1):rInd(end)
            Xc = ElementPos(n,1);
            X = Xc + elementWidth2*[-1 -1 1 1];
            Yc = ElementPos(n,2);
            Y = Yc + elementLength2*[-1 1 1 -1];
            Zc = ElementPos(n,3);
            Z = [0 0 0 0];
            patch (X,Y,Z, 'facecolor', 'r', 'edgecolor','none','FaceAlpha',0.2),  view(-215,20);
            text(Xc,Yc,Zc, num2str(n), 'fontsize', 8, 'fontweight', 'normal')
        end
        
        % plot column
        cInd = find(ElementPos(:,2));
        for n = cInd(1):cInd(end)
            Xc = ElementPos(n,1);
            X = Xc + elementLength2*[-1 -1 1 1];
            Yc = ElementPos(n,2);
            Y = Yc + elementWidth2*[-1 1 1 -1];
            Zc = ElementPos(n,3);
            Z = [0 0 0 0];
            patch (X,Y,Z, 'facecolor','b', 'edgecolor', 'none','FaceAlpha',0.2),  view(-215,20);
            text(Xc,Yc,Zc, num2str(n), 'fontsize', 8, 'fontweight', 'normal')
        end
        
    otherwise
        disp ([mfilename ': Trans.type = ' num2str(Trans.type) ' not yet supported'])
        return
end

hold off

%% Finish the plot appearance if running function as standalone
if nargin==0
    axis equal
    
    % extend axes to a view larger than the largest element
    A = axis;
    dX = (A(2)-A(1))/10;    dY = (A(4)-A(3))/10;    dZ = (A(6)-A(5))/10;
    A  = [A(1)-dX, A(2)+dX,  A(3)-dY, A(4)+dY,  A(5)-dZ, A(6)+dZ ];
    axis(A);
    
    % label and reverse axes
    xlabel (['X (' axisname ')']), ylabel (['Y (' axisname ')']), zlabel (['Z (' axisname ')']), rotate3d on
    title (['Transducer geometry for ' Trans.name], 'interpreter', 'none')
    set (gca, 'ZDir', 'reverse')
    set (gca, 'Ydir', 'reverse')
    
    showCoordAxes (figHandle);
    % view(180, -90) % view from +Z side onto the array, i.e., from the acosutic medium
    
    h2 = figure (12); h2.Position= [20 20 800 600];
    showTransIR (h2);
    
    %   h3 = subplot (2,2,3);
    h3 = figure (13); h3.Position= [850 20 800 600];
    showTransImpedance (13);
end



end
