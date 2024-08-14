function [Region] = computeRegions(PData)
%
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% [Region] = computeRegions(PData) Calculates PData.Region structure attributes of numPixels and
% PixelsLA, using information in the Shape structure.  The pixels specified are those pixels in
% the intersection of the Shape with the PData array region.  If no PData.Region structures are
% provided, a single Region is created that specifies the entire PData region, and is given the
% shape name of 'PData'.
%
% 2D Shapes currently supported:
%   Annulus
%   Circle
%   Parallelogram
%   PData
%   Rectangle
%   Sector
%   SectorFT
%   Trapezoid
% 3D Shapes currently supported:
%   Pyramid
%   Frustum
%   Cone
%   Slice

% Last update
% 01/02/2020 Added Frustum shape to 3D shapes
% 12/14/2018 VTS-1004 correct Coord spelling in error messages
% 12/02/2015 Added Shape definition 'Slice' and 'andWithPrev' processing.
% 10/03/2015 Added Shape definition 'Section' for sections of PData.
% 4/09/2015 Added Shape definition for 'Cone'.
% 1/05/2015 Added checkerboarding of PData.

% Define some non-structure values for faster execution.
numRows = PData.Size(1);
numCols = PData.Size(2);
if size(PData.Size,2)==3
    numSecs = PData.Size(3);
else
    numSecs = 1;
end
xOrg = PData.Origin(1);
yOrg = PData.Origin(2);
zOrg = PData.Origin(3);

% Evaluate pixel deltas
%  - if PDelta attribute exists, use it to set pdelta values.
if isfield(PData,'PDelta') && ~isempty(PData.PDelta)
    if isfield(PData,'Coord')
        switch PData.Coord
            case 'rectangular'
                if size(PData.PDelta,2) ~= 3
                    error('computeRegions: PData(x).PDelta number of columns must be 3.');
                end
                pdeltaX = PData.PDelta(1);
                pdeltaY = PData.PDelta(2);
                pdeltaZ = PData.PDelta(3);
            case 'polar'
                if size(PData.PDelta,2) == 3
                    if PData.PDelta(3) ~= 0
                        error('computeRegions: PData(x).Coord = ''polar'', but PData(x).PDelta(3)~=0.');
                    end
                elseif size(PData.PDelta,2) ~= 2
                    error('computeRegions: PData(x).Coord = ''polar'', but size(PData(x).PDelta,2)~=2 or 3.');
                end
                pdeltaT = PData.PDelta(1);
                pdeltaR = PData.PDelta(2);
            case 'cylindrical'
                if size(PData.PDelta,2) ~= 3
                    error('computeRegions: PData(x).Coord = ''cylindrical'', but size(PData(x).PDelta,2)~=3.');
                end
                pdeltaT = PData.PDelta(1);
                pdeltaR = PData.PDelta(2);
                pdeltaZ = PData.PDelta(3);
            case 'spherical'
                if size(PData.PDelta,2) ~= 3
                    error('computeRegions: PData(x).Coord = ''spherical'', but size(PData(x).PDelta,2)~=3.');
                end
                pdeltaA = PData.PDelta(1);
                pdeltaE = PData.PDelta(2);
                pdeltaR = PData.PDelta(3);
            otherwise
                error('computeRegions: Unrecognized PData(x).Coord string.');
        end
    else
        % Default to rectangular coordinates if no Coord attribute
        if size(PData.PDelta,2) ~= 3
            error('computeRegions: size(PData(x).PDelta,2) for ''rectangular'' Coord. must be 3.');
        end
        PData.Coord = 'rectangular';
        pdeltaX = PData.PDelta(1);
        pdeltaY = PData.PDelta(2);
        pdeltaZ = PData.PDelta(3);
    end
elseif isfield(PData,'pdelta') && ~isempty(PData.pdelta)
% if no PDelta and PData.pdelta is set, it implies all non-specified pdelta? are equal to its value.
    if ~strcmp(PData.Coord,'rectangular')
        error('computeRegions: For PData(n).pdelta set, PData(n).Coord must be ''rectangular''.');
    end
    pdelta = PData.pdelta;
    pdeltaX = pdelta;
    pdeltaY = pdelta;
    pdeltaZ = pdelta;
else
    % if PDelta and pdelta are missing, at the least, pdeltaX and pdeltaZ must be specified.
    if ~isfield(PData,'pdeltaX') || isempty(PData.pdeltaX)
        error('computeRegions: PData(x).pdeltaX must be specified if PData(x).PDelta or pdelta missing.');
    elseif ~isfield(PData,'pdeltaZ') || isempty(PData.pdeltaZ)
        error('computeRegions: PData(x).pdeltaZ must be specified if PData(x).PDelta or pdelta missing.');
    end
    % Read specified pdeltaX,Y,Z.
    pdeltaX = PData.pdeltaX;
    if isfield(PData,'pdeltaY')&&(~isempty(PData.pdeltaY)), pdeltaY = PData.pdeltaY; end
    pdeltaZ = PData.pdeltaZ;
    PData.Coord = 'rectangular';
end

% Check for region to compute. If no Region defined in PData, create Region for the entire
% PData array with Region.Shape.Name = 'PData', UNLESS SFormat exists (older region definition).
% FOR BACKWARDS COMPATIBILITY (SFormat should not be used for new scripts) if SFormat exists,
% create Region structures based on its definition.
if isfield(PData,'Region')&&(~isempty(PData.Region))
    Region = PData.Region;
    numRegs = size(PData.Region,2);
else % try to program Region(s) based on SFormat for backwards compatibility.
    try
        SFormat = evalin('caller','SFormat');
    catch  % SFormat not found
        SFormat = [];
    end
    if ~isempty(SFormat)
        if isfield(PData,'sFormat'), SF = SFormat(PData.sFormat);
        else SF = SFormat(1);
        end
        PData.Coord = 'rectangular';  % only rectangular coord supported for SFormat.
        % Create Region structures based on SFormat
        switch SF.scanFormat
            case 'RLIN'
                if SF.numRays == 1
                    Region = struct(...
                        'Shape',struct('Name','Rectangle', ...
                                       'Position',[0,0,SF.startDepth], ...
                                       'width',SF.rayDelta, ...
                                       'height',SF.endDepth-SF.startDepth), ...
                        'numPixels',[],...
                        'PixelsLA',[]);
                    numRegs = 1;
                else
                    Region = repmat(struct(...
                        'Shape',struct('Name','Rectangle', ...
                                       'Position',[0,0,SF.startDepth], ...
                                       'width',SF.rayDelta, ...
                                       'height',SF.endDepth-SF.startDepth), ...
                        'numPixels',[],...
                        'PixelsLA',[]),1,SF.numRays);
                    numRegs = SF.numRays;
                    for nr = 1:numRegs
                        Region(nr).Shape.Position(1) = SF.FirstRayLoc(1) + (nr-1)*SF.rayDelta;
                    end
                end
            case 'SLIN'
                if SF.numRays == 1
                    Region = struct(...
                        'Shape',struct('Name','Parallelogram', ...
                                       'Position',[0,0,SF.startDepth], ...
                                       'width',SF.rayDelta, ...
                                       'height',SF.endDepth-SF.startDepth, ...
                                       'angle',SF.theta), ...
                        'numPixels',[],...
                        'PixelsLA',[]);
                    numRegs = 1;
                else
                    Region = repmat(struct(...
                        'Shape',struct('Name','Rectangle', ...
                                       'Position',[0,0,SF.startDepth], ...
                                       'width',SF.rayDelta, ...
                                       'height',SF.endDepth-SF.startDepth, ...
                                       'angle',SF.theta), ...
                        'numPixels',[],...
                        'PixelsLA',[]),1,SF.numRays);
                    numRegs = SF.numRays;
                    for nr = 1:numRegs
                        Region(nr).Shape.Position(1) = SF.FirstRayLoc(1) + (nr-1)*SF.rayDelta;
                    end
                end
            case 'CLIN'
                if SF.numRays == 1
                    Region = struct(...
                        'Shape',struct('Name','Sector', ...
                                       'Position',[0,0,0], ...
                                       'r1',SF.radius+SF.startDepth, ...
                                       'r2',SF.radius+SF.endDepth, ...
                                       'angle',SF.rayDelta, ...
                                       'steer',0), ...
                        'numPixels',[],...
                        'PixelsLA',[]);
                    numRegs = 1;
                else
                    Region = repmat(struct(...
                        'Shape',struct('Name','Sector', ...
                                       'Position',[0,0,0], ...
                                       'r1',SF.radius+SF.startDepth, ...
                                       'r2',SF.radius+SF.endDepth, ...
                                       'angle',SF.rayDelta, ...
                                       'steer',0), ...
                        'numPixels',[],...
                        'PixelsLA',[]),1,SF.numRays);
                    numRegs = SF.numRays;
                    for nr = 1:numRegs
                        Region(nr).Shape.steer = SF.theta + (nr-1)*SF.rayDelta;
                    end
                end
            case 'VAPX'
                if SF.numRays == 1
                    Region = struct(...
                        'Shape',struct('Name','SectorFT', ...
                                       'Position',[0,0,-SF.radius], ...
                                       'z',SF.startDepth, ...
                                       'r',SF.radius+SF.endDepth, ...
                                       'angle',SF.rayDelta, ...
                                       'steer',0), ...
                        'numPixels',[],...
                        'PixelsLA',[]);
                    numRegs = 1;
                else
                    Region = repmat(struct(...
                        'Shape',struct('Name','SectorFT', ...
                                       'Position',[0,0,-SF.radius], ...
                                       'z',SF.startDepth, ...
                                       'r',SF.radius+SF.endDepth, ...
                                       'angle',SF.rayDelta, ...
                                       'steer',0), ...
                        'numPixels',[],...
                        'PixelsLA',[]),1,SF.numRays);
                    numRegs = SF.numRays;
                    for nr = 1:numRegs
                        Region(nr).Shape.steer = SF.theta + (nr-1)*SF.rayDelta;
                    end
                end
            case 'PYRM'
                if (SF.numRays ~= 1)||((SF.FirstRayLoc(1)+SF.FirstRayLoc(2)+SF.FirstRayLoc(3))~=0)
                    error('computeRegions: ''PYRM'' SFormat must have 1 ray at FirstRayLoc(0,0,0).');
                end
                Region = struct(...
                    'Shape',struct('Name','Pyramid',...
                                   'Position',[0,0,-SF.radius],...
                                   'angle',SF.theta,...
                                   'z1',SF.radius,...
                                   'z2',SF.radius+SF.endDepth,...
                                   'Steer',[0,0]), ...
                    'numPixels',[],...
                    'PixelsLA',[]);
                numRegs = 1;
            case 'CONE'
                if (SF.numRays ~= 1)||((SF.FirstRayLoc(1)+SF.FirstRayLoc(2)+SF.FirstRayLoc(3))~=0)
                    error('computeRegions: ''CONE'' SFormat must have 1 ray at FirstRayLoc(0,0,0).');
                end
                Region = struct(...
                    'Shape',struct('Name','Cone',...
                                   'Position',[0,0,-SF.radius],...
                                   'angle',SF.theta,...
                                   'z1',SF.radius,...
                                   'z2',SF.radius+SF.endDepth,...
                                   'Steer',[0,0]), ...
                    'numPixels',[],...
                    'PixelsLA',[]);
                numRegs = 1;
        end
    else
        Region = struct('Shape',struct('Name','PData'),'numPixels',[],'PixelsLA',[]);
        numRegs = 1;
    end
end

% Calculate Region structures depending on 'Region.Shape.Name'.
for nr = 1:numRegs
    if isfield(Region(nr),'Shape')&&(~isempty(Region(nr).Shape))
        Shape = Region(nr).Shape;
    else
        fprintf(2,'computeRegions: PData.Region(%d) has no Shape attribute.\n',nr);
    end
    switch Shape.Name

        case {'Rectangle','Parallelogram'}  % Parallelogram/Rectangle shape has the following fields:
            % Name          string       % 'Parallelogram' or 'Rectangle'
            % Position      [1x3 double] % x,y,z coordinate of center of top segment
            % width         [double]     % width of top segment
            % height        [double]     % z dimension
            % angle         [double]     % angle between z axis and left side (default = 0)
            % - Check for necessary attributes.
            if isfield(PData,'Coord')&&~strcmp(PData.Coord,'rectangular')
                error('computeRegions: ''Rectangle|Parallelogram'' shape in Region(%d) only available in ''rectangular'' coordinates.',nr);
            end
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Parallelogram|Rectangle'' shape in PData.Region(%d).',nr);
            end
            if ~isfield(Shape,'width')||isempty(Shape.width)
                error('computeRegions: Missing width attribute for ''Parallelogram|Rectangle'' shape in PData.Region(%d).',nr);
            end
            if ~isfield(Shape,'height')||isempty(Shape.height)
                error('computeRegions: Missing height attribute for ''Parallelogram|Rectangle'' shape in PData.Region(%d).',nr);
            end
            if ~isfield(Shape,'angle')||isempty(Shape.angle)
                Shape.angle = 0;  % default is rectangle shape
            end
            % - Construct X and Z arrays containing PData pixel coordinates.
            X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
            X = repmat(X,numRows,1);
            Z = (zOrg:pdeltaZ:(zOrg+(numRows-1)*pdeltaZ))';
            Z = repmat(Z,1,numCols);
            P = logical(X); % predefine logical array same size as X and Z
            % Construct logical array same size as PData that
            % indicates which points in PData are part of image.
            xleft = Shape.Position(1) - Shape.width/2;
            % Check to see if xleft is very close to a previous xright. If so set equal to prevent gaps.
            if exist('xright','var')&&abs(xleft-xright)<0.01, xleft=xright; end
            xright = xleft + Shape.width;
            top = Shape.Position(3); % top is same as z position.
            bottom = top + Shape.height;
            angle = Shape.angle;
            for i = 1:numRows
                zDistFromTop = (zOrg+(i-1)*pdeltaZ)-top;
                P(i,:) = (X(i,:) > xleft + (zDistFromTop*tan(angle)))&...
                         (X(i,:) <= xright + (zDistFromTop*tan(angle)))& ...
                         (Z(i,:) > top)&(Z(i,:) <= bottom);
            end
            PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
            if ~isempty(PixelsLA)  % check for no points in this region.
                Region(nr).numPixels = size(PixelsLA,1);
                Region(nr).PixelsLA = int32(PixelsLA);
            else
                Region(nr).numPixels = 0;
            end

        case {'Sector','Circle','Annulus'}  % Sector, Circle and Annulus shape
            % Name          string       % 'Sector'|'Circle'|'Annulus
            % Position      [1x3 double] % x,y,z coordinate of apex or center
            % r             [double]     % radius for circle shape
            % r1            [double]     % inner radius for sector or annulus
            % r2            [double]     % outer radius for sector or annulus
            % angle         [double]     % angle between left and right side of sector
            % steer         [double]     % angle between z axis and centerline of sector [default 0]
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Sector|Circle|Annulus'' shape in PData.Region(%d).',nr);
            end
            if (strcmp(Shape.Name,'Sector')||strcmp(Shape.Name,'Annulus'))&&(~isfield(Shape,'r1')||isempty(Shape.r1))
                error('computeRegions: Missing r1 attribute for ''Sector|Annulus'' shape in PData.Region(%d).',nr);
            end
            if (strcmp(Shape.Name,'Sector')||strcmp(Shape.Name,'Annulus'))&&(~isfield(Shape,'r2')||isempty(Shape.r2))
                error('computeRegions: Missing r2 attribute for ''Sector|Annulus'' shape in PData.Region(%d).',nr);
            end
            if strcmp(Shape.Name,'Sector')||strcmp(Shape.Name,'Annulus')
                if Shape.r1 >= Shape.r2, error('computeRegions: r1 must be less than r2 for ''Sector'' or ''Annulus'' shape in PData.Region(%d).',nr); end
            end
            if strcmp(Shape.Name,'Sector')&&(~isfield(Shape,'angle')||isempty(Shape.angle))
                error('computeRegions: Missing angle attribute for ''Sector'' shape in PData.Region(%d).',nr);
            end
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
            if strcmp(PData.Coord,'rectangular')
                % - Construct X and Z arrays containing PData pixel coordinates.
                X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
                X = repmat(X,numRows,1);
                Z = (zOrg:pdeltaZ:(zOrg+(numRows-1)*pdeltaZ))';
                Z = repmat(Z,1,numCols);
                % - Compute a polar coordinate array from X and Z with apex at Shape.Position.
                [Theta,R] = cart2pol(Z-Shape.Position(3),X-Shape.Position(1));
            elseif strcmp(PData.Coord,'polar')
                if ~isequal(Shape.Position,PData.Origin)
                    error('computeRegions: With PData.Coord=''polar'', Shape.Position must equal PData.Origin for ''Sector'' family shapes.');
                end
                % - Construct Theta and R arrays based on PData deltas and size.
                Theta = -pdeltaT*(numCols-1)/2:pdeltaT:pdeltaT*(numCols-1)/2;
                Theta = repmat(Theta,numRows,1);
                R = (0:pdeltaR:((numRows-1)*pdeltaR))';
                R = repmat(R,1,numCols);
            else error('computeRegions: only ''rectangular'' or ''polar'' PData.Coord supported for ''Sector'',''Annulus'' or ''Circle'' shape.');
            end
            % - Define angles of sides from z axis.
            theta1 = Shape.steer - Shape.angle/2;
            % Check to see if theta1 is very close to a previous theta2. If so set equal to prevent gaps.
            if exist('theta2','var')&&abs(theta1-theta2)<0.001, theta1=theta2; end
            theta2 = Shape.steer + Shape.angle/2;
            r1 = Shape.r1;
            r2 = Shape.r2;
            if theta1>(-pi)
                if theta2<pi, P = (Theta>=theta1)&(Theta<=theta2)&(R>r1)&(R<=r2);
                else P = ((Theta>theta1)&(Theta<=pi)&(R>r1)&(R<=r2))| ...
                         ((Theta<(theta2-2*pi))&(Theta>(-pi))&(R>r1)&(R<=r2));
                end
            else
                P = ((Theta>theta1+2*pi)&(Theta<=pi)&(R>r1)&(R<=r2))| ...
                    ((Theta<(theta2))&(Theta>(-pi))&(R>r1)&(R<=r2));
            end
            PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
            if ~isempty(PixelsLA)  % check for no points in this region.
                Region(nr).numPixels = size(PixelsLA,1);
                Region(nr).PixelsLA = int32(PixelsLA);
            else
                Region(nr).numPixels = 0;
            end

        case 'SectorFT'  % Sector shape with flat top parallel to x axis
            % Name          string       % 'Sector'
            % Position      [1x3 double] % x,y,z coordinate of apex
            % z             [double]     % z coor. of flat top of sector
            % r             [double]     % radius to bottom arc
            % angle         [double]     % angle between left and right side
            % steer         [double]     % angle between z axis and centerline of sector [default 0]
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'z')||isempty(Shape.z)
                error('computeRegions: Missing z attribute for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'r')||isempty(Shape.r)
                error('computeRegions: Missing r attribute for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'angle')||isempty(Shape.angle)
                error('computeRegions: Missing angle attribute for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if (Shape.angle > pi)
                error('computeRegions: angle > pi for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'steer')||isempty(Shape.steer)
                Shape.steer = 0;
            elseif (Shape.steer+Shape.angle/2 > pi)||(Shape.steer-Shape.angle/2<(-pi))
                error('computeRegions: steer+angle/2 > pi or steer-angle/2 < -pi for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            % Check that r is greater than distance from apex to center of top.
            if Shape.r < (Shape.z-Shape.Position(3))/cos(Shape.steer)
                error('computeRegions: r must be greater than distance from apex to center of top for ''SectorFT'' shape in PData.Region(%d).\n',nr);
            end
            if strcmp(PData.Coord,'rectangular')
                if (numSecs > 1)  % Are we creating an XZ shape in a 3D PData?
                    % Is there an XZ plane with Y==0?
                    yzero = -1;   % yzero is row index of PData that gives plane at Y==0
                    for i = 1:PData.Size(1) % for all Y positions, check if Y==0
                        if abs(PData.Origin(2) - (i-1)*PData.PDelta(2))<0.001, yzero=i; break, end
                    end
                    if yzero == -1
                        error('computeRegions: ''SectorFT'' shape in 3D PData requires that an XZ plane exists at Y=0.\n');
                    end
                    % - Construct X and Z arrays containing coordinates of pixels in PData Y=0 plane.
                    X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
                    X = repmat(X,numSecs,1);
                    Z = (zOrg:pdeltaZ:(zOrg+(numSecs-1)*pdeltaZ))';
                    Z = repmat(Z,1,numCols);
                else  % shape is in 2D PData
                    % - Construct X and Z arrays containing PData pixel coordinates.
                    X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
                    X = repmat(X,numRows,1);
                    Z = (zOrg:pdeltaZ:(zOrg+(numRows-1)*pdeltaZ))';
                    Z = repmat(Z,1,numCols);
                end
                % - Define a polar coordinate array the same size as PData with R=0 at apex.
                [Theta,R] = cart2pol(Z-Shape.Position(3),X-Shape.Position(1));
            elseif strcmp(PData.Coord,'polar')
                if ~isequal(Shape.Position,PData.Origin)
                    error('computeRegions: With PData.Coord=''polar'', Shape.Position must equal PData.Origin for ''SectorFT'' shape.');
                end
                % - Construct Theta and R arrays based on PData deltas and size.
                Theta = -pdeltaT*(numCols-1)/2:pdeltaT:pdeltaT*(numCols-1)/2;
                Theta = repmat(Theta,numRows,1);
                R = (0:pdeltaR:((numRows-1)*pdeltaR))';
                R = repmat(R,1,numCols);
                % - Define rectangular coord arrays the same size as PData.
                [X,Z] = pol2cart(Theta-(pi/2),R);
                Z = -Z;
            else error('computeRegions: only ''rectangular'' or ''polar'' PData.Coord supported for ''SectorFT'' shape.');
            end
            % - Define angles of sides from z axis.
            thetal = Shape.steer - Shape.angle/2;
            if thetal<(-pi), thetal = thetal + 2*pi; end
            if exist('thetar','var')&&(abs(thetar-thetal)<0.001) % if new thetal essentially equal to old thetar
                thetal = thetar; % This insures there is no tiny gap between regions defined serially.
            end
            thetar = Shape.steer + Shape.angle/2;
            if thetar>pi, thetar = thetar - 2*pi; end
            z = Shape.z; % - zOrg;
            r = Shape.r;
            P = (Theta>thetal)&(Theta<=thetar)&(Z>z)&(R<=r);
            PixelsLA = find(P);  % find the non-zero linear indices of P.
            if numSecs > 1 % for 2d shape in 3D volume, compute 3D linear indices
                [iz,ix] = ind2sub([numSecs,numCols],PixelsLA);
                PixelsLA = (ix-1)*numRows + yzero + (iz-1)*(numRows*numCols) - 1;
            else
                PixelsLA = PixelsLA - 1; % make indices from 0 instead of 1
            end
            if ~isempty(PixelsLA)  % check for no points in this region.
                Region(nr).numPixels = size(PixelsLA,1);
                Region(nr).PixelsLA = int32(PixelsLA);
            else
                Region(nr).numPixels = 0;
            end

        case {'Trapezoid','Triangle'}
            % Trapezoid shape, with top and bottom parallel to x axis, or triangle.
            % - If triangle, either top or bottom length is 0.
            % Name          string       % 'Trapezoid' or 'Triangle'
            % Position      [1x3 double] % x,y,z coordinate of center of top
            % height        [double]     % height (bottom to top)
            % top           [double]     % length of top segment
            % bottom        [double]     % length of bottom segment (base)
            % steer
            if isfield(PData,'Coord')&&~(strcmp(PData.Coord,'rectangular')||strcmp(PData.Coord,'polar'))
                error('computeRegions: ''Trapezoid|Triangle'' shape in Region(%d) only available in ''rectangular'' or ''polar'' coordinates.',nr);
            end
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Trapezoid'' or ''Triangle'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'height')||isempty(Shape.height)
                error('computeRegions: Missing height attribute for ''Trapezoid'' or ''Triangle'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'top')||isempty(Shape.top)
                error('computeRegions: Missing top attribute for ''Trapezoid'' or ''Triangle'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'bottom')||isempty(Shape.bottom)
                error('computeRegions: Missing bottom attribute for ''Trapezoid'' or ''Triangle'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'steer')||isempty(Shape.steer)
                Shape.steer = 0;
            end
            if strcmp(Shape.Name,'Triangle')&&((Shape.top*Shape.bottom)~=0)
                error('computeRegions: ''Triangle'' shape in PData.Region(%d) must have either top or bottom length equal to zero.\n',nr);
            end
            if strcmp(PData.Coord,'rectangular')
                % - Construct X and Z arrays containing PData pixel coordinates.
                X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
                X = repmat(X,numRows,1);
                Z = (zOrg:pdeltaZ:(zOrg+(numRows-1)*pdeltaZ))';
                Z = repmat(Z,1,numCols);
            elseif strcmp(PData.Coord,'polar')
                % - Construct Theta and R arrays based on PData deltas and size.
                Theta = -pdeltaT*(numCols-1)/2:pdeltaT:pdeltaT*(numCols-1)/2;
                Theta = repmat(Theta,numRows,1);
                R = (0:pdeltaR:((numRows-1)*pdeltaR))';
                R = repmat(R,1,numCols);
                [X,Z] = pol2cart(Theta,R);
            else
                error('computeRegions: only ''rectangular'' or ''polar'' PData.Coord supported for ''Sector'',''Annulus'' or ''Circle'' shape.');
            end
            % Construct logical array same size as PData that
            % indicates which points in PData are part of region.
            xtleft = Shape.Position(1) - Shape.top/2;
            ztleft = Shape.Position(3); % top is same as z position.
            xtright = xtleft + Shape.top;
            ztright = ztleft;
            xbleft = Shape.Position(1) - Shape.bottom/2 + Shape.height*tan(Shape.steer);
            zbleft = Shape.Position(3) + Shape.height;
            xbright = xbleft + Shape.bottom;
            zbright = zbleft;
            % Compute cross product between side line segments and PData points to decide if point
            % falls on left or right of side line.
            R = (xtleft - xbleft)*(Z - ztleft) - (ztleft - zbleft)*(X - xtleft); % positive if on right
            L = (xtright - xbright)*(Z - ztright) - (ztright - zbright)*(X - xtright); % neg. if on left
            P = (R > 0)&(L < 0)&(Z > ztleft)&(Z < zbleft);
            PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
            if ~isempty(PixelsLA)  % check for no points in this region.
                Region(nr).numPixels = size(PixelsLA,1);
                Region(nr).PixelsLA = int32(PixelsLA);
            else
                Region(nr).numPixels = 0;
            end

        case 'Pyramid'
            % Pyramidal shape truncated by two cross sections parallel to xy plane.
            %                                      .
            %                 +y              .   /|
            %                  |  +x   .         / |
            %                  |/|/             /  |
            %             .    / |             .   |
            %        .        /|/|      .      |   |
            % -z -+--------- / /------.-------------------------- +z
            %        .       |/|/         `   .|   |
            %             .  | /               |   .
            %               /|/|               |  /
            %             -x   |    .          | /
            %                  |          .    |/
            %                 -y               .
            %
            %
            % Name          string       % 'Pyramid'
            % Position      [1x3 double] % 0,0,z coordinate of apex
            % angle         [1x2double]  % centerline to outer surface angle
            % z1            [double]     % optional distance to first cross section
            % z2            [double]     % optional distance to second cross section
            % Steer         [1x2 double] % optional centerline azimuth and elevation
            %
            % Validate attributes
            if isfield(PData,'Coord')&&~strcmp(PData.Coord,'rectangular')
                error('computeRegions: ''Pyramid'' shape in Region(%d) only available in ''rectangular'' coordinates.',nr);
            end
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.Position(1)~=0)||(Shape.Position(2)~=0)
                error('computeRegions: Position(1) and Position(2) must be 0 for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            elseif Shape.Position(3)>zOrg
                error('computeRegions: Position(3) must be <= PData.Origin(3) for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'angle')||isempty(Shape.angle)
                error('computeRegions: Missing angle attribute for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.angle < 0)||(Shape.angle >= pi/2)
                error('computeRegions: angle attribute must be > 0 and < pi/2 for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'z1')||isempty(Shape.z1)
                Shape.z1 = zOrg;
            elseif Shape.z1 < Shape.Position(3)
                error('computeRegions: Shape.z1 value must be > Shape.Position(3) for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'z2')||isempty(Shape.z2)
                Shape.z2 = zOrg + pdeltaZ*numSecs;
            elseif Shape.z2 <= Shape.z1
                error('computeRegions: Shape.z2 value must be > Shape.z1 for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'Steer')||isempty(Shape.Steer)
                Shape.Steer = [0,0];
            elseif (Shape.Steer(1)-Shape.angle<=-pi/2)||(Shape.Steer(1)+Shape.angle>=pi/2)
                error('computeRegions: Shape.Steer(1) value must be between -pi/2+Shape.angle and pi/2-Shape.angle for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.Steer(2)-Shape.angle<=-pi/2)||(Shape.Steer(2)+Shape.angle>=pi/2)
                error('computeRegions: Shape.Steer(1) value must be between -pi/2+Shape.angle and pi/2-Shape.angle for ''Pyramid'' shape in PData.Region(%d).\n',nr);
            end
            % Create X and Y coordinate arrays for the first PData section.
            X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
            Y = (yOrg:-pdeltaY:(yOrg-(numRows-1)*pdeltaY))';
            X = repmat(X,numRows,1);
            Y = repmat(Y,1,numCols);
            % Process PData sections one at a time.
            z0 = zOrg - Shape.Position(3); % distance from apex of pyramid to zOrg (PData.Origin(3)).
            zSec = zOrg;
            xOff = 0; % xOff is the x offset to the centerline of the pyramid
            yOff = 0;
            Region(nr).numPixels = 0;
            Region(nr).PixelsLA = [];
            for m = 1:numSecs
                if (zSec < Shape.z1)||(zSec > Shape.z2)
                    zSec = zSec + pdeltaZ;
                    continue
                end
                % if steering, offset coordinate data in section
                if (Shape.Steer(1)~=0)||(Shape.Steer(2)~=0)
                    xOff = (zSec+z0)*tan(Shape.Steer(1));
                    yOff = sqrt(xOff*xOff + (zSec+z0)*(zSec+z0)) * tan(Shape.Steer(2));
                end
                r = (zSec+z0)*tan(Shape.angle); % r is cntrln to edge distance
                xlft = xOff - r;
                xrht = xOff + r;
                ytop = yOff + r;
                ybtm = yOff - r;
                P = (X>xlft)&(X<=xrht)&(Y>ybtm)&(Y<=ytop);
                PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
                if ~isempty(PixelsLA)  % check for no points in this region.
                    Region(nr).numPixels = Region(nr).numPixels + size(PixelsLA,1);
                    linAdrOff = (m-1)*(numRows*numCols);
                    Region(nr).PixelsLA = cat(1,Region(nr).PixelsLA, (PixelsLA + linAdrOff));
                end
                zSec = zSec + pdeltaZ;
            end
            Region(nr).PixelsLA = int32(Region(nr).PixelsLA);

        case 'Frustum'
            % Frustum shape. This shape is the volume between two rectangles positioned
            %    along the z axis and parallel to the xy plane.
            %                                      .
            %          +y                     .   /|
            %           |  +x          .       l2/ |
            %           |  /    /|              /| |
            %           | /  l1/ |             . | |
            %           |/    /|/|      .      | |/|
            % -z -------/----/-X---------------|-/--------------- +z
            %          /|    |/|/         `   .|/| |
            %         / |  w1| /             w2| | .
            %        /  |    |/|               | |/
            %      -x   |      |     .         | /
            %           |      |            .  |/|
            %           |      |<----height----->|
            %          -y
            %
            % Name          string       % 'Frustum'
            % Position      [1x3 double] % [x y z] of X center point of first rectangle
            % l1            [double]     % length of rectangle 1 (parallel to x axis)
            % w1            [double]     % width of rectangle 1 (parallel to y axis)
            % height        [double]     % distance to second cross section rectangle
            % l2            [double]     % length of rectangle 2
            % w2            [double]     % width of rectangle 2
            % Steer         [1x2 double] % optional centerline azimuth and elevation
            %                              from center of first rectangle.
            % Validate attributes
            if isfield(PData,'Coord')&&~strcmp(PData.Coord,'rectangular')
                error('computeRegions: ''Frustum'' shape in Region(%d) only available in ''rectangular'' coordinates.',nr);
            end
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'l1')||isempty(Shape.l1)
                error('computeRegions: Missing or empty Shape.l1 value for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'w1')||isempty(Shape.w1)
                error('computeRegions: Missing or empty Shape.w1 value for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'height')||isempty(Shape.height)||Shape.height<=0
                error('computeRegions: Shape.height value must exist and be > 0 for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'l2')||isempty(Shape.l2)
                error('computeRegions: Missing or empty Shape.l2 value for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'w2')||isempty(Shape.w2)
                error('computeRegions: Missing or empty Shape.w2 value for ''Frustum'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'Steer')||isempty(Shape.Steer)
                Shape.Steer = [0 0]; % azimuth and elevation set to zero.
            elseif Shape.Steer(1)<=(-pi/2)||Shape.Steer(1)>=(pi/2)
                error('computeRegions: Shape.Steer(1) value for ''Frustum'' shape in PData.Region(%d) must be between -pi/2 and pi/2.\n',nr);
            elseif Shape.Steer(2)<=(-pi/2)||Shape.Steer(2)>=(pi/2)
                error('computeRegions: Shape.Steer(2) value for ''Frustum'' shape in PData.Region(%d) must be between -pi/2 and pi/2.\n',nr);
            end
            % Create X and Y coordinate arrays for the first PData section.
            X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
            Y = (yOrg:-pdeltaY:(yOrg-(numRows-1)*pdeltaY))';
            X = repmat(X,numRows,1);
            Y = repmat(Y,1,numCols);
            % Process PData sections one at a time.
            zSec = zOrg;
            zt = Shape.Position(3);
            zb = Shape.Position(3) + Shape.height;
            xOff = 0;
            yOff = 0;
            Region(nr).numPixels = 0;
            Region(nr).PixelsLA = [];
%             figure;
%             hold on;
            for m = 1:numSecs
                if (zSec < zt)||(zSec > zb)
                    zSec = zSec + pdeltaZ; % increment to next section
                    continue
                end
                xl = -Shape.l1/2 + (zSec-zt)/Shape.height * (Shape.l1-Shape.l2)/2;
                yw = -Shape.w1/2 + (zSec-zt)/Shape.height * (Shape.w1-Shape.w2)/2;
                % if steering, compute offset to centerline of rectangle.
                if (Shape.Steer(1)~=0)||(Shape.Steer(2)~=0)
                    xOff = (zSec-zt)*tan(Shape.Steer(1));
                    yOff = sqrt(xOff*xOff + (zSec-zt)*(zSec-zt)) * tan(Shape.Steer(2));
                end
                % Compute rectangle boundaries, including steering offset
                xm = xl + Shape.Position(1);
                xp = xm + 2*abs(xl) + xOff;
                xm = xm  + xOff;
                ym = yw + Shape.Position(2);
                yp = ym + 2*abs(yw) + yOff;
                ym = ym + yOff;
%                 plot(m,xm,'*');
%                 plot(m,ym,'+');
                P = (X>xm)&(X<=xp)&(Y>ym)&(Y<=yp);
                PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
                if ~isempty(PixelsLA)  % check for no points in this region.
                    Region(nr).numPixels = Region(nr).numPixels + size(PixelsLA,1);
                    linAdrOff = (m-1)*(numRows*numCols);
                    Region(nr).PixelsLA = cat(1,Region(nr).PixelsLA, (PixelsLA + linAdrOff));
                end
                zSec = zSec + pdeltaZ;
            end
            Region(nr).PixelsLA = int32(Region(nr).PixelsLA);
            
        case 'Cone'
            % Truncated conical shape, with top and bottom parallel to x,y plane.
            % Name          string       % 'Cone'
            % Position      [1x3 double] % 0,0,z coordinate of apex
            % angle         [double]     % centerline to outer surface angle
            % z1            [double]     % optional distance to first cross section
            % z2            [double]     % optional distance to second cross section
            % Steer         [double]     % optional centerline azimuth and elevation
            %
            % Validate attributes
            if isfield(PData,'Coord')&&~strcmp(PData.Coord,'rectangular')
                error('computeRegions: ''Cone'' shape in Region(%d) only available in ''rectangular'' coordinates.',nr);
            end
            if ~isfield(Shape,'Position')||isempty(Shape.Position)
                error('computeRegions: Missing Position attribute for ''Cone'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.Position(1)~=0)||(Shape.Position(2)~=0)
                error('computeRegions: Position(1) and Position(2) must be 0 for ''Cone'' shape in PData.Region(%d).\n',nr);
            elseif Shape.Position(3)>0
                error('computeRegios: Position(3) must be <= 0 for ''Cone'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'angle')||isempty(Shape.angle)
                error('computeRegions: Missing angle attribute for ''Cone'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.angle < 0)||(Shape.angle >= pi/2)
                error('computeRegions: angle attribute must be > 0 and < pi/2 for ''Cone'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'z1')||isempty(Shape.z1)
                Shape.z1 = zOrg;
            elseif Shape.z1 < Shape.Position(3)
                error('computeRegions: Shape.z1 value must be > Shape.Position(3) for ''Cone'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'z2')||isempty(Shape.z2)
                Shape.z2 = zOrg + pdeltaZ*numSecs;
            elseif Shape.z2 <= Shape.z1
                error('computeRegions: Shape.z2 value must be > Shape.z1 for ''Cone'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'Steer')||isempty(Shape.Steer)
                Shape.Steer = [0,0];
            elseif (Shape.Steer(1)-Shape.angle<=-pi/2)||(Shape.Steer(1)+Shape.angle>=pi/2)
                error('computeRegions: Shape.Steer(1) value must be between -pi/2+Shape.angle and pi/2-Shape.angle for ''Cone'' shape in PData.Region(%d).\n',nr);
            elseif (Shape.Steer(2)-Shape.angle<=-pi/2)||(Shape.Steer(2)+Shape.angle>=pi/2)
                error('computeRegions: Shape.Steer(1) value must be between -pi/2+Shape.angle and pi/2-Shape.angle for ''Cone'' shape in PData.Region(%d).\n',nr);
            end
            % Create X and Y coordinate arrays for a single PData section.
            X = xOrg:pdeltaX:(xOrg+(numCols-1)*pdeltaX);
            if isempty(X), X = xOrg; end
            Y = (yOrg:-pdeltaY:(yOrg-(numRows-1)*pdeltaY))';
            if isempty(Y), Y = yOrg; end
            X = repmat(X,numRows,1);
            Y = repmat(Y,1,numCols);
            % Process PData sections one at a time.
            z0 = zOrg - Shape.Position(3); % distance from apex of pyramid to zOrg (PData.Origin(3)).
            zSec = zOrg;
            xOff = 0;
            yOff = 0;
            Region(nr).numPixels = 0;
            Region(nr).PixelsLA = [];
            for m = 1:numSecs
                if (zSec < Shape.z1)||(zSec > Shape.z2)
                    zSec = zSec + pdeltaZ;
                    continue
                end
                % if steering, offset coordinate data in section
                if (Shape.Steer(1)~=0)||(Shape.Steer(2)~=0)
                    xOff = (zSec+z0)*tan(Shape.Steer(1));
                    yOff = sqrt(xOff*xOff + (zSec+z0)*(zSec+z0)) * tan(Shape.Steer(2));
                end
                % Convert X and Y to polar coordinates
                [~,R] = cart2pol(X-xOff,Y-yOff);
                r = (zSec+z0)*tan(Shape.angle);
                P = (R<=r);
                PixelsLA = find(P) - 1;  % find the non-zero linear indices of P, starting from 0.
                if size(PixelsLA,1)==1, PixelsLA = PixelsLA'; end
                if ~isempty(PixelsLA)  % check for no points in this region.
                    Region(nr).numPixels = Region(nr).numPixels + size(PixelsLA,1);
                    linAdrOff = (m-1)*(numRows*numCols);
                    Region(nr).PixelsLA = cat(1,Region(nr).PixelsLA, (PixelsLA + linAdrOff));
                end
                zSec = zSec + pdeltaZ;
            end
            Region(nr).PixelsLA = int32(Region(nr).PixelsLA);

        case 'Section'
            % Region consists of one or more contiguous sections of PData.
            if ~isfield(Shape,'Sections')||isempty(Shape.Sections)
                error('computeRegions: Missing Sections attribute for ''Section'' shape in PData.Region(%d).\n',nr);
            elseif length(Shape.Sections) ~= 2
                error('computeRegions: Shape.Sections in PData.Region(%d) must have two elements - [start section, no. of sections].\n',nr);
            elseif (Shape.Section(1) + (Shape.Section(2)-1)) > numSecs
                error('computeRegions: Shape.Sections in PData.Region(%d) specifies section(s) greater than PData.Size(3).\n',nr);
            end
            secSize = numRow * numCols;
            Region(nr).numPixels = secSize * Shape.Section(2);
            startPix = secSize*(Shape.Section(1)-1);
            endPix = startPix + secSize*(Shape.Section(2)) - 1;
            Region(nr).PixelsLA = int32(startPix:endPix);

        case 'Slice'
            % Region is a slice of a 3D PData volume. The slice includes all the voxels in a plane
            % through the PData volume.  The 'andWithPrev' attribute can be used to include only the
            % voxels that are common between the slice and the previous PData.Region().
            if ~isfield(Shape,'Orientation')||isempty(Shape.Orientation)
                error('computeRegions: Missing Orientation attribute for ''Slice'' shape in PData.Region(%d).\n',nr);
            end
            if ~isfield(Shape,'oPAIntersect')||isempty(Shape.oPAIntersect)
                error('computeRegions: Missing oPAIntersect attribute for ''Slice'' shape in PData.Region(%d).\n',nr);
            elseif length(Shape.oPAIntersect) ~= 1
                error('computeRegions: Shape.oPAIntersect in PData.Region(%d) must have single element.\n',nr);
            end
            switch Shape.Orientation
                case 'xz'
                    j = round((yOrg - Shape.oPAIntersect)/pdeltaY) + 1; % row address of xz plane
                    if (j<1)||(j>numRows)
                        error('computeRegions: oPAIntersect of %d in Slice for PData.Region(%d) doesn''t intersect PData volume.\n',Shape.oPAIntersect,nr);
                    end
                    PixelsLA = zeros(numSecs*numCols,1);
                    for m = 1:numSecs
                        startPix = (m-1)*numRows*numCols + j;  % first pixel on row
                        endPix = startPix + numRows*(numCols-1);
                        PixelsLA((1+numCols*(m-1)):(numCols*m)) = startPix:numRows:endPix;
                    end
                    Region(nr).PixelsLA = int32(PixelsLA - 1);
                    Region(nr).numPixels = size(PixelsLA,1);
                case 'yz'
                    j = round((Shape.oPAIntersect - xOrg)/pdeltaX) + 1;  % col address of yz plane
                    if (j<1)||(j>numCols)
                        error('computeRegions: oPAIntersect of %d in Slice for PData.Region(%d) doesn''t intersect PData volume.\n',Shape.oPAIntersect,nr);
                    end
                    PixelsLA = zeros(numSecs*numRows,1);
                    for m = 1:numSecs
                        startPix = (m-1)*numRows*numCols + (j-1)*numRows + 1;
                        endPix = startPix + numRows - 1;
                        PixelsLA((1+numRows*(m-1)):(numRows*m)) = startPix:endPix;
                    end
                    Region(nr).PixelsLA = int32(PixelsLA - 1);
                    Region(nr).numPixels = size(PixelsLA,1);
                case 'xy'
                    j = round((Shape.oPAIntersect - zOrg)/pdeltaZ) + 1;  % page address of xy plane
                    if (j<1)||(j>numSecs)
                        error('computeRegions: oPAIntersect of %d in Slice for PData.Region(%d) doesn''t intersect PData volume.\n',Shape.oPAIntersect,nr);
                    end
                    PixelsLA = zeros(numRows*numCols,1);
                    startPix = (j-1)*numRows*numCols + 1;
                    endPix = startPix + numRows*numCols - 1;
                    PixelsLA(1:(numRows*numCols)) = startPix:endPix;
                    Region(nr).PixelsLA = int32(PixelsLA - 1);
                    Region(nr).numPixels = size(PixelsLA,1);
                otherwise
                    error('computeRegions: Unknown Shape.Orientation in PData.Region(%d).\n',nr);
            end
        case 'PData'    % PData shape or no scan format specified.
            % Is there a shape.Fill specified for rectangular coords?
            if (strcmp(PData.Coord,'rectangular'))&&(isfield(Shape,'Fill')&&~isempty(Shape.Fill))
                switch Shape.Fill
                    case 'checkerboard1'
                        if 2*floor(numRows/2) ~= numRows % if numRows is odd
                            X = 0:((numRows*numCols)-1);
                            Y = X(1:2:size(X,2));
                            Region(nr).numPixels = size(Y,2);
                            Region(nr).PixelsLA = int32(Y); % linear indices, starting from 0
                        else  % numRows is even
                            X = 0:((numRows*numCols)-1);
                            Z = [];
                            for i = 0:(numCols-1)
                                j = i*numRows + 1 + (i - 2*floor(i/2));
                                Y = X(j:2:((i+1)*numRows));
                                Z = [Z,Y];
                            end
                            Region(nr).numPixels = size(Z,2);
                            Region(nr).PixelsLA = int32(Z); % linear indices, starting from 0
                        end
                    case 'checkerboard2'
                        if 2*floor(numRows/2) ~= numRows % if numRows is odd
                            X = 0:((numRows*numCols)-1);
                            Y = X(2:2:size(X,2));
                            Region(nr).numPixels = size(Y,2);
                            Region(nr).PixelsLA = int32(Y); % linear indices, starting from 0
                        else  % numRows is even
                            X = 0:((numRows*numCols)-1);
                            Z = [];
                            for i = 0:(numCols-1)
                                j = i*numRows + 1 + 1-(2*floor(i/2)-i);
                                Y = X(j:2:((i+1)*numRows));
                                Z = [Z,Y];
                            end
                            Region(nr).numPixels = size(Z,2);
                            Region(nr).PixelsLA = int32(Z); % linear indices, starting from 0
                        end
                    otherwise
                        fprintf(2,'computeRegions: Unrecognized Shape.Fill attribute for PData.Region(%d).\n',nr);
                end
            else
                % Specify Region structure for the entire PData array (rectangular, polar or cylindrical coords).
                Region(nr).numPixels = numRows*numCols*numSecs;
                Region(nr).PixelsLA = int32(0:(numRows*numCols*numSecs-1)); % linear indices, starting from 0
            end
        case 'Custom'   % User defined Region specification.
            if ~isfield(Region(nr),'numPixels')||isempty(Region(nr).numPixels)
                fprintf(2,'computeRegions: Missing numPixels specification for PData.Region(%d) with ''custom'' Shape.\n',nr);
            end
            if ~isfield(Region(nr),'PixelsLA')||isempty(Region(nr).PixelsLA)
                fprintf(2,'computeRegions: Missing PixelsLA specification for PData.Region(%d) with ''custom'' Shape.\n',nr);
            end
            if ~isa(Region(nr).PixelsLA, 'int32')
                fprintf(2,'PData.Region(%d).PixelsLA must be of class ''int32''.\n',nr);
            end
        otherwise
           fprintf(2,'computeRegions: Unrecognized Shape.Name specification for PData.Region(%d).\n',nr);

    end  % end of switch for Shape.Name

    % If 'andWithPrev' set, find the common pixels/voxels with a previously computed Region.
    if isfield(Shape,'andWithPrev')&&~isempty(Shape.andWithPrev)
        p = Shape.andWithPrev;
        if (p~=0)
            if nr==1, error('computeRegions: ''andWithPrev'' set for Region(%d) but no previous Region.\n',nr); end
            if p>=nr, error('computeRegions: ''andWithPrev'' value in Region(%d) must be less than current Region index.\n',nr); end
            i = 0;
            j = 1;
            k = 1;
            PixelsLATemp = Region(nr).PixelsLA;
            P = sort(Region(p).PixelsLA); % make sure lin adrs are in ascending order
            C = sort(Region(nr).PixelsLA);
            while (j<=Region(p).numPixels)&&(k<=Region(nr).numPixels)
                if P(j)==C(k)
                    i = i+1;
                    PixelsLATemp(i) = C(k);
                    j = j+1;
                    k = k+1;
                elseif P(j)<C(k)
                    j = j+1;
                else
                    k = k+1;
                end
            end
            Region(nr).PixelsLA = PixelsLATemp(1:i);
            Region(nr).numPixels = i;
        end
    end
end


