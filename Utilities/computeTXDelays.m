function [Delay,FocalPt] = computeTXDelays(TX,varargin)
%
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% Compute transmit delays using TX.Origin, TX.focus and TX.Steer.
% Four cases can be computed, based on Trans.type index:
%   type=0: All element y and z values are zero (linear array).
%   type=1: All element y values are zero (curved linear array).
%   type=2: non-zero x and y values (2D array with all z=0, or 3D array).
%   type=3: Annular array of ring elements, all centered on z axis.
%   type=4: row/column array, with element positions along x and y axes.
%
% Delays are computed for each element of the transducer, assuming there are as many
% transmitters as elements. Delays are adjusted so that the minimum delay is always zero.
% The computed delays are then assigned to the output "Delay" array, with
% the same length and element mapping as TX.Apod.
%
% Geometry
%           +y
%            |  +x
%            |/|/
%            - |
%           /|/|
% -z < - - / /-------------------------------- +z
%          |/|/
%          | -
%         /|/|
%       -x   |
%            |
%           -y
%
%
% Revision history:
% June 2020 VTS-1748 add support for Trans.type = 4: row-column array
% Dec 2018 VTS-1034 Fix bug for plane wave (TX.focus = 0) on 2D arrays
% June 2018 VTS-839 extend use of TX.FocalPt to all types
% June 13, 2018: VTS-804 add support for both "all-element" and "active element"
% versions of Apod array mapping to elements
% Jan 12, 2018: add support for annular array (Trans.type = 3),
% and XYZ FocalPt specification for type 2 using optional field TX.FocalPt
% or TX.FocalPtMm

% Define non-structure elements to speed up execution.  If a script
% uses an array of TX structures and defines a field in some of but and not
% in others, then it will be present but empty in those others.  Therefore
% we have to check for not empty before using the field.
if isfield(TX,'Origin') && ~isempty(TX.Origin)
    Origin = TX.Origin;
else
    Origin = [0,0,0];
end
if isfield(TX,'focus') && ~isempty(TX.focus)
    focus = TX.focus;
else
    focus = 0;
end
if isfield(TX,'Steer') && ~isempty(TX.Steer)
    azimuth = TX.Steer(1);
    elevation = TX.Steer(2);
else
    azimuth = 0;
    elevation = 0;
end

% Get needed objects from caller workspace.
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');

% get Trans.type
if ~isfield(Trans,'type') || isempty(Trans.type)
    error('computeTXDelays: Trans.type must be specified in Trans structure.');
else
    transType = Trans.type;
end

% check for mm units and convert to wavelengths if needed
if ~isfield(Trans,'units') || isempty(Trans.units)
    error('computeTXDelays: Trans.units must be specified as ''mm'' or ''wavelengths''.');
end

% create wavelengths conversion factor in case it is needed
speedOfSound = 1.540;  % default speed of sound in mm/usec
if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
    if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
        speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
    end
end
scaleToWvl = Trans.frequency/speedOfSound;

if strcmp(Trans.units, 'mm')
    % convert the mm values to wavelengths in Trans.ElementPos
    Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
    Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
    Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
    if transType == 3
        % for type 3 the fourth entry is also a distance, not an angle
        Trans.ElementPos(:,4) = Trans.ElementPos(:,4) * scaleToWvl;
    end
elseif ~strcmp(Trans.units, 'wavelengths')
    error('computeTXDelays: unrecognized Trans.units.  Must be ''mm'' or ''wavelengths''.');
end

if isfield(TX,'FocalPtMm') && ~isempty(TX.FocalPtMm)
    % convert to wavelengths for use in this function
    TX.FocalPt = TX.FocalPtMm * scaleToWvl;
end

if isfield(TX,'FocalPt') && ~isempty(TX.FocalPt)
    userFocalPt = 1;  % flag that says we are using TX.FocalPt
    FocalPt = TX.FocalPt;
    % create focus distance from FocalPt position to xdcr origin (Not
    % TX.Origin!!)
    focus = sqrt(sum(FocalPt.*FocalPt));
else
    userFocalPt = 0;
    FocalPt = [0,0,0];
end

% Check for flag set to 'Transmit On All Elements'.  If this flag is set, the minimum delay will be
% set based on all elements transmitting; otherwise, it is set based only on the active aperture.
toae = 0;
if size(varargin,2) > 0
    if strcmpi(varargin{1},'TOAE'), toae = 1; end
end

% Compute transmit delays according to array type and focus value. For a positive TX.focus, the
%   distance from each transducer element to the focal point is calculated in D.  For TX.focus = 0,
%   the flat wavefront is steered around the beam origin or 0 pt.  For TX.focus negative, the distance
%   from the beam origin behind the transducer to each element is computed in D.

% The delays are calculated for all elements.  Then TX.Apod and Trans
% Aperture or Connector are used to select the elements that will be
% included in the output TX.Delay array which must be the same size as Apod
switch transType
    case {0, 1, 2}
        if all(Trans.ElementPos(:,2)==0) % test for all y values zero.
            % All element y values zero, no elevation steering (elevation ignored).
            if all(Trans.ElementPos(:,3)==0) % test for all z values zero.
                % All element z values zero => linear array.
                % This is Trans.type 0.  Both element steering angles will
                % be ignored.
                if userFocalPt
                    % user has defined an explicit focal point location in x
                    % z coordinates. Ignore y coordinate for Trans.type = 0.
                    % Compute distance to focal point from each element.
                    X = Trans.ElementPos(:,1)' - FocalPt(1);
                    Z = Trans.ElementPos(:,3)' - FocalPt(3);
                    D = sqrt(X.*X + Z.*Z);
                else
                    if focus > 0.0  % focus>0 => normal focused transmit beam.
                        FocalPt(1) = Origin(1) + focus * sin(azimuth);
                        FocalPt(2) = 0.0;
                        FocalPt(3) = focus * cos(azimuth);
                        % Compute distance to focal point from each active element.
                        X = Trans.ElementPos(:,1)' - FocalPt(1);
                        D = sqrt(X.*X + FocalPt(3)*FocalPt(3));
                    elseif focus == 0.0   % focus=0 => flat focus
                        % For flat focus, the wavefront is steered around the origin.
                        D = (Trans.ElementPos(:,1)'- Origin(1)) * sin(azimuth);
                    else
                        % Focus < 0 => diverging beam from point behind transducer.
                        FocalPt(1) = Origin(1) + focus * sin(azimuth);
                        FocalPt(2) = 0.0;
                        FocalPt(3) = focus * cos(azimuth);
                        % Compute distance to focal point from each active element.
                        X = Trans.ElementPos(:,1)' - FocalPt(1);
                        D = sqrt(X.*X + FocalPt(3)*FocalPt(3));
                    end
                end
            else
                % Some element z values non-zero => curved array, Trans.type = 1.  For the curved array,
                % the steering direction is defined as the angle of the beam from the
                % perpendicular to the surface of the array at the Origin point. For this
                % calculation, the azimuth directions of the element normals must be defined
                % in Trans.ElementPos(:,4).
                if userFocalPt
                    % user has defined an explicit focal point location in x
                    % z coordinates. Ignore y coordinate for Trans.type = 1.
                    % Compute distance to focal point from each element.
                    X = Trans.ElementPos(:,1)' - FocalPt(1);
                    Z = Trans.ElementPos(:,3)' - FocalPt(3);
                    D = sqrt(X.*X + Z.*Z);
                else
                    % Focal point not specified by user, so create it from
                    % TX fields Origin, Focus, and Steer
                    % - Find the closest array element to the beam origin.
                    if size(Trans.ElementPos,2)<4 %|| ~any(Trans.ElementPos(:,4))
                        error('computeTXDelays: Normals to curved array elements must be specified in Trans.ElementPos(:,4)');
                    end
                    D = (Trans.ElementPos(:,1) - Origin(1)).^2 + (Trans.ElementPos(:,3) - Origin(3)).^2;
                    [d1,i] = min(D);
                    if d1 < 0.005  % if element position essentially same as Origin
                        angle = Trans.ElementPos(i,4);
                    else
                        D(i) = 100;  % cancel out previous minimum
                        [d2,j] = min(D);
                        dtheta = Trans.ElementPos(j,4) - Trans.ElementPos(i,4);
                        angle = Trans.ElementPos(i,4) + dtheta*d1/(d1+d2);
                    end
                    if focus > 0.0  % focus>0 => normal focused transmit beam.
                        FocalPt(1) = Origin(1) + focus * sin(angle + azimuth);
                        FocalPt(2) = 0.0;
                        FocalPt(3) = Origin(3) + focus * cos(angle + azimuth);
                        % Compute distance to focal point from each active element.
                        X = Trans.ElementPos(:,1)' - FocalPt(1);
                        Z = Trans.ElementPos(:,3)' - FocalPt(3);
                        D = sqrt(X.*X + Z.*Z);
                    elseif focus == 0.0 % focus=0 => flat focus
                        % For flat focus, the wavefront is steered around the effective center
                        % of the array (x=0), using a linear with arc length change in delay.
                        radius = Trans.radius;
                        D = radius*Trans.ElementPos(:,4) * sin(azimuth);
                    else
                        % Focus < 0 => diverging beam from point behind transducer.
                        FocalPt(1) = Origin(1) + focus * sin(azimuth);
                        FocalPt(2) = 0.0;
                        FocalPt(3) = focus * cos(azimuth);
                        % Compute distance to focal point from each active element.
                        X = Trans.ElementPos(:,1)' - FocalPt(1);
                        D = sqrt(X.*X + FocalPt(3)*FocalPt(3));
            %             disp('TX.focus < 0 not allowed for curved array.')
                    end
                end
            end
        else
            % x, y and possibly z values => 2D array.
            if userFocalPt
                % user has defined an explicit focal point location in x y
                % z coordinates.
                % Compute distance to focal point from each element.
                X = Trans.ElementPos(:,1)' - FocalPt(1);
                Y = Trans.ElementPos(:,2)' - FocalPt(2);
                Z = Trans.ElementPos(:,3)' - FocalPt(3);
                D = sqrt(X.*X + Y.*Y + Z.*Z);
            else
                % Focal point not specified by user, so create it from
                % TX fields Origin, Focus, and Steer
                if focus > 0.0 % focus>0 => normal focused transmit beam.
                    % For matrix arrays, it is assumed that the beam Origin is located near an element position.
                    % If the beam is steered, the steering angles are with respect to the normal of the closest element.
                    %   The elements azimuth and elevation angles, if provided, are specified in
                    %   Trans.ElementPos(:,4) and Trans.ElementPos(:,5).
                    if size(Trans.ElementPos,2) == 3  % check to see if element orientation angles provided
                        Trans.ElementPos(:,4) = 0;
                        Trans.ElementPos(:,5) = 0;
                    elseif size(Trans.ElementPos,2) == 4 % only azimuth provided
                        Trans.ElementPos(:,5) = 0;
                    end
                    % Focal point not specified by user, so we have to create it
                    % Find the element that is closest to the beam origin. Only the x and y positions are compared.
                    [mn,index] = min(sqrt((Trans.ElementPos(:,1)-Origin(1)).^2 + (Trans.ElementPos(:,2)-Origin(2)).^2));
                    if mn>Trans.spacing
                        fprintf(['computeTXDelays: For 2D arrays with a focused transmit, the beam origin'...
                               'should be near to an element position, otherwise steering may be inaccurate.\n']);
                    end
                    FocalPt(1) = Origin(1) + focus * sin(azimuth + Trans.ElementPos(index,4));
                    FocalPt(2) = Origin(2) + focus * sin(elevation + Trans.ElementPos(index,5));
                    FocalPt(3) = Origin(3) + focus * cos(elevation + Trans.ElementPos(index,5)) * cos(azimuth + Trans.ElementPos(index,4));
                    % Compute distance to focal point from each element.
                    X = Trans.ElementPos(:,1)' - FocalPt(1);
                    Y = Trans.ElementPos(:,2)' - FocalPt(2);
                    Z = Trans.ElementPos(:,3)' - FocalPt(3);
                    D = sqrt(X.*X + Y.*Y + Z.*Z);
                elseif focus == 0.0
                    if any(Trans.ElementPos(:,3) ~=0 )  % test for not all z values zero.
                        if (azimuth~=0.0)||(elevation~=0.0)
                            error('computeTXDelays: 2D array with non-zero z values and focus = 0 must have no steering angles.');
                        end
                        D = zeros(1,Resource.Parameters.numTransmit);
                    else
                        D = (Trans.ElementPos(:,1)' - Origin(1)) * sin(azimuth) + (Trans.ElementPos(:,2)' - Origin(2)) * sin(elevation);
                    end
                else
                    % focus<0  => diverging beam from point behind transducer.
                    FocalPt(1) = Origin(1) + focus * cos(elevation) * sin(azimuth);
                    FocalPt(2) = Origin(2) + focus * sin(elevation);
                    FocalPt(3) = Origin(3) + focus * cos(elevation) * cos(azimuth);
                    % Compute distance to focal point from each active element.
                    X = Trans.ElementPos(:,1)' - FocalPt(1);
                    Y = Trans.ElementPos(:,2)' - FocalPt(2);
                    Z = Trans.ElementPos(:,3)' - FocalPt(3);
                    D = sqrt(X.*X + Y.*Y + Z.*Z);
                end
            end
        end
    case 3 % annular array, with different ElementPos definitions than types 0, 1, 2
        if focus > 0
            % For annular array, transmit delays are based on distance from
            % focal point on the z axis to rc, zc position of each element.
            % There is no steering and origin is assumned to be at the
            % center of the face of the array, so axis of the annular array
            % is the Z axis.
            if userFocalPt == 0
                % Focal point not specified by user, so we have to create
                % it from focus
                FocalPt(1) = 0.0;
                FocalPt(2) = 0.0;
                FocalPt(3) = focus;
            end
            % Compute distance to focal point from each active element.
            if FocalPt(1)~=0 || FocalPt(2)~=0
                error('computeTXDelays: Trans.type 3 for annular arrays requires focal point on Z-axis with X=0 and Y=0.');
            end
            Zp = FocalPt(3) - Trans.ElementPos(:,4)'; % path length component in the Z direction
            Rp = Trans.ElementPos(:,3)'; % path length component in the R direction
            D = sqrt(Zp.*Zp + Rp.*Rp); % actual path length to each element
        else
            error('computeTXDelays: Trans.type 3 for annular arrays requires TX.focus > 0');
        end
    case 4 % row/column array, with element positions along x and y axes. Elements along x axis have
        % Trans.elementLengths parallel to y axis and Trans.ElementPos(:,2)=0. Elements along y axis
        % have Trans.elementLengths parallel to x axis and Trans.ElementPos(:,1)=0;
        if any(and(Trans.ElementPos(:,1),Trans.ElementPos(:,2)))
            error('computeTXDelays: row/col arrays can''t have both Trans.ElementPos x and y non-zero.');
        end
        X = zeros(1,Trans.numelements);
        Y = X;
        if userFocalPt
            % user has defined an explicit focal point location in x y z coordinates.
            % Compute distance to focal point from each element.
            Z = FocalPt(3) - Trans.ElementPos(:,3)';
            % For x elements, ignore y position of FocalPt
            indices = find(Trans.ElementPos(:,1));
            X(indices) = FocalPt(1) - Trans.ElementPos(indices,1)';
            % For y elements, ignore x position of FocalPt
            indices = find(Trans.ElementPos(:,2));
            Y(indices) = FocalPt(2) - Trans.ElementPos(indices,2)';
            D = sqrt(X.*X + Y.*Y + Z.*Z);
        else
            % Focal point not specified by user, so create it from
            % TX fields Origin, Focus, and Steer
            if focus > 0.0
                % Positive transmit focus with Trans.Origin in plane of array and focal
                % point x and y beneath footprint of array elements. If the beam is steered,
                % the steering angles are with respect to a line from the beam origin parallel
                % to the z axis.
                FocalPt(1) = Origin(1) + focus * cos(elevation) * sin(azimuth);
                FocalPt(2) = Origin(2) + focus * sin(elevation);
                FocalPt(3) = Origin(3) + focus * cos(elevation) * cos(azimuth);
                Z = Trans.ElementPos(:,3)' - FocalPt(3);
                % For x elements, ignore y position of FocalPt
                indices = find(Trans.ElementPos(:,1));
                X(indices) = FocalPt(1) - Trans.ElementPos(indices,1)';
                % For y elements, ignore x position of FocalPt
                indices = find(Trans.ElementPos(:,2));
                Y(indices) = FocalPt(2) - Trans.ElementPos(indices,2)';
                D = sqrt(X.*X + Y.*Y + Z.*Z);
            elseif focus == 0
                % Plane wave transmit. In this case, azimuth angle steering, TX.Steer(1), is 
                % applied to elements along the x axis, and elevation steering, TX.Steer(2),
                % is applied to elements along the y axis. The TX.Origin point is ignored.
                % For flat focus, the wavefront is steered around the z axis.
                D = zeros(1,Trans.numelements);
                indices = find(Trans.ElementPos(:,1));
                D(indices) = Trans.ElementPos(indices,1)' * sin(azimuth);
                indices = find(Trans.ElementPos(:,2));
                D(indices) = Trans.ElementPos(indices,2)' * sin(elevation);
            elseif focus < 0
                % Negative focus implies diverging wave from focal point behind transducer.
                FocalPt(1) = Origin(1) + focus * cos(elevation) * sin(azimuth);
                FocalPt(2) = Origin(2) + focus * sin(elevation);
                FocalPt(3) = Origin(3) + focus * cos(elevation) * cos(azimuth);
                Z = Trans.ElementPos(:,3)' - FocalPt(3);
                % For x elements, ignore y position of FocalPt
                indices = find(Trans.ElementPos(:,1));
                X(indices) = FocalPt(1) - Trans.ElementPos(indices,1)';
                % For y elements, ignore x position of FocalPt
                indices = find(Trans.ElementPos(:,2));
                Y(indices) = FocalPt(2) - Trans.ElementPos(indices,2)';
                D = sqrt(X.*X + Y.*Y + Z.*Z);                
            end
        end

    otherwise
        error('computeTXDelays: unrecognized Trans.type value of %d.\n', transType);
end

if FocalPt(3) > 0
    normMethod = 1;
else
    normMethod = 2;
end

% VTS-804 Determine the size of the Apod array and create Delay array to
% match
sizeApod = size(TX.Apod,2); % number of entries in the Apod arrayt
if sizeApod == Trans.numelements
    % This is an all-element Apod array; either HVMux probe using dynamic
    % aperture on non-HVMux probe with all-element mapping (or with no inactive
    % elements in Trans.Connector)
    if isfield(Trans,'HVMux')
        % Check for parallel elements through HVMux switching, and force
        % same delay within each parallel group of elements.
        if isfield(Trans,'ConnectorES')
            % Create the active Aperture array from Trans.Connector qualified by
            % Apod, for HVMux scripts using dynamic mux programming
            LApod = (TX.Apod ~= 0)';
            Aperture = Trans.ConnectorES .* LApod;
            % Aperture now lists the channel each active element is mapped to, and
            % zero for inactive elements
            % If HVMux is selecting elements in parallel, apply delay from lowest
            % element number within the parallel group
            [~,DupIndices] = ismember(Aperture,Aperture); % duplicate channels will get the lowest element index
            D = D(DupIndices); % set delays of parallel elements equal to delay of lowest element.
        else
            % Trans.ConnectorES is requires for HVMux scripts using dynamic
            % mux programming
            error('computeTXDelays: Trans.ConnectorES field is required when using dynamic mux programming.');
        end
    end
    % Initialize output Delay array to all zeros
    Delay = zeros(1,Trans.numelements);
    % set only delays of elements in aperture with non-zero Apod
    ActvIndices = find(TX.Apod);
    if toae
        % normalize over all elements, not just those that are active
        switch normMethod
            case 1
                D = max(D) - D;
                Delay(ActvIndices) = D(ActvIndices);
            case 2
                D = D - min(D);
                Delay(ActvIndices) = D(ActvIndices);
        end
    else
        % toae is false, so normalize only over the active elements
        switch normMethod
            case 1
                Delay(ActvIndices) = max(D(ActvIndices)) - D(ActvIndices);
            case 2
                Delay(ActvIndices) = D(ActvIndices) - min(D(ActvIndices));
        end
    end
elseif sizeApod < Trans.numelements
    % This is an "active aperture" Apod array- either the legacy HVMux
    % programming with precomputed apertures, or a non-muxed probe with
    % some inactive elements in Trans.Connector that have been excluded
    % from Apod
    % Find the Aperture being used and check for correct Apod length
    if isfield(Trans,'HVMux')
        % This is an HVMux probe so TX.aperture is needed to select Aperture
        if ~isfield(TX,'aperture') || isempty(TX.aperture)
            error('computeTXDelays: TX.aperture must be specified for HVMux script using active element Apod.');
        elseif TX.aperture > size(Trans.HVMux.ApertureES,2) || TX.aperture < 1
            error('computeTXDelays: TX.aperture index into Trans.HVMux.ApertureES is out of range.');
        end
        Aperture = Trans.HVMux.ApertureES(:,TX.aperture);
        if sizeApod ~= nnz(Aperture)
            error('computeTXDelays: length of TX.Apod does not match number of active elements in Trans.HVMux.ApertureES.');
        end
    else
        % this is not an HVMux script, so use Trans.Connector
        if isfield(Trans,'ConnectorES')
            Aperture = Trans.ConnectorES(:, 1);
        else
            % we already know Apod is smaller than Trans.numelements, so we
            % cannot create a default Aperture array- this is an error
            error('computeTXDelays: Trans.ConnectorES is required field when using active element Apod mapping.');
        end
        if sizeApod ~= nnz(Aperture)
            error('computeTXDelays: length of TX.Apod does not match number of active elements in Trans.ConnectorES.');
        end

    end
    % We now have valid Aperture and Apod arrays; create Delay from D with
    % same element mapping as Apod
    Elements = find(Aperture); % list of active elements in Aperture
    Delay = zeros(1,sizeApod);  % initialize Delay to all zeros
    if toae==0
        DelayTemp(1:sizeApod) = D(Elements);
        % Normalize delays to delays of active elements in the Aperture array.
        Indices = find(TX.Apod); % get indices of active transducer elements in the Apod array.
        switch normMethod
            case 1
                Delay(Indices) = max(DelayTemp(Indices)) - DelayTemp(Indices);
            case 2
                Delay(Indices) = DelayTemp(Indices) - min(DelayTemp(Indices));
        end
    else
        % toae is true so normalize over all elements
        switch normMethod
            case 1
                D = max(D) - D;
                Delay(1:sizeApod) = D(Elements);
            case 2
                D = D - min(D);
                Delay(1:sizeApod) = D(Elements);
        end
    end
else
    % size of Apod exceeds Trans.numelements; this is an error condition
    error('computeTXDelays: Unrecognized TX.Apod array with length greater than Trans.numelements.');
end
