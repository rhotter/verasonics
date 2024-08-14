function movePoints
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% movePoints moves the x location of points in the Media.Pts array in a
% sinusoidal pattern.

persistent tval
delta = 2*pi/30; % change in tval between successive calls
amp = 10;        % amplitude of displacement in wavelengths

if isempty(tval)
    tval = 0;
end

if evalin('base','exist(''Media'',''var'')')
    Media = evalin('base','Media');
else
    disp('Media object not found in workplace.');
    return
end

% Compute amount of change in x position from previous.
old = amp*sin(tval);
tval = tval + delta;
if tval >= 2*pi
    tval = tval - 2*pi;
end
new = amp*sin(tval);
dif = new - old;

% Modify x position of all media points
Media.MP(:,1) = Media.MP(:,1) + dif;

assignin('base', 'Media', Media);

return
