function moveDopplerPoints()
% Move scatterers along line to simulate Doppler string phantom.
%
% Parameters:
%
% Return value:
%
% Copyright 2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

% Keep track of time step between each call
persistent t
if isempty(t)
    t = 1;
end

if evalin('base','exist(''Media'',''var'')')
    Media = evalin('base','Media');
else
    disp('Media object not found in workplace.');
    return
end

% Compute amount of change in x position from previous (velocity)
%Vx = -3;  %flow to the left (red)
Vx = 3; %flow to the right (blue)  

% Modify x position of all media points
Media.MP(:,1) = (-45:5:45) + Vx*t;
t=t+1;
if t==5
    t=0;
end
    
assignin('base', 'Media', Media);

return
