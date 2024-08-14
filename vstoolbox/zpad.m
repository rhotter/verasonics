function [xp,yp]=zpad(x,y);
%pad vectors to equal length
L = max(length(x),length(y));
[xp,yp]=deal(x,y);
xp(end+1:L)=0;
yp(end+1:L)=0;
end
