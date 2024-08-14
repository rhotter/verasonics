function y=sincvs(x);
%sincvs:  sin(x)/x

% john flynn 7/22/2012

y=(1/pi)*sin(pi*x)./x;
zind = x==0;
%detect nan/inf separately.
%want to preserve incoming NaNs.
iind=isinf(y);
y(zind)=1;
y(iind)=1;

