function [y]=filtfiltvs(b,a,x,varargin);
%filtfiltvs: zero-phase filtering, similar behaviour to filtfilt.
% usage:
% [y]=filtfiltvs(b,a,x);
% Motivated by papoulis gerchberg alg, with "dual" space being reversed time
% rather than freq. domain.
% Not supporting initial conditions specification (could be modified to do
% so).
%Alg iteration applies F-B filtering on this vector, enforcing zeros and
%orig. data each iteration:
% | ---zeros--- , ----decay/warmup---, orig.data , ----decay/warmup----, ---zeros---- |
%
%ref:
% A. Papoulis "Probability Random Variables and Stochastic Processes"

%john flynn 7/17/2012

Kvai = length(varargin);kvai = 1; %use kvai and template below:

if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; Niter = tempVar;

if isempty(Niter),
    Niter = 12;
end

flipped = size(x,2)==1;
if flipped x=x(:); end

n = length(x);

zzz = zeros(1,2*n);
Nz = length(zzz);
xa = [zzz,x,zzz]; %zero pad the data to assumed transient fade

N = length(xa);

xex = xa;
for j=1:Niter,

    xa(Nz+[1:n])=x; %enforce original data, keep "extrapolated"

    yf =  filter(b,a,xa);
    yfr = filter(b,a,yf(end:-1:1));
    yfr = yfr(end:-1:1);

    %enforce outer halves of transient decay data to be zero:
    yfr(1:Nz/2)=0;
    yfr(end-Nz/2+1:end)=0;

    xa= yfr;

end
%final iteration we kept the filtered center data for output.
y = yfr(Nz+[1:n]); %remove zero padding

if flipped , y = y.';end

end% main %%%%%%%%%%%%%

