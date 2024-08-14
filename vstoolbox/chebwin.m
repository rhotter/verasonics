function w = chebwin(M , ain );
%chebwin: Dolph-Chebychev window (Verasonics version)
% usage: w = chebwin(N,a);
% inputs:
%	N - window size
%	a - sidelobe level in db (def. is 100)
%
% ref: "Fast Transforms", Elliot & Rao
% also, R. Lyons article "Computing Chebyshev Window Sequences"
% on DSPrelated.com, under "tips and tricks, Jan 2008.

% John Flynn 2/21/1994 , 7/17/2012
%modified from dolph.m

if nargin<2,
    ain = [];
end
if isempty(ain),
    ain = 100;
end

a = 10^abs(ain/20);
N = M-1;
k=0:N-1;

beta = (cosh(1/N * acosh(a)));

g= (beta*cos(pi*k/N));

flags =   g <=1.0;

D_k = (( -1).^k).* cos_h( N*acos_h(g,flags),flags) ...
    ./cosh(N*acosh(beta));

w = real(ifft(D_k,N)).';
w(1)=w(1)/2;
w(M) = w(1);

w=w(:)/max(w);

end %%%%%%%%%%%%%%%%%%%%%main

function [y]=acos_h(x,flags)

y=[];
y(flags) = acos(x(flags));
y(~flags) = acosh(x(~flags));

end
%%%%%%%%%%%%%%%%%%
function y=cos_h(x,flags);
y=[];
y(flags) = cos(x(flags));
y(~flags)= cosh(x(~flags));
end


