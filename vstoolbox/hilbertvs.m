function xa = hilbertvs(x)
%hilbertvs:  hilbert/analytic extraction for vector x
% signature to match Mathworks' hilbert.m
% Ref:
% Stefan Hahn, "Hilbert Transforms in Signal Processing" Artech

%john flynn 7/17/2012

if min(size(x))~=1,
    error('supporting vector input only')
end
A = 2;
B = 0;

L = length(x);

X = fft(x);

G = A*ones(size(X));

if iseven(L)
    k=L/2 + 1 ;
    kz = [1 k];
else
    k= (L+1)/2 + 1 ;
    kz =   1;
end
G(k:end) = B;
G(kz) = (A+B)/2;

xa = ifft(G.*X);

end %main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = iseven(x)
y = mod(x,2)==0;
end
