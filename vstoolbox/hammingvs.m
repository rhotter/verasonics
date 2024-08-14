function w=hammingvs(M,alph);
%hammingvs:  hamming window
%usage:  w=hammingvs(N);
%usage:  w=hammingvs(N,alpha);
% default alpha is =0.54;
%ref:  "Discrete Time Signal Processing", Oppenheim et al

%john flynn 7/21/2012

if nargin<2,alph=[];end
if isempty(alph), alph=0.54;end

if M==1,
    w=1;
else
    N=M-1;
    n = [0:N]';  %need symmetric window so endpoints  match
    w = alph-(1-alph)*cos((2*pi/N)*n);
end

end %main %%%%%%%%%%%%%%%%%%%%%%%

