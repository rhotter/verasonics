function [b,a]=subspoly(Fb,Fa,Gb,Ga);
%subspoly: substitute polynomial in another.
% Perform the substitution of discrete-time
% transfer functions F(z^-1) with G(Z^-1):
%
% F(z^-1)|
%        |z^-1  = G(Z^-1)
%See Oppenheim & Schafer

%john flynn 7/17/2012

[Fb,Fa]=zpad(Fb,Fa);
[Gb,Ga]=zpad(Gb,Ga);

[zF,pF,kF] = tf2zp(Fb,Fa);
M = size(zF,1);
N = size(pF,1);

if M==N,
    cancelCommon = 1;
else
    error('this version requires equal polynomial orders in num/denom');
end

if cancelCommon,
    [zNumFG,pNumFG,kNumFG]=subspolyone(Gb,Ga,zF,cancelCommon);%subs numerator factors
    [zDenFG,pDenFG,kDenFG]=subspolyone(Gb,Ga,pF,cancelCommon);%subs denom factors
    zTot = [zNumFG];
    pTot = [zDenFG];
    kTot = kF*kNumFG/kDenFG;
else
    %this form has problem with multiple numerically equal zeros/pole
    %inexact cancellation.
    [zNumFG,pNumFG,kNumFG]=subspolyone(Gb,Ga,zF);%subs numerator factors
    [zDenFG,pDenFG,kDenFG]=subspolyone(Gb,Ga,pF);%subs denom factors
    zTot = [zNumFG;pDenFG];
    pTot = [pNumFG;zDenFG];
    kTot = kF*kNumFG/kDenFG;

end


[b,a]= zp2tf(zTot,pTot,kTot);

end %main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,p,k]=subspolyone(Gb,Ga,zF,numOnly)
%subspolyone: substitute G into a first-order term(s) with root(s) zF
%zF may be a vector of roots representing multiple 1st order terms.
%Note: Gb and Ga must be same length.
%
% evaluates this:
%   1-zF*Gb(Z^-1)/Ga(Z^-1)
% Returns in [z,p,k] zero-pole product form
%.
% if numOnly==1:
% (useful for combining with related rational polynomials):
% then returns [z,p,k] zero-pole product form for:
% Ga(Z^-1) - zF*Gb(Z^-1)

[zn,pn]=deal({});
kn=[];
if numOnly,
    GaNull = Ga*0;GaNull(1) = 1;
end
for n=1:length(zF)
    num = Ga - zF(n)*Gb;
    if numOnly,
        [zn{n},pn{n},kn(n)]=tf2zp(num,GaNull);
    else
        [zn{n},pn{n},kn(n)]=tf2zp(num,Ga);
    end
end
k=prod(kn);
z=cat(1,zn{:});
p=cat(1,pn{:});

end

%%%%%%%%%%%%%%%%%%%%%

