function [b,a]=filterTranslate(Fb,Fa,nup,nuv,style);
%translate filter prototype to different band
%usage:
% [b,a]=filterTranslate(bp,ap,nup,nuv,style);
% bp,ap: prototype transfer fcn coeff
% nup: prototype norm. freq (0.5==Nyquist, 1.0==Fs)
% nuv: target norm. freq (0.5==Nyquist, 1.0==Fs)

%john flynn 7/19/2012
thetap = nup*2*pi;
if nargin<3,style=[];end
if isempty(style),
    switch length(nuv),
        case 1
            style = 'low';
        case 2
            style='bandpass';
        otherwise
            error('bad switch')
    end
else
    style=lower(style);
    if (length(nuv)~=2)&strcmp(style,'bandpass')
        error('Nu must be vector if style is bandpass'),
    end
    if (length(nuv)~=1)&~strcmp(style,'bandpass')
        error('Nu must be scalar unless style is bandpass'),
    end
end

%transformations.  Opp&Shaf p. 230
%calc the function G(Z^-1)
switch style
    case 'low'
        omegad = (nuv )*2*pi;
        alph=sin((thetap-omegad)/2)/sin((thetap+omegad)/2);
        GofZm1Num = [-alph 1]; %Numerator
        GofZm1Den = [1 -alph]; %Denom

    case 'high'
        omegad = (nuv )*2*pi;
        alph=-cos((thetap+omegad)/2)/cos((thetap-omegad)/2);
        GofZm1Num = -[alph 1]; %Numerator
        GofZm1Den = [1 alph]; %Denom

    case 'stop'
        error('under construction')

    case 'bandpass'
        nuupper=max(abs(nuv));
        nulower=min(abs(nuv));
        omega2 = (nuupper )*2*pi;
        omega1 = (nulower )*2*pi;
        alph = cos((omega2+omega1)/2) / cos((omega2-omega1)/2) ;
        k = tan(thetap/2)*cot((omega2-omega1)/2);
        GofZm1Num = -[ (k-1)/(k+1)  , -2*alph*k/(k+1)  , 1          ]; %Numerator
        GofZm1Den =  [ 1            , -2*alph*k/(k+1)  , (k-1)/(k+1)]; %Denom

    otherwise
        error('bad switch')
end

%substitute the variable z^-1 in F(z^-1) with G(Z^-1)
[b,a]=subspoly(Fb,Fa,GofZm1Num,GofZm1Den);

end %main %%%%%%%%%%%%%%%%%%%%%%

