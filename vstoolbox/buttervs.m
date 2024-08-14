function [b,a]=buttervs(Nord,nuv,style);
%buttervs: butterworth designer (verasonics version)
%usage:
% [b,a]=buttervs(Nord,nu,'low');
% [b,a]=buttervs(Nord,[nuLower nuUpper]); %bandpass [
% nu == 1 corresponds to Fs/2 (matlab compatible freq spec).
%.
% Method: Use digital LP prototype (loaded from table).
% Currently limited to finite number of orders avail.
% Then use digital domain transformation to get desired digital filter.
%.
% Note: not using the bilinear transform.
% Ref:
% "DSP" by oppenheim & schafer (older text) and "DSP" by roberts & mullis

%john flynn  7/17/2012

persistent Prototypes
if isempty(Prototypes)
    load butter_prototypes
end
if nargin<3,style = [];end

% get the prototype filter TF
nup = Prototypes(Nord).nuML/2;
Fb = Prototypes(Nord).b;
Fa = Prototypes(Nord).a;

[b,a]=filterTranslate(Fb,Fa,nup,nuv/2,style);

end %main



