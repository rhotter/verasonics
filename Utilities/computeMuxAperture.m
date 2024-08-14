function [aperNum,varargout] = computeMuxAperture(Apod,TransL)
% computeMuxAperture - Convert apodization array to an Aperture array and associated aperNum index value.
%   For the Trans structure provided in the input argument "TransL",
%   compute a HVMux.Aperture column array from the apodization array input
%   argument "Apod". If this Aperture does not match any existing column of
%   Trans.HVMux.Aperture, it is appended to the Trans.HVMux.Aperture array
%   and an aperture number is assigned. If the Aperture already exists, the
%   aperture number of the existing array is returned.
%
%   If the Aperture array has duplicate channels (more than one element connected to
%   the same channel), a warning is issued, but only once after a clear all.  If
%   there are duplicate channels in the Aperture not supported by the MUX circuitry,
%   or the Apod array is not supported by the MUX circuitry, an error is generated.
%
%   Inputs:
%      Apod - Apodization array that is the same size as the number of elements in the transducer.
%      TransL - The Trans structure for the probe being used.  (In the user
%         SetUp script the Trans structure must be fully defined prior to the
%         first call to computeMuxAperture.)
%   Output:
%      aperNum - The number of the aperture (col. no.) created in Trans.HVMux.Aperture.
%   varargout outputs:
%      {1}  - Aper: A column vector of doubles of length Trans.numelements with the channel number
%             corresponding to the element.
%      {2}  - BP: A column vector of bytes containing the bit pattern for the MUX switches. The first
%             byte specifies the number of bytes in the pattern.

persistent  warn

MuxAperture = TransL.HVMux.ApertureES;

if size(Apod,2) ~= TransL.numelements
    error('computeMuxAperture: Apod must be a row vector of doubles whose size = Trans.numelements.');
end

if strcmpi('preset', TransL.HVMux.type)
    error('computeMuxAperture: Probe HVMux type is ''preset''; dynamic mux programming cannot be used.');
end


% Define Aperture
Aper = +logical(abs(Apod))' .* TransL.ConnectorES;
% Check Aper array for compatibility with MUX circuitry.
[~,~,Chnls] = find(Aper);
if length(Chnls)~=length(unique(Chnls))
    if isempty(warn)
        fprintf(2,['computeMuxAperture: Warning - multiple elements in MUX aperture connected to same\n'...
                   'channel.  This warning is only given once after a clear all.\n']);
        warn = 0;
    end
end

% Check if duplicate Aperture exists and if not, assign new aperture number.
aperNum = 0;
for j = 1:size(MuxAperture,2)
    if isequal(Aper, MuxAperture(:,j))
        aperNum = j;
        break;
    end
end
if aperNum == 0 % if duplicate Aperture not found
    aperNum = size(MuxAperture,2) + 1;
    MuxAperture(:,aperNum) = Aper;
    TransL.HVMux.ApertureES = MuxAperture;
    assignin('caller','Trans',TransL);
end
varargout{1} = Aper;
varargout{2} = [];

end

