function Waveform = computeTGCWaveform(TGC)
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Compute TGC waveform from TGC structure attributes: CntrlPts, rangemax.
%    TGC.CntrlPts is a vector with 8 control pts.
%    TGC.rangeMax is the maximum range in wavelengths.
% Function first computes rangeMax in number of 800nsec intervals, which is
%    the TGC generator sample period.  It then divides this distance by 7 to
%    find the depth of the control points (1st point is at 0 range).  These
%    points, and a duplicate of the eighth point at max range (511) are
%    used to compute a piecewise cubic Hermite waveform at the TGC sample
%    period.
%
% Updated June 30, 2016 to test for a rangeMax value that exceeds the 512
% sample length of the TGC.Waveform array.  In this case, the control
% points will be evenly distributed over the 512 sample waveform, and the
% value of the eighth control point will be held by the HW system for the
% remainder of the acquisition interval.

% Get center frequency, and convert wavelengths to samples at 800nsec.
fc = evalin('caller', 'Trans.frequency')*1.0e06;
maxSamples = min(510.9, round(2*double(TGC.rangeMax)/800e-09/fc));
% Note the limit of 510.9 ensures that the eighth control point position
% will always be less than the ninth one at 511, since the interopoation
% routine requires that they be monotonic.

% Construct P and C arrays to specify known points on waveform.
P = zeros(1,9);
P(1:8) = 0:(maxSamples/7):maxSamples;
P(9) = 511;
C = zeros(1,9);
C(1:8) = double(TGC.CntrlPts);
C(9) = double(TGC.CntrlPts(8));
X = 0:511;

% Compute the interpolated points.
Waveform = round(interp1(P,C,X,'pchip'));
% figure
% plot(X(1:maxSamples),Waveform(1:maxSamples))

return
