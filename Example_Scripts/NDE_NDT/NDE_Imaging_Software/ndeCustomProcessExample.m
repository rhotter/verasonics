function [imgData] = ndeCustomProcessExample( RData, imagingParameters)
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: ndeCustomProcessExample.m - Example of custom processing function
% for use with the NDE Imaging Software
%
% Description: An example processing function for use with the "custom
% reconstruction" class of custom mode in the NDE Imaging
% software. Images using singular value decomposition of the frequency
% response matrix.  Includes optional custom user parameters for the number
% of singular values and number of frequency points to evaluate.
%
% Last update:
% 05/12/2021

%get imaging parameters
FrequencyMhz = imagingParameters.FrequencyMhz;

LongitudinalVelMps = imagingParameters.LongitudinalVelMps;

NumElements = imagingParameters.NumElements;

sampleLength = imagingParameters.EndSample(1) - imagingParameters.StartSample(1) +1;

SampleRateMHz = imagingParameters.SampleRateMHz;

frequencyPts = (1:sampleLength)/sampleLength*SampleRateMHz;

ElementPosXMm = imagingParameters.ElementPosXMm;

ElementPosZMm = imagingParameters.ElementPosZMm;

XCoordsMm = imagingParameters.XCoordsMm;

ZCoordsMm = imagingParameters.ZCoordsMm;

Connector = imagingParameters.Connector;


%get custom user parameters
defaultSvals = 4;
%get number of singular values to evaluate
if size(imagingParameters.UserParameter,1) > 0
    numSvals = abs(round(imagingParameters.UserParameter(1).ParamValue));
    
    if numSvals < 1 || numSvals > NumElements
        numSvals = defaultSvals;
        
    end
    
else
    %if not defined default to this value
    numSvals = defaultSvals;
end

%get bandwidth
defaultFPoints = 1;
if size(imagingParameters.UserParameter,1) > 1
     numFPoints = abs(round(imagingParameters.UserParameter(2).ParamValue));
     
else
     numFPoints = defaultFPoints;
  
end

%get index of center frequency
[~, centreFreqPt] = min(abs(frequencyPts - FrequencyMhz));

centreFreqPt = centreFreqPt(1);

%convert RF data from int16 to double
RData = double(RData);

%initialize image data matrix
ImgData = zeros(size(XCoordsMm));


if numFPoints == 1
    
    freqPtArray = centreFreqPt;
    
else
    
    FRange = round((numFPoints-1)/2);
    
    freqPtArray = (centreFreqPt-FRange):1:(centreFreqPt-FRange+numFPoints-1);
    
end

for freqPt = freqPtArray
    
    %building frequency response matrix
    FRMat = zeros(NumElements);
    
    for TX = 1:NumElements
        
        for RX = 1:NumElements
            
            %get channel number for receive element
            RChannel = Connector(RX);
            
            AScan = RData(imagingParameters.StartSample(TX):imagingParameters.EndSample(TX), RChannel );
            
            AScanFFT = fft(AScan);
            
            FRMat(TX, RX) = AScanFFT(freqPt);
            
        end
        
    end
    
    %return singular values and vectors
    [SVecs, ~, SVals] = svd(FRMat);
    
    %wave number
    k = 2*pi*frequencyPts(freqPt)/LongitudinalVelMps*1e3;
    
    for xIndex = 1:size(XCoordsMm,1)
        
        for zIndex = 1:size(XCoordsMm,2)
            
            %absolute distances from pixel to elements
            R = sqrt((ElementPosXMm - XCoordsMm(xIndex, zIndex)).^2 + (ElementPosZMm - ZCoordsMm(xIndex, zIndex)).^2).';
            
            %greens functions between pixel and elements
            G = conj(-1i/4*besselh(0,k*R));
            G = G./norm(G);
            
            ImgData(xIndex, zIndex) = ImgData(xIndex, zIndex) + G'* SVecs(:,1:numSvals)*SVals(:,1:numSvals)'*conj(G);
            
        end
        
    end
    
end

%set image data output
imgData = abs(ImgData);

end
