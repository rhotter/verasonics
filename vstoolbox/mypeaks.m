function [isPeak,npeaks,peaksStruct,peakVals]=mypeaks(x,negCurvatureThresh,minpeakheight);
%mypeaks: find peaks in column data of matrix
%[isPeak,numPeaksVec,peaksStruct,peakVals]= mypeaks(X);
%[isPeak,numPeaksVec,peaksStruct,peakVals]= mypeaks(X,negCurvatureThresh);
%single column X only:
%[isPeak,numPeaksVec,peaksStruct,peakVals]= mypeaks(X,negCurvatureThresh,minPeakHeight);

%john flynn VS 31dec2009

sizex = size(x);
if nargin<2,
    negCurvatureThresh = [];
end
if nargin<3,
    minpeakheight = [];
end

flipped=sizex(1)==1;
if flipped,
    x=x.';
end

[Nr,Nc]=size(x);
switch Nc
    case 1
        dimClass = 'vector';
    otherwise
        dimClass = 'matrix';
end


if isempty(negCurvatureThresh)
    d3 = 1;
else
    d3 = [-diff(x,2,1)]>=negCurvatureThresh;
end

if isempty(minpeakheight),
    d4=1;
else
    d4 = x>minpeakheight;
end

switch dimClass
    case 'matrix'

        d1=x(2:end-1,:)>x(1:end-2,:);
        %inequality handles the case of adjacent identical valued entries forming a peak:
        d2 = x(2:end-1,:)>=x(3:end,:);
        d1_And_d2 = d1&d2;
        isPeak=[zeros(1,Nc);d1_And_d2 & d3 ;zeros(1,Nc)]&d4;
        npeaks =  sum(isPeak,1);

        [pindR,pindC] = find(isPeak);
        upindC = unique(pindC);

        peakVals = cell(1,Nc);
        for k=1:Nc,
            peakVals{k} = [];
        end
        peaksStruct = [];
        for kc=1:Nc,
            c= kc;%upindC(kc); %a specific column
            c_k = find(pindC==c);
            pri = pindR(c_k);

            peakvals = x(pri,c);
            if isempty(peakvals),peakvals = [];end
            peakVals{c}= peakvals;

            peakstr.peakRowInds = pri;
            pci = ones(size(pri))*c;
            if isempty(pci),pci = [];end
            peakstr.peakColInds = pci;
            peakstr.peakVals = peakvals;

            if isempty(peaksStruct),
                peaksStruct=peakstr;
            else
                peaksStruct(end+1) = peakstr;
            end
        end


    case 'vector'

        d1_And_d2 = logical(zeros(size(x)));

        for kc=1:Nc, %loop over columns

            d=diff(x(:,kc));
            dnz=find(abs(d)~=0);

            %"curvature" measure valid for flat regions
            ddgtz=d(dnz)>0;
            ddltz=d(dnz)<=0;

            ddi = ddgtz(1:end-1 )&ddltz(2:end );

            ddii  = find(ddi)+1;
            d1_And_d2(ddii ,:) = true;

        end

        isPeak=[zeros(1,Nc); ones(Nr-2,1)&d3 ;zeros(1,Nc)]& d1_And_d2 &  d4;

        npeaks = sum( isPeak,1);%num peaks in each column

        if nargout>=3, %create peaks info struct
            npset=unique(npeaks);

            for k=1:length(npset), %loop over num peaks
                %indicate columns with numpeaks of interest:
                peakstr.numPeak=npset(k);

                npks_k = npeaks==npset(k);
                pci = find(npks_k);

                % whos pci
                peakstr.peakColInds = pci;
                peakstr.numCol=sum(npks_k);

                if npset(k)==0,
                    %columns with no peaks
                    peakstr.peakRowInds = [];
                    peakstr.peakLocations = ...
                        [NaN+zeros(peakstr.numCol,1),peakstr.peakColInds(:)]';

                else

                    %columns with peaks
                    [row_k,col_k]=find(isPeak(:,npks_k));
                    peakstr.peakRowInds = reshape(row_k,npset(k),length(row_k)/npset(k));
                    aaaa=row_k(:);
                    bbbb= peakstr.peakColInds(col_k(:)) ;
                    %whos aaaa bbbb
                    peakstr.peakLocations = [aaaa(:),bbbb(:)]';

                end

                peaksStruct(k) = peakstr;

            end

        end


        if nargout>=4,

            peakVals = x(isPeak);
        end

    otherwise
        error('bad switch')
end



if flipped
    isPeak=isPeak.';
end

end %main


