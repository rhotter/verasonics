function [Ap1024]=Rnd1024Aperture()
% This function creates a Random Aperture taking into account the
% constrains of the MUX switches for the 1024 MUX UTA.
% concentration of elements towards the center of the probe than in the
% perifery. No input is required for this function . the output
% of the function is a vector containing 1024 elements with "ones" for the
% elements activated for that particular aperture. Keep in mind that this
% function was developed for the 3 and 8 MHz probes from Vermon which have
% four subapertues (1-256, 257-512, 513-768, 769-1024). The function
% verifies that no two elements connected to the same MUX switch can be
% activated at the same time.
%
% Copyright 2001-2022 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

el256=0;
el256rep=0;

el512=0;
el512rep=0;

el768=0;
el768rep=0;

el1024=0;
el1024rep=0;

for i=1:256
    de=0;

    elem(i)=randi([1 1024],1,1);
    while de ~= 1
        if elem(i) <=256
            if sum(elem == elem(i))>1 || sum(elem == elem(i)+256)~=0  || sum(elem == elem(i)+512)~=0 || sum(elem == elem(i)+768)~=0
                elem(i)=randi([1 256],1,1);
                el256rep=el256rep+1;
            else
                de=1;
                el256=el256+1;
            end
        end

        if elem(i)>256 && elem(i) <=512
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-256)~=0  || sum(elem == elem(i)+256)~=0 || sum(elem == elem(i)+512)~=0
                elem(i)=randi([257 512],1,1);
                el512rep=el512rep+1;
            else
                de=1;
                el512=el512+1;
            end
        end

        if elem(i)>512 && elem(i) <=768
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-512)~=0  || sum(elem == elem(i)-256)~=0 || sum(elem == elem(i)+256)~=0
                elem(i)=randi([513 768],1,1);
                el768rep=el768rep+1;
            else
                de=1;
                el768=el768+1;
            end
        end
        if elem(i)>768 && elem(i) <=1024
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-768)~=0  || sum(elem == elem(i)-512)~=0 || sum(elem == elem(i)-256)~=0
                elem(i)=randi([769 1024],1,1);
                el1024rep=el1024rep+1;
            else
                de=1;
                el1024=el1024+1;
            end
        end

    end
end

Ap1024=zeros(1,1024);
for i=1:256
Ap1024(elem(i))=1;
end
ap1024=reshape(Ap1024,32,32);

end



%% To test if there are repeated numbers

 %if this number is equal to length(elem) then there where not repeated numbers

% for i=1:256
%
%     if elem(i) <=256
%         elem256(i)=elem(i);
%     end
%
%     if elem(i)>256 && elem(i) <=512
%         elem256(i)=elem(i)-256;
%     end
%
%     if elem(i)>512 && elem(i) <=768
%         elem256(i)=elem(i)-512;
%     end
%
%     if elem(i)>768 && elem(i) <=1024
%         elem256(i)=elem(i)-768;
%     end
% end
%
%
% [length(unique(elem)) length(unique(elem256))]
% length(unique(mod(elem,256)))
