function [Ap1024]=FocusRnd1024ApeGradient(Trans)
% This function creates a Random Aperture that has a slight higher
% concentration of elements towards the center of the probe than in the
% perifery. The distribution is determined by a series of radii (in mm)
% that determine the number of elements that must be alocated to each area.
% The number of elements in each of the aread is specified in the variable
% "farea". The input for the function is the Trans structure and the output
% of the function is a vector containing 1024 elements with "ones" for the
% elements activated for that particular aperture. Keep in mind that this
% function was developed for the 3 and 8 MHz probes from Vermon which have
% four subapertues (1-256, 257-512, 513-768, 769-1024). The function
% verifies that no two elements connected to the same MUX switch can be
% activated at the same time.
%
% Copyright 2001-2022 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

radius=[1 3 4 7 ];
x1=0;y1=0; % sets the point (coordinates) where the radii are taken into consideration.

x=Trans.ElementPos(:,1);
y=Trans.ElementPos(:,2);
X=reshape(x,32,32);
Y=reshape(y,32,32);


for i=1:1024
    dist(i)=sqrt((x(i)-x1)^2+(y(i)-y1)^2);
end
Dist=reshape(dist,32,32);

delr{1}=find(dist<=radius(1));

tempr2=find(dist<=radius(2));
delr{2}=(tempr2(find(ismember(tempr2,delr{1})==0)));

tempr3=find(dist<=radius(3));
delr{3}=(tempr3(find(ismember(tempr3,tempr2)==0)));

tempr4=find(dist<=radius(4));
delr{4}=(tempr4(find(ismember(tempr4,tempr3)==0)));

farea=[0 16 100 148 256]; %this gives a percentage of usage of the different radius of  0.5714    0.3387    0.2909    0.1742 for radius of 1 3 4 7 WITH 256 ELEMENTS


fel256=0;
fel256rep=0;
fel512=0;
fel512rep=0;
fel768=0;
fel768rep=0;
fel1024=0;
fel1024rep=0;

el256=0;
el256rep=0;
el512=0;
el512rep=0;
el768=0;
el768rep=0;
el1024=0;
el1024rep=0;

j=0;

elem=1;
while size(elem,2)<farea(end)
    j=j+1;
    for i=farea(j)+1:farea(j+1)
        de=0;
        di=0;


        while di~=1
            elem(i)=delr{j}(randi([1 size(delr{j},2)],1,1));

            if elem(i) <=256
                if sum(elem == elem(i))>1 || sum(elem == elem(i)+256)~=0  || sum(elem == elem(i)+512)~=0 || sum(elem == elem(i)+768)~=0
                    fel256rep=fel256rep+1;
                    %                         if fel256rep>5e3
                    %                             radius=radius+1;
                    %                         end

                else
                    di=1;
                    fel256=fel256+1;
                end
            end

            if elem(i)>256 && elem(i) <=512
                if sum(elem == elem(i))>1 || sum(elem == elem(i)-256)~=0  || sum(elem == elem(i)+256)~=0 || sum(elem == elem(i)+512)~=0
                    fel512rep=fel512rep+1;
                    %                         if fel512rep>5e3
                    %                             radius=radius+1;
                    %                         end


                else
                    di=1;
                    fel512=fel512+1;
                end
            end
            if elem(i)>512 && elem(i) <=768
                if sum(elem == elem(i))>1 || sum(elem == elem(i)-512)~=0  || sum(elem == elem(i)-256)~=0 || sum(elem == elem(i)+256)~=0
                    fel768rep=fel768rep+1;
                    %                         if fel768rep>5e3
                    %                             radius=radius+1;
                    %                         end

                else
                    di=1;
                    fel768=fel768+1;
                end
            end
            if elem(i)>768 && elem(i) <=1024
                if sum(elem == elem(i))>1 || sum(elem == elem(i)-768)~=0  || sum(elem == elem(i)-512)~=0 || sum(elem == elem(i)-256)~=0
                    fel1024rep=fel1024rep+1;
                    %                         if fel1024rep>5e3
                    %                             radius=radius+1;
                    %                         end

                else
                    di=1;
                    fel1024=fel1024+1;
                end
            end
%             if fel256rep>9e3 || fel512rep>9e3 || fel768rep>9e3 || fel1024rep>9e3
% %                 radius = radius +1;
%                 %                     fel256rep=0;
%                 %                     fel512rep=0;
%                 %                     fel768rep=0;
%                 %                     fel1024rep=0;
%             end
        end
    end
    %    j=j+1;

end


Ap1024=zeros(1,1024);
for i=1:length(elem)
    Ap1024(elem(i))=1;
end
ap1024=reshape(Ap1024,32,32);


end
% return

%% To test if there are repeated numbers

%  if this number is equal to length(elem) then there where not repeated numbers

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
%
%
%

