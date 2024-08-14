function Xz = upsamplevs(X,N,phase);
%upsamplevs: zero-fill signal vector or matrix of signal columns, like upsample.m
%usage:

%john flynn 7/18/2012

if nargin<3,phase = [];end
if isempty(phase),phase=0;end


%normal processing:
sX=size(X);
ndim = length(sX);
if ndim>2, error('only supporting vector or matrix input') ,end

if isvec(X),
    L=length(X);
    flipped = sX(2)==1; %want row vec
    if flipped, X=X.' ;end
    Xz=zeros(N,L);
    Xz(phase+1,:)=X;
    Xz=Xz(:).';

    if flipped, Xz=Xz.';end

else
    %matrix

    if 1,
        Xz=zeros( [N,sX] );
        Xz(phase+1,:)=X(:);
        Xz=reshape(Xz,sX.*[N,1]);
    else
        Z = zeros([sX,N-1]);
        Xa=cat(3,X,Z);
        Xas=shiftdim(Xa,2);
        Xz=reshape(Xas,sX.*[Nz N ]);
    end

end



end %main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=isvec(x)
b = min(size(x))==1;
end

