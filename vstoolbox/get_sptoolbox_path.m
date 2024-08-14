function sptbxpath = get_sptoolbox_path(fcand);
%get_sptoolbox_path: get location of signal proc toolbox
%get_sptoolbox_path(spfunc) - subdir of spfunc (must be in SP TBX !! ).

if nargin<1,
    fcand = [];
end
if isempty(fcand),
    findSubdir = 1;
    fcand = 'chebwin';
else
    findSubdir = 0;
end

ppp = which( fcand , '-all' );
plist = strfind(ppp,['toolbox',filesep,'signal']);
%remove user directories from candidate list:
ind = [];
for k=1:length(plist),
    pk = plist{k};
    if ~isempty(pk),
        ind(end+1) = k;
    end
end
if length(ind)>1,
    warning(['more than one SP toolbox path found']);
end
ppp = ppp{ind(1)};

sptbxpath = fileparts(ppp);
if findSubdir,
    sptbxpath = strrep(sptbxpath,['signal',filesep,'signal'],['signal']);
end
