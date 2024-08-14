function dfunListTot = spdepfind(exludeScripts,excludeList)
%spdepfind:  dependencies for sptoolbox; prints results to screen
%usage:
%
% >> dependentFunctionsList = spdepfind(exludeScripts,1); %exclude scripts
% >> dependentFunctionsList = spdepfind(exludeScripts,[],excludedFuncListCellArr);

%john flynn 2012

if nargin<1
    exludeScripts = [];
end
if isempty(exludeScripts),
    exludeScripts = 1;
end
if nargin<2
    excludeList = {};
end

%all .m files
mlist = dir('*.m');
mlist = {mlist.name};
if exludeScripts,
    %remove sequence programming scripts from search:
    mlist(strmatch('SETUP',upper(mlist)))=[];
end
mlist = setdiff(mlist,excludeList );
%find sptoolbox

sptbxpath = get_sptoolbox_path;
dfunListTot = {};
disp(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp([mfilename,': checking Signal Processing Toolbox dependencies in directory: '])
disp(pwd)
disp( ' ' )
for fnamek= mlist, %loop over function m-files to check

    dfunList = {};
    %depfun on file

    [dl]= depfun(fnamek, '-toponly' ,'-quiet' );


    %check depfuns for sp toolbox
    ind =  strmatch(sptbxpath,dl );

    dl_sptbx = dl(ind);
    if ~isempty(dl_sptbx),
        disp('- - - - - - - - - - - - - - - - - - - - -')
        disp(['SPTbx Check on file: ',fnamek])

        for j=1:length(dl_sptbx),
            dfunspec=dl_sptbx{j};
            [ppp,dfun]=fileparts(dfunspec );
            dfunList{end+1}=dfun;
            disp(['****** SP Tbx DEPENDENCY:  ===>',dfun])
        end

        dfunListTot =[ dfunListTot , dfunList ];
    end

end
disp(datestr(now))
disp([mfilename,': finished checking Signal Processing Toolbox dependency in directory: '])
disp(pwd)
disp( ' ')
if isempty(dfunListTot),
    disp('No dependencies found in this directory.')
else
    disp('*** Found one or more dependencies in this directory.')

end
disp( ' ')

disp('------------------')


end %main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
