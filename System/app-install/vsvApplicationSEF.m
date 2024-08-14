function vsvApplicationSEF()
%vsvApplicationSEF Single enttry function for Verasonics applications 
%   
%   To use the SoniVue-Matlab application framework you need to include
%   this function into your project as the single entry function (SEF). The
%   SEF then calls and loads the AppInstallData that was installed using
%   the appInstaller
%
% Version 1.0 | 2020-08-26 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    % First Check whether launcherData.mat exists if it doesn't exist then
    % we haven't created it. this most likely means that this is the first
    % time we are running the application. 
    %
    % There is still the chance that the application 
    
    
    import vsv.apps.install.AppInstallData;
    
    fpath = mfilename('fullpath');
    launcherFile = AppInstallData.appDataFileName();
     
    % info = vsv.apps.install.AppInstallData
    info = AppInstallData.findInstallationFile( fpath, launcherFile );
    
    if ~isempty(info)

        info.startApplication( );
        
    else
        warndlg('Installation is corrupt (cannot find installation info file) please reinstall');
    end         
    
    
end


