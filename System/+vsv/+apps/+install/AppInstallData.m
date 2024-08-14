classdef (Sealed) AppInstallData < handle
%APPINSTALLDATA Installation information data for vantage applications
%   
%   The AppInstallData is supposed to store information about the
%   application that was installed. The information can then be used to
%   compare the current sw version with the installed version of the app.
%   
%   This can be used to validate whether the current software version is
%   compatible with the version as installed.
%
%   The class also provides several static functions to determine the file
%   location of the AppInfo data based on the appInfo as returne by the
%   matlab.apputil.getInstalledAppInfo. The file name of the app info that
%   contains an instance of this class is constant and is
%   'VsAppInfoFile.mat'. This also means that the installed app or the
%   vantage software, including user software, should not have a mat-file
%   in their path called 'VsAppInfoFile.mat'
%
%   The AppInstallData also contains a function called startApplication
%   that will start the application as defined by the AppInstallData. This
%   way the SEF for the .prj file only needs to contain very limited code
%   and can rely on the AppInstallData to start an application
%
%   The code for the SEF can look like this:
%
%     import vsv.apps.install.AppInstallData;
%     
%     fpath = mfilename('fullpath');
%     launcherFile = AppInstallData.appDataFileName();
%      
%     info = AppInstallData.findInstallationFile( fpath, launcherFile );
%     
%     if ~isempty(info)
% 
%         info.startApplication();
%         
%     else
%         warndlg('Installation is corrupt (cannot find installation info file) please reinstall');
%     end 
%     
%
% Version 1.0 | 2020-08-26 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

% @NOTE:
% This class is public and is used by an install application, all
% modifications must be carefully reviewed in future releases! 
%
% This class should not make any references to other packages and classes. 

    properties(SetAccess=private)
       
        % the file path where the application was installed and was called
        % the first time when installed
        InstallFolder (1,:) char
        
        % the installatin version
        InstallVersion (1,3) int16 = zeros(1, 3, 'int16');
        
        % the function name of the main entry start function. This should
        % be the full package name such as vsv.apps.launcher.startLaunher()
        StartFunction (1,:) char;
                
        % Allows to add additional path information for vsv external
        % functions
        UserPath {mustBeCellString} = {};
        
        % a user defined activate function
        UserActivate (1,:) char;
        
        % The ID of the installed app
        InstallAppID (1,:) char;
        
        % indicates whether the app install data is in a valid state
        IsValid (1,1) logical = true;
        
        % must be a cell array of verbosity options: { 'nothing', 'print', ...
        % 'warning', 'warndlg' } the user can provide several options for
        % example if he/she wants to print to the console and show the warn
        % dialog at the same time
        Verbosity {mustBeValidVerbosity} = {'print', 'warndlg'};
        
    end
    
    properties(SetAccess=private)
        
        % last error message 
        LastError (1, :) char = '';
        
    end
    
    properties(Access=private)
        
        % internal file location when loaded using from file
        FileLocation (1, :) char = '';
        FileVarName (1, :) char = '';
        
    end
    
    properties(SetAccess=private, Dependent)
        InstallVersionStr (1,:) char;
    end
    
    %% static helper functions 
    % help to locate the install info file 
    methods(Static)

        function file = appDataFileName( )
        % Assumes appInfo as returned by matlab.apputil.getInstalledAppInfo
        % and contains a location and id
        %
        % the returned file will be the file where the additional app info
        % will be stored. The app data that is storred there will be used
        % to compare to the installed version vs current running and
        % executing version
        %
            file = 'VsAppInfoFile.mat';
            
        end
        
        function file = appDataFileFromInfo( appInfo)
        % Assumes appInfo as returned by matlab.apputil.getInstalledAppInfo
        % and contains a location and id
        %
        % the returned file will be the file where the additional app info
        % will be stored. The app data that is storred there will be used
        % to compare to the installed version vs current running and
        % executing version
        %
        
            filename = vsv.apps.install.AppInstallData.appDataFileName();
            file = fullfile( appInfo.location, filename );
             
        end
        
        function [ data, msg] = fromFile(file)
        % returnes the app install data from the given file
        %
        % If file contains several objects that are of type AppInstallData
        % this will return the first parameter found. Msg will be empty
        % in this case
        %
        % @param file - @type char the file to load from
        % @return data - @type AppInstallData the data loaded from file. If
        %                file does not exist or if the file does not contain a
        %                AppInstallData object this will return []
        % @return msg - @type char, a message string indicating if
        %               something went wrong. If data was loaded
        %               successfully, msg will be empty
        
            
            msg  = '';
            data = [];
            
            if exist(file, 'file')
                
                fileData = load( file );
                values  = struct2cell(fileData);
                finames = fieldnames(fileData);
                
                nVals = numel(values);
                found = false;
                
                for i = 1:nVals
                    
                    vi = values{i};
                    if isa( vi, 'vsv.apps.install.AppInstallData' )
                        found = true;
                        data = vi;
                        data.FileLocation = file;
                        data.FileVarName  = finames{i};
                        break;
                    end 
                    
                end
                
                if ~found
                    msg = 'File Does not contain an app install data object';
                end
                
            else
                msg  = 'Given file does not exist';
            end
            
        end
        
        function info = findInstallationFile( fpath, launcherFile )
        % Find the installation file which can be located in one of the
        % parent folders
        %
        
            info = [];

            while ~isempty(fpath) && ~strcmp(fpath, filesep)
                filePath = fullfile( fpath, launcherFile );

                info = vsv.apps.install.AppInstallData.fromFile( filePath );
                if isempty(info)
                    
                    oldFilePath = fpath;
                    % this should give the parent folder
                    [fpath, ~, ~] = fileparts( fpath ); 
                    if strcmp(oldFilePath, fpath)
                        break;
                    end
                else
                    break; % found it
                end

            end

        end

    end
    
    methods(Static, Hidden)
        
        function mustBeValidAppDataVerbosity(verbosity)
            mustBeValidVerbosity(verbosity);
        end
        
    end
    
    %% public methods
    methods
        
        
        function this = AppInstallData( varargin )
        %APPINSTALLDATA Construct an instance of AppInstallData
        %   
        %
        % Create an instance using parameter value pairs
        %
        % @Usage
        %
        %   installData = vsv.apps.install.AppInstallData( ...
        %                     'InstallFolder',   cd , ...
        %                     'InstallVersion',  '1.0', ...
        %                     'StartFunction',   'myFunction', ...
        %                     'FunctionPackage', 'vsv.apps.launcher', ...
        %                     'InstallAppID',    'MyApp');
        
            params = varargin(1:2:end);
            vals = varargin(2:2:end);
            
            nParams = numel(params);
        
            for i = 1:nParams
                this.(params{i}) = vals{i};
            end
            
        end
        
        function [success, errMsg] = activateInstallVersion(this)
        % activates the software from which this info was instaleld.
        %
        % @return success - @type logical, true if activation was
        %                   successful, false otherwise
        % @return errMSg - @type char, an messeage specifying the error in
        %                  case success == false, otherwise errMsg == ''
        
            oldPath = cd;
            success = true;
            errMsg  = '';
            
            if exist( this.InstallFolder, 'dir' )
                cd( this.InstallFolder );
            else
                success = false;
                errMsg = [ 'Cannot access installation folder: ' newline ...
                              this.InstallFolder ];
            end
            
            if success
                
                try
                    evalin('base', 'activate();')
                catch
                    errMsg = [ 'Cannot run activate in installation folder: ' newline ...
                              this.InstallFolder ];
                    success = false;
                end
                cd( oldPath );
                
            end
            
        end
        
        
        function ver = get.InstallVersionStr(this)
        % returns the version of the installed software as a char array
        %
        % @return ver - @type char, the version as a char array
        
            ver = [ num2str(this.InstallVersion(1) ) '.' ...
                    num2str(this.InstallVersion(2) ) '.' ...
                    num2str(this.InstallVersion(3) ) ];
        end
        
        function isNewer = isNewerVersionThan(this, version)
        % check whether the underlying version is newer then the given
        % version
        %
        % @param version - @type numeric, 1, two, or three-elment version
        %                   vector
        % @return isNewer - @type logical, true if the version as defined
        %                   by this appInstallData is newer then the given
        %                   version
        
            nVersion = numel(version);
            
            switch nVersion
                case {1 2 3}
            
                    if ~this.doesVersionMatch(version)
                        thisVersion = this.InstallVersion(1:nVersion);

                        isNewerSubVersion = thisVersion > version;
                        isEqualSubVersion = thisVersion == version;
                        
                        isNewer = true;
                        for i = 1:numel(isNewerSubVersion)
                            if isNewerSubVersion(i)
                                break;
                            elseif ~isEqualSubVersion(i)
                                isNewer = false;
                                break;
                            end
                        end
                    else
                        isNewer = false;
                    end
                    
                otherwise
                    error( 'AppInstallData:doesVersionMatch:invalidVersionLenght', ...
                        'version must have length 1,2, or 3');
            end
            
            
        end
        
        function doesMatch = doesVersionMatch(this, version)
        % returns true if given version matches the underlying version
        % 
        % @param version - @type numeric array, max length 3 min length is
        %                  1, if version length is less than 3 only the
        %                  first version numbers will be compared. This
        %                  allwos to check whether major version numbers
        %                  match. E.g., if version == [4 3] this will check
        %                  whether underlying InstallVersion(1:2) ==
        %                  version. If all matches this will return true
        % @return doesMatch - @type logical, true if version matches
        
            nVersion = numel(version);
            switch nVersion
                case {1 2 3}
                    doesMatch = all( this.InstallVersion(1:nVersion) == version );
                otherwise
                    error( 'AppInstallData:doesVersionMatch:invalidVersionLenght', ...
                        'version must have length 1,2, or 3');
            end
            
        end
        
        
        function startApplication( this )
        % This will run the application as described by 
        %
        % 
            
            % see whether the package is available for the starter function
            % if this is the case we can just run that function
            % metaPackage = meta.package.fromName( this.FunctionPackage );
            
            % now check whether Vantage is activated and vsv is 
            if this.isVantageActivated()
                
                % check whether the start function is a matlab function if
                % so everything should work and just execute that function
                if isFunction( this.StartFunction )

                    % this means all is activated and ready to go
                    this.runApplicationInBase();

                else
                    % Ok this is critical, there could be the reason that
                    % the user function was not included, try that first
                    this.runUserActivation();
                    % check again whether the user function now exists and
                    % is a function , if this is the case we run the matlab
                    % application
                    if isFunction( this.StartFunction ) 
                        this.runMatlabApp();
                    else % this was not the solution we have to check the versions
                        
                        % this may mean conflicting versions or at least the
                        % function could not be found. In the rare case
                        % what could happen is that the user deleted files
                        % from the app installation folder 
                        systemVersion = getRunningSWVersion();
                        if this.isNewerVersionThan( systemVersion )
                            this.throwWarning( [ ...
                              'A vantage version was activated but it appears ' ... 
                              'that the version of the application is newer and that the activated version does not contain the necessary code.' ...
                              newline newline 'restart matlab and then start the application with out activating Vantage']);
                        else % this means the system version is newer or equal to the version 
                             % of the installed application this is a
                             % problem. In this case this is a real failure
                             % and should require to reinstall the software
                             % this can actually only happen if files from
                             % the installation folder were removed
                             this.IsValid = false;
                             this.updateInstallDataFile();
                             this.throwWarning( [ ...
                                  'A vantage version was activated but it appears ' ... 
                                  'that it does not match the version of the application.' ...
                                  newline newline ...
                                  ]);
                        end
                    end
                end
                
            else
                
                [success, errMsg] = this.activateInstallVersion();
                if success
                    % run the matlab app to reactivate and add the path for
                    % the application!!! to not run the start function
                    this.runMatlabApp();
                else
                    this.IsValid = false;
                    this.updateInstallDataFile();
                    this.throwWarning([ 'Cannot not activate Installed version: ' ...
                              newline newline ...
                              errMsg] );
                end
                
            end
        end
        
        
    end
    
    % provide access for testing
    methods(Access=private)
        
        function updateInstallDataFile(this)
            if ~isempty(this.FileLocation) && ~isempty(this.FileVarName)
                
                eval([ this.FileVarName ' = this;']);
                if exist( this.FileLocation, 'file')
                    delete(this.FileLocation);
                end
                
                save( this.FileLocation, this.FileVarName );
            end
        end
    end
    
    methods(Access=private)
        
        
        
        function activateUserAndRunApp(this, errMsg  )
            this.runUserActivation();
            if isFunction( this.StartFunction ) 
                this.runMatlabApp();
            else
                this.IsValid = false;
                this.updateInstallDataFile();
                this.throwNoUserFunctionError( errMsg );
            end    
        end
        
        function throwNoUserFunctionError(this, errMsg )
            this.throwWarning( [ ...
                      'After activating Vantage and user defined directories, ' newline ...
                      'still cannot find Function: ' this.StartFunction ...
                      newline newline ...
                      errMsg] );
        end
        
        function runUserActivation(this)
            cd( this.InstallFolder );
            this.addUserPath(  );
            this.runUserActivate(  );
        end
        
        function isactive = isVantageActivated(~)
        % returns true if a vantage version was activated already
        %
        % @return isactive - @type logical true if a vantage version was
        %                    activated already
            isactive = getenv('VERASONICS_IS_ACTIVATED') == '1';
        end
        
        function success = runUserActivate(this)
            
            userActivateFnc = this.UserActivate;
            success = true;
            if ~isempty(userActivateFnc)
                this.executeFunctionSecure(userActivateFnc, 'The user defined activate function');
            end
            
        end
        
        function executeFunctionSecure(this, fncName, fncDesc)
        % execute the function in a secure manner and throw a warning in
        % case the function does not work
        %
        %   There are two possibilities, first the function does not
        %   exists. We report that if it does not exist. In case it does
        %   exist there can be the possibility that the functions thorows
        %   an error
        %
        % @param fncName - @type char the name of the function
        % @param fncDesc - @type char a description of the function used
        %                  in the warning and error messages. E.g., fncDesc
        %                  = 'The user defined activate function'

            try 
                if isFunction( fncName )
                    evalin('base', [fncName ';'] );
                else
                    this.throwWarning( ...
                            [ fncDesc ' could not be found:'  ...
                            newline 'Function: ' fncName] );
                end
            catch errID
                report = getReport(errID);
                this.IsValid = false;
                this.throwWarning( ...
                    [ fncDesc ' threw the following error:'  newline ...
                    '' newline ...
                    report ] );
            end
        end
        
        function success = addUserPath(this)
        % This will add the pathes as defined by the user
        %

            userPath = this.UserPath;
            
            nPath = numel(userPath);
            allExists = true;
            success = false;
            
            for i = 1:nPath
                if ~exist( userPath{i}, 'dir' )
                    allExists = false;
                    break
                end
            end
            
            if ~allExists
                this.IsValid = false;
                this.throwWarning( ...
                       ['One or more of the user defined pathes does not exist: ' newline ...
                       userPath{i}]  );
            else
                
                id = 'MATLAB:mpath:nameNonexistentOrNotADirectory';

                status = warning( 'query', id);
                warning( 'off', id);

                if ~isempty(userPath)
                    % cell expansion
                    addpath( userPath{:} );
                end

                % restore warnign state
                warning( status.state, id);
                success = true;
            end
            
        end
        
        function throwWarning(this, message )
        % will throw the given warning message based on the defined
        % verboisty level 
        
            verbosity = this.Verbosity;
            
            nVerbose = numel(verbosity);
            
            for i = 1:nVerbose
                verbose = verbosity{i};
                switch verbose
                    case "nothing"

                    case "print"
                        disp(message);
                    case "warning"
                        warning( 'AppInstallData:warning', message );
                    case "warndlg"
                        warndlg( message );
                    case "save"
                        this.LastError = message;
                        this.updateInstallDataFile();
                end
            end
            
        end
        
        function runApplicationInBase(this)
            
            startFnc = this.StartFunction;
            
            if ~isempty(startFnc)
                this.executeFunctionSecure(startFnc, 'The single entry application start function');
            else
                this.IsValid = false;
                this.throwWarning( 'The single entry application start function is empty' );
            end
            
        end
        
        function runMatlabApp(this)
        % this will run the app from the app menu. This will put the path back 
        % to the matlab path if lost for some reason. 
        %
        % This should be run after activate instead of runApplication.
        % Because this will then ensure that the correct pathes are set.
        %
            ID = this.InstallAppID;
            
            try
                matlab.apputil.run( ID );
            catch e
                errmsg = [ 'No app found for ID: ' ID, newline ...
                            e.message ];
                this.IsValid = false;
                this.updateInstallDataFile();
                this.throwWarning( errmsg);
            end
            
        end
        
    end
    
end

%% Private non-class member helper functions

function version = getRunningSWVersion()
% Returns the version of the current system. This uses the version
% of runAcq
%
% @return version - @type numeric 

    cmmdStr = 'vsv.vsx.command.getSWVersion();';
    try
        %  we have to make this eval to not include this into the exported
        %  code!
        version = eval( cmmdStr );
    catch
        version = [];
    end
            
end
        
function mustBeCellString(cell)
%MUSTBECELLSTRING will throw an error if cell is not a cell string
%   A cell string is a cell array that only contains char arrays
%
% Version 1.0 | 2020-05-15 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    if ~iscellstr(cell) %#ok 
        error('vsv:mustBeCellString', 'Given cell array is not a cellstring')
    end
    
end

function mustBeValidVerbosity(verbose)

    mustBeCellString(verbose);
    
    nVerbose = numel(verbose);
    for i = 1:nVerbose
        ele = verbose{i};
        switch ele
            case { 'nothing', 'print', 'warning', 'warndlg', 'save' }
                isOK = true;
            otherwise
                isOK = false;
                break;
        end
    end
    
    if ~isOK
        error('vsv:mustBeValidVerbosity', ...
            [ 'All verbosity values must be a one of the ' newline ... 
              ' { ''nothing'', ''print'', ''warning'', ''warndlg'' }' ] );
    end
    
end

function isf = isFunction(fnc)

    try
        nargin(fnc); % nargin errors when FUNNAME is not a function
        isf = true;
    catch 
        isf = false;
    end
end
