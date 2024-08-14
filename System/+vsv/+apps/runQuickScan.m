function WsData = runQuickScan( )
% RUNQUICKSCAN will start the QuickScan application
%
% The QuickScan app is a graphical user interface, which allows you to
% quickly perform ultrasound scans in a more clinical-like setting. The app
% provides many features, which you would usually find on a clinical device
% including data capturing, review, and measurement tools. 
%
% You can find more details in the accompanying app note. 
%
% The function will start a script selector first, which allows you to
% select a High Image Quality script that fits your current
% hardware configuration (i.e., connected Transducer). If no hardware is
% connected  you will be able to select from all possible scripts to run in
% simulation.
%
% The method will preserve your current workspace. There is no need for
% calling 'clear all' or saving your current workspace variables. Even if
% the function will crash for some unexpected reason, your variables should
% be restored when the function exits. 
%
%
%
% @Usage:
%   % start the application
%   wsData = vsv.apps.runQuickScan();
%
%   % the returned object will allow you to access the workspace variables
%   % that were produced by the last running scrip. This data is equivalent
%   % to what would have been placed in the main workspace using the
%   % standard research GUI.
%   dataStruct = wsData.convertToStruct();
%
%   % if you would prefer to copy all the data to the main workspace you
%   % can use the following function. However, please not that this will
%   % delete all variables that are currently in the workspace!
%   wsData.restoreWorkspace();
%
% @return wsData - @type vsv.vsx.WorkspaceData, stored workspace data as
%                  produced by the script. 
%
%
% @Trouble shooting
%   If the app reports an error for any reason it will popup a dialog.
%   There can be several reasons for this. One reason could be a corrupt
%   mat file (e.g., it got accidently overwritten with invalid data). In
%   this case you can try to recompute the Matfile using the dialog. 
%
%   If this still doesn't work, please try to use clear all in the main
%   workspace and run the function again or try to restart Matlab. 
%
%   In case the error still persists we kindly ask you to write down the
%   steps that you performed that led to the error message and send a
%   bug report to the Verasonics customer service. Please also include
%   error messages and screenshots if that is possible. 
%   
%
% Version 1.0 | 2020-07-07 
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
    
    vsxApplication = vsv.apps.VsQuickScanApp;
    WsData = vsxApplication.runApplication( );
    
    if nargout < 1
        clear WsData
    end
    
end

