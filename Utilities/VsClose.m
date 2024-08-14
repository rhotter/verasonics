function VsClose
% close Matlab figures and Verasonics Viewer
close all

% First ask the Java Garbage Collector to clean up.
java.lang.System.gc();

% Get the list of windows that existed at some point in time.
VsWindow = java.awt.Window.getWindows();

for n = 1:length(VsWindow)
    if ~isempty(strfind(VsWindow(n),'VantageWindow'))
        VsWindow(n).dispose();
        VsWindow(n) = [];
    end
end

end
