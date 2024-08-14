function vsCustomGUI()

    delete(findobj('tag','UI'));

    ScrnSize = get(0,'ScreenSize');
    Bkgrnd = [0.8 0.8 0.8];
    
    figure( 'Visible','on',...  %'Units','normalized',...
            'Position',[ScrnSize(3)-500,(ScrnSize(4)-620)/2,450,620],... %'Position',[0.7,0.25,0.25,0.50],...
            'Name','VSX Control',...
            'Color',Bkgrnd,...
            'NumberTitle','off',...
            'MenuBar','none', ...
            'Resize','on', ...
            'DefaultUicontrolBackgroundColor',Bkgrnd,...
            'CloseRequestFcn', @closefunc, ...
            'tag','UI');
    
    setupTPC();


end

function closefunc(source,varargin)

    assignin('base', 'vsExit', 1);

    delete(source);
    hFigure = findobj('tag','UI');
    if ~isempty(hFigure)
        close(hFigure(contains(get(hFigure,'Tag'),'Msgbox')));
    end

end

function setupTPC()
% High Voltage slider.
%   The slider's min and max range is determined by the limits of the Verasonics
%   TPC and a max voltage limit that can be set in the user's setup script.
%
%   The user script may override the "Max" attribute of this slider to impose
%   a high voltage maximum limit that is less than the capability of the detected
%   Verasonics TPC.

    TPC = evalin('base','TPC');
    hvps2 = 0; 
    % - For more than one active profile, determine the number to use for the 2nd slider. Profile 5
    %   has priority over profiles 2-4.
    if TPC(5).inUse
        hvps2 = 5;
    else
        for i = 2:4
            if TPC(i).inUse ~= 0
                hvps2 = i;
                break,
            end
        end
    end
    
    %VSX requires this parameter
    assignin('base', 'hv2GUIprofile', hvps2);
end