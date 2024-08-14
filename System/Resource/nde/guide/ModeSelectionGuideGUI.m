function varargout = ModeSelectionGuideGUI(varargin)
% MODESELECTIONGUIDEGUI MATLAB code for ModeSelectionGuideGUI.fig
%      MODESELECTIONGUIDEGUI, by itself, creates a new MODESELECTIONGUIDEGUI or raises the existing
%      singleton*.
%
%      H = MODESELECTIONGUIDEGUI returns the handle to a new MODESELECTIONGUIDEGUI or the handle to
%      the existing singleton*.
%
%      MODESELECTIONGUIDEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODESELECTIONGUIDEGUI.M with the given input arguments.
%
%      MODESELECTIONGUIDEGUI('Property','Value',...) creates a new MODESELECTIONGUIDEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModeSelectionGuideGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModeSelectionGuideGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModeSelectionGuideGUI

% Last Modified by GUIDE v2.5 05-May-2021 18:12:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModeSelectionGuideGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ModeSelectionGuideGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ModeSelectionGuideGUI is made visible.
function ModeSelectionGuideGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModeSelectionGuideGUI (see VARARGIN)

% Choose default command line output for ModeSelectionGuideGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ModeSelectionGuideGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ModeSelectionGuideGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ButtonAddMode.
function ButtonAddMode_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonAddMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonRemoveCustomMode.
function ButtonRemoveTab_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemoveCustomMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in ModeSwitchDropDown.
function ModeSwitchDropDown_Callback(hObject, eventdata, handles)
% hObject    handle to ModeSwitchDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ModeSwitchDropDown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ModeSwitchDropDown


% --- Executes during object creation, after setting all properties.
function ModeSwitchDropDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModeSwitchDropDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonNewCustomMode.
function ButtonNewCustomMode_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonNewCustomMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonRemoveCustomMode.
function ButtonRemoveCustomMode_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemoveCustomMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonOK.
function ButtonOK_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonCancel.
function ButtonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CheckCustomProcessing.
function CheckCustomProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to CheckCustomProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckCustomProcessing



function EditFunctionName_Callback(hObject, eventdata, handles)
% hObject    handle to EditFunctionName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditFunctionName as text
%        str2double(get(hObject,'String')) returns contents of EditFunctionName as a double


% --- Executes during object creation, after setting all properties.
function EditFunctionName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFunctionName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonImport.
function ButtonImport_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in BNuttonTest.
function BNuttonTest_Callback(hObject, eventdata, handles)
% hObject    handle to BNuttonTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonEdit.
function ButtonEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonRemoveParam.
function ButtonRemoveParam_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemoveParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
