function varargout = MultiAxesMeasureGUI(varargin)
% MULTIAXESMEASUREGUI MATLAB code for MultiAxesMeasureGUI.fig
%      MULTIAXESMEASUREGUI, by itself, creates a new MULTIAXESMEASUREGUI or raises the existing
%      singleton*.
%
%      H = MULTIAXESMEASUREGUI returns the handle to a new MULTIAXESMEASUREGUI or the handle to
%      the existing singleton*.
%
%      MULTIAXESMEASUREGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIAXESMEASUREGUI.M with the given input arguments.
%
%      MULTIAXESMEASUREGUI('Property','Value',...) creates a new MULTIAXESMEASUREGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MultiAxesMeasureGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MultiAxesMeasureGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MultiAxesMeasureGUI

% Last Modified by GUIDE v2.5 22-Jun-2021 07:17:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MultiAxesMeasureGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MultiAxesMeasureGUI_OutputFcn, ...
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


% --- Executes just before MultiAxesMeasureGUI is made visible.
function MultiAxesMeasureGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MultiAxesMeasureGUI (see VARARGIN)

% Choose default command line output for MultiAxesMeasureGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MultiAxesMeasureGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MultiAxesMeasureGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in PopupDisplay.
function PopupDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to PopupDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopupDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopupDisplay


% --- Executes during object creation, after setting all properties.
function PopupDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopupDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PopupMeasureType.
function PopupMeasureType_Callback(hObject, eventdata, handles)
% hObject    handle to PopupMeasureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopupMeasureType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopupMeasureType


% --- Executes during object creation, after setting all properties.
function PopupMeasureType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopupMeasureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonMeasure.
function ButtonMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
