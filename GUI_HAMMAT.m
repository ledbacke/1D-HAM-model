function varargout = GUI_HAMMAT(varargin)
% GUI_HAMMAT MATLAB code for GUI_HAMMAT.fig
%      GUI_HAMMAT, by itself, creates a new GUI_HAMMAT or raises the existing
%      singleton*.
%
%      H = GUI_HAMMAT returns the handle to a new GUI_HAMMAT or the handle to
%      the existing singleton*.
%
%      GUI_HAMMAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_HAMMAT.M with the given input arguments.
%
%      GUI_HAMMAT('Property','Value',...) creates a new GUI_HAMMAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_HAMMAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_HAMMAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_HAMMAT

% Last Modified by GUIDE v2.5 16-Mar-2016 15:16:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_HAMMAT_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_HAMMAT_OutputFcn, ...
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


% --- Executes just before GUI_HAMMAT is made visible.
function GUI_HAMMAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_HAMMAT (see VARARGIN)

handles.Height = 0.12;
set(handles.Edit_Height,'string',handles.Height);
handles.Width  = 0.12;
set(handles.Edit_Width,'string',handles.Width);
handles.NumberLayers = 1;
set(handles.Edit_NumberLayers,'string',handles.NumberLayers);
handles.NumberLayersPrevious = 1;
handles.ArrayToMexGeom = zeros(1,3);

handles.ThicknessLayer = [0.01];
set(handles.uitableThicknessLayer,'Data',handles.ThicknessLayer);
handles.Exponent = 0;
set(handles.edit_Exponent,'string',handles.Exponent);
handles.Mesh = [10];
set(handles.uitable_Mesh_Cell_Layer,'Data',handles.Mesh);

% MENU MATERIAL DATA
handles.DryDensity  = [550];
set(handles.uitable_DryDensity,'Data',handles.DryDensity);
handles.ThermalChar = [1500,0.11,0];
set(handles.uitable_ThermalChar,'Data',handles.ThermalChar);
handles.HygricChar  = [-0.043387,0.050052,-0.000048,3.59e-12];
set(handles.uitable_HygricChar,'Data',handles.HygricChar);

handles.InitialConditions = [23,60];
set(handles.uitable_InitialConditions,'Data',handles.InitialConditions);
handles.ConditionsLeft  = [23;93;101325;8;0];
set(handles.uitable_ConditionsLeft,'Data',handles.ConditionsLeft);
handles.ConditionsRight = [23;93;101325;8;0];
set(handles.uitable_ConditionsRight,'Data',handles.ConditionsRight);

handles.StartTime = 0;
handles.StopTime  = 60*60*3600*365;
handles.TimeStep  = 15*60;
handles.ArrayToMexParametersTime = zeros(1,3);

handles.EstimateDiffusionCoefficient = 0;
set(handles.checkbox_EstimateDiffusionCoefficient,'Value',handles.EstimateDiffusionCoefficient);
set(handles.textWhichLayer,'ForegroundColor',[0.8,0.8,0.8]);
set(handles.editWhichLayer,'Enable','off'); 
handles.EstimateDiffCoeff = -1;
handles.MeasuredWeight = zeros(2,1);

% Choose default command line output for GUI_HAMMAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_HAMMAT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_HAMMAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%MENU GEOMETRY
%-------------------------------------------------------------------------------------
function Panel_Geometry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Panel_Geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
function Edit_Height_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_Height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_Height as text
%        str2double(get(hObject,'String')) returns contents of Edit_Height as a double

Height = str2double(get(hObject,'String'));
if isnan(Height)|| Height <= 0
   set(handles.Edit_Height,'string', 0.12);
end
handles.Height = Height;
guidata(hObject,handles);
function Edit_Height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_Height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_Width_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_Width as text
%        str2double(get(hObject,'String')) returns contents of Edit_Width as a double
Width= str2double(get(hObject,'String'));
if isnan(Width) || Width <= 0
    Width = 0.1198;
    set(handles.Edit_Width,'string', Width);
end
handles.Width = Width;
guidata(hObject,handles);
function Edit_Width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_NumberLayers_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_NumberLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_NumberLayers as text
%        str2double(get(hObject,'String')) returns contents of Edit_NumberLayers as a double

% Set Default Value
 
% Set value if wrong input
NumberLayers = str2num(get(hObject,'String'));
if isnan(NumberLayers)|| NumberLayers <= 0
    NumberLayers = 1;
    set(handles.Edit_NumberLayers,'string', NumberLayers);
end
NumberLayers = str2num(get(hObject, 'String'));

%Store Value
handles.NumberLayers = NumberLayers;
ExtraRows = NumberLayers-handles.NumberLayersPrevious;
handles.NumberLayersPrevious =NumberLayers; 
guidata(hObject,handles);

%Make extra line in Tables without losing data
Old_data_Table_ThicknessLayer   = get(handles.uitableThicknessLayer,'Data');
Old_data_Table_NumberGridLayers = get(handles.uitable_Mesh_Cell_Layer,'Data');
Old_data_Table_DryDensity       = get(handles.uitable_DryDensity,'Data');
Old_data_Table_ThermalChar      = get(handles.uitable_ThermalChar,'Data');   
Old_data_Table_HygricChar       = get(handles.uitable_HygricChar,'Data');
Old_data_Table_InitialConditions= get(handles.uitable_InitialConditions,'Data');

if NumberLayers > size(Old_data_Table_ThicknessLayer,1); 
    data_Table_ThicknessLayer   = [Old_data_Table_ThicknessLayer; zeros(ExtraRows,1)];
    data_Table_NumberGridLayers = [Old_data_Table_NumberGridLayers; zeros(ExtraRows,1)];
    data_Table_DryDensity       = [Old_data_Table_DryDensity; zeros(ExtraRows,1)];
    data_Table_ThermalChar      = [Old_data_Table_ThermalChar; zeros(ExtraRows,3)];
    data_Table_HygricChar       = [Old_data_Table_HygricChar; zeros(ExtraRows,4)];
    data_Table_InitialConditions= [Old_data_Table_InitialConditions; zeros(ExtraRows,2)];
    
elseif NumberLayers < size(Old_data_Table_ThicknessLayer,1)
    data_Table_ThicknessLayer   = Old_data_Table_ThicknessLayer(1:NumberLayers,:);
    data_Table_NumberGridLayers = Old_data_Table_NumberGridLayers(1:NumberLayers,:);
    data_Table_DryDensity       = Old_data_Table_DryDensity(1:NumberLayers,:);
    data_Table_ThermalChar      = Old_data_Table_ThermalChar(1:NumberLayers,1:3);
    data_Table_HygricChar       = Old_data_Table_HygricChar(1:NumberLayers,1:4);
    data_Table_InitialConditions= Old_data_Table_InitialConditions(1:NumberLayers,1:2);
end
set(handles.uitableThicknessLayer,'Data',data_Table_ThicknessLayer);
set(handles.uitable_Mesh_Cell_Layer,'Data',data_Table_NumberGridLayers);
set(handles.uitable_DryDensity,'Data',data_Table_DryDensity);
set(handles.uitable_ThermalChar,'Data',data_Table_ThermalChar);
set(handles.uitable_HygricChar,'Data',data_Table_HygricChar);
set(handles.uitable_InitialConditions,'Data',data_Table_InitialConditions);
guidata(hObject, handles);
function Edit_NumberLayers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_NumberLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Thickness of each material layer
function uitableThicknessLayer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitablethicknesslayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', 0.0122)
set(hObject, 'ColumnName',{'Thickness layer [m]'})
function uitableThicknessLayer_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableThicknessLayer (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.ThicknessLayer = get(handles.uitableThicknessLayer,'Data');
guidata(hObject, handles);

% MENU MESH
%-----------------------------------------------------------------------------------------------------
%Number of mesh cells for each material layer
function uipanel_Mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Panel_Geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function uitable_Mesh_Cell_Layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_Mesh_Cell_Layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1,1))
set(hObject, 'ColumnName',{'Number of cells for layer '})
function uitable_Mesh_Cell_Layer_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_Mesh_Cell_Layer (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.Mesh = get(handles.uitable_Mesh_Cell_Layer,'Data');
guidata(hObject, handles);

%GridExponent
function edit_Exponent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Exponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Exponent as text
%        str2double(get(hObject,'String')) returns contents of edit_Exponent as a double

Exponent_GridConcentration = str2double(get(hObject, 'String'));
handles.Exponent = Exponent_GridConcentration;
guidata(hObject,handles);
function edit_Exponent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Exponent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%pushbutton Plot Mesh
function pushbutton_PlotMesh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlotMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Calculate Mesh with function from C++
mex Mesh.cpp;
[Length,DeltaX] = Mesh(handles.Mesh,handles.ThicknessLayer,handles.Exponent);


%fh = figure();
varargout = PlotMesh();



% MENU MATERIAL DATA
%-----------------------------------------------------------------------------------------------------
%Dry Density [kg/m�]
function uitable_DryDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_DryDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1,1))
set(hObject, 'ColumnName',{'Dry density of layers '})
function uitable_DryDensity_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_DryDensity (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.DryDensity = get(handles.uitable_DryDensity, 'Data');
guidata(hObject, handles);

%Thermal Char
function uitable_ThermalChar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_ThermalChar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1,3))
set(hObject, 'ColumnName',{'Heat Capacity [J/(kg.K)]','Lamda dry [W/(m.K)]','Lamda wet [%]'})
function uitable_ThermalChar_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_DryDensity (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.ThermalChar = get(handles.uitable_ThermalChar,'Data');
guidata(hObject, handles);

%Hygric Char
function uitable_HygricChar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_HygricChar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1,4));
set(hObject, 'ColumnName',{'A','B','C','Diff Coeff [kg/(m.s.Pa)]'});
function uitable_HygricChar_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_DryDensity (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.HygricChar = get(handles.uitable_HygricChar,'Data');
guidata(hObject, handles);


% MENU INITIAL AND BOUNDARY CONDITIONS
%-----------------------------------------------------------------------------------------------------

%Initial Conditions
function uitable_InitialConditions_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_InitialConditions (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.InitialConditions = get(handles.uitable_InitialConditions, 'Data');
guidata(hObject, handles);

function uitable_InitialConditions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_InitialConditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(1,2));
set(hObject, 'ColumnName',{'T[�C]','RH[%]'});


%Outside Conditions
function checkbox_OutsideLeft_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_OutsideLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_OutsideLeft

function checkbox_OutsideRight_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_OutsideRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_OutsideRight


% Conditions Left
function uitable_ConditionsLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_ConditionsLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(5,1))
set(hObject, 'RowName',{'T[�C]','RH [%]','Atm.pres [Pa]','Alpha[W/m�K]','Beta[kg/m�sPa]'},'ColumnName',{'Value'})
% --- Executes when entered data in editable cell(s) in uitable_ConditionsLeft.
function uitable_ConditionsLeft_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_ConditionsLeft (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.ConditionsLeft = get(handles.uitable_ConditionsLeft,'Data');
guidata(hObject, handles);

%Conditions Right
function uitable_ConditionsRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_ConditionsRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(5,1))
set(hObject, 'RowName',{'T[�C]','RH [%]','Atm.pres [Pa]','Alpha[W/m�K]','Beta[kg/m�sPa]'},'ColumnName',{'Value'})
% --- Executes when entered data in editable cell(s) in uitable_ConditionsLeft.
function uitable_ConditionsRight_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_ConditionsLeft (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.ConditionsRight = get(handles.uitable_ConditionsRight,'Data');
guidata(hObject, handles);



% MENU NUMERICAL CALCULATION
%-----------------------------------------------------------------------------------------------------

%Starttime
function edit_StartTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StartTime as text
%        str2double(get(hObject,'String')) returns contents of edit_StartTime as a double
StartTime = str2double(get(hObject,'String'));
handles.StartTime = StartTime;
guidata(hObject, handles);
function edit_StartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Stoptime
function edit_StopTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StopTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StopTime as text
%        str2double(get(hObject,'String')) returns contents of edit_StopTime as a double
StopTime = str2double(get(hObject,'String'));
handles.StopTime = StopTime;
guidata(hObject, handles);
function edit_StopTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StopTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%TimeStep
function edit_TimeStep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TimeStep as text
%        str2double(get(hObject,'String')) returns contents of edit_TimeStep as a double
TimeStep = str2double(get(hObject,'String'));
handles.TimeStep = TimeStep;
guidata(hObject, handles);
function edit_TimeStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Store data in Matrix Input_NumCal


% MENU OPTIONS CALCULATION
%-----------------------------------------------------------------------------------------------------
% --- Executes on button press in checkbox_EstimateDiffusionCoefficient.
function checkbox_EstimateDiffusionCoefficient_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_EstimateDiffusionCoefficient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.EstimateDiffusionCoefficient = get(hObject,'Value');
if handles.EstimateDiffusionCoefficient == 1
    [FileName, PathName] = uigetfile({'*.xls';'*.xlsx'},'Measured Weight','MultiSelect', 'off');
    set(handles.textWhichLayer,'ForegroundColor',[0,0,0]);
    set(handles.editWhichLayer,'Enable','on');
    if isequal(FileName,0) | isequal(PathName,0)
        disp('User selected Cancel')
        disp('Estimate Diffusion Coefficient not possible')
        set(handles.checkbox_EstimateDiffusionCoefficient,'Value',0.0);
        set(handles.textWhichLayer,'ForegroundColor',[0.8,0.8,0.8]);
        set(handles.editWhichLayer,'Enable','off'); 
        handles.EstimateDiffCoeff = -1;
        guidata(hObject, handles);
    else
        %fh = figure();
         Name = fullfile(PathName,FileName);
         set(handles.editWhichLayer,'Enable','on');
         set(handles.textWhichLayer,'ForegroundColor',[0,0,0]);
         varargout = MeasuredWeight(Name);
         Measurement = xlsread(FileName);
         set(handles.MeasuredWeight,'Data',Measurement);
         guidata(hObject, handles);
    end
else
    set(handles.textWhichLayer,'ForegroundColor',[0.8,0.8,0.8]);
    set(handles.editWhichLayer,'Enable','off');
    handles.EstimateDiffCoeff = handles.editWhichLayer;
    guidata(hObject, handles);
end

%Define Layer for which to estimate diffusion coefficient
function textWhichLayer_CreateFcn(hObject, eventdata, ~)
% hObject    handle to textWhichLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
function editWhichLayer_Callback(hObject, eventdata, handles)
% hObject    handle to editWhichLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.editWhichLayer = str2double(get(hObject,'String'));
handles.EstimateDiffCoeff = handles.editWhichLayer;
function editWhichLayer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWhichLayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% Pushbutton Start Calculation.
function pushbuttonStartCalculation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStartCalculation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ArrayToMexGeom = [handles.Height,handles.Width,handles.NumberLayers];
handles.ArrayToMexParametersTime = [handles.StartTime,handles.StopTime,handles.TimeStep];
[CalculatedWeight,CalculatedDiffCoeff] =     HAM_MATLAB(handles.ArrayToMexGeom,...
                                                        handles.ThicknessLayer,...
                                                        handles.Exponent,...    
                                                        handles.Mesh,...
                                                        handles.DryDensity,...
                                                        handles.ThermalChar,...
                                                        handles.HygricChar,...                                                        
                                                        handles.InitialConditions,...                                                        
                                                        handles.ConditionsLeft,...
                                                        handles.ConditionsRight,...
                                                        handles.ArrayToMexParametersTime,...
                                                        handles.EstimateDiffCoeff,...
                                                        handles.MeasuredWeight);
