function varargout = parameters(varargin)
% PARAMETERS MATLAB code for parameters.fig
%      PARAMETERS, by itself, creates a new PARAMETERS or raises the existing
%      singleton*.
%
%      H = PARAMETERS returns the handle to a new PARAMETERS or the handle to
%      the existing singleton*.
%
%      PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERS.M with the given input arguments.
%
%      PARAMETERS('Property','Value',...) creates a new PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parameters

% Last Modified by GUIDE v2.5 30-Sep-2016 20:06:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @parameters_OutputFcn, ...
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

% --- Executes just before parameters is made visible.
function parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parameters (see VARARGIN)

% Choose default command line output for parameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using parameters.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = parameters_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    load test.mat
    I0=read(V,frame); % ԭʼ��Ƶ֡ͼ��
    I0=rgb2gray(I0);
    I0=imresize(I0,resize); %����   
    % �Ҷ�ģʽ����
    if strcmp(darkfish,'yes')
        I = back - I0; %�뱳������Ϊ��������������back-I����֮I-back
    else
        I = I0 - back;
    end
    I=I-region; 
   imshow(I0);
   I(I<minLight) = 0; % 2015.9.8 ������С��minLight�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺ�
    level=graythresh(I);
    I=im2bw(I,level);
    
    if strelSize >= 1
        se=strel('square',strelSize); % strelSizeӦ�ø������С������һ�����Ϊ3��4
        %http://www.cnblogs.com/tornadomeet/archive/2012/03/20/2408086.html
        % 2015.6.11 ������imclose������imopen
%         I=imo   pen(I,se); % �ȸ�ʴ�����ͣ�����С���塢����ϸ���������塢ƽ��������߽磬�˴���Ҫ�Ǹɵ���β��
        I=imclose(I,se); % �����ͺ�ʴ������ڲ��ն��������ٽ����塢ƽ��������߽磬�˴�������
    end
    I=bwareaopen(I,minS,4); %��4�ھӲ�����ͨ���ţ�ȥ�����ص�����minS����
    
    % ����Եõ�������ͨ���Ž�����ʽ�ָ�
    [L,num] = bwlabel(I, 4); %����ͨ����
    S0 = regionprops(L,'Area'); %����ͨ�������
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %����ͨ����λ��
    C0 = cat(1,C0.Centroid);
    len0 = regionprops(L,'MajorAxisLength'); % ��Բ����
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % ��Բ����
    wid0 = cat(1,wid0.MinorAxisLength);
    Ec0 = regionprops(L,'Eccentricity'); % ������
    Ec0 = cat(1,Ec0.Eccentricity);
    B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
    B0=cat(1,B0.BoundingBox);
    global numfish;
    numfish=num;
set(handles.edit3,'string',num2str(numfish));
set(handles.edit1,'string',num2str(minS));
set(handles.edit2,'string',num2str(maxS));
set(handles.edit5,'string',num2str(max(S0)));
set(handles.edit6,'string',num2str(min(S0)));
set(handles.edit7,'string',num2str(mean(S0)));
set(handles.edit8,'string',num2str(frame));

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global minArea;
minArea=get(handles.edit1,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global maxArea;
maxArea=get(handles.edit2,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    load test.mat
    global fframe;
    frame1=fframe;
    I0=read(V,frame1); % ԭʼ��Ƶ֡ͼ��
    I0=rgb2gray(I0);
    I0=imresize(I0,resize); %����   
    % �Ҷ�ģʽ����
    if strcmp(darkfish,'yes')
        I = back - I0; %�뱳������Ϊ��������������back-I����֮I-back
    else
        I = I0 - back;
    end
    I=I-region; 
    imshow(I);
    I(I<minLight) = 0; % 2015.9.8 ������С��minLight�ĵ��Ϊȫ�ڣ�֮������graythresh������Ч���ܺ�
    level=graythresh(I);
    I=im2bw(I,level);
    global maxArea;
    global minArea;
    minS=str2double(minArea);
    maxS=str2double(maxArea);
    if strelSize >= 1
        se=strel('square',strelSize); % strelSizeӦ�ø������С������һ�����Ϊ3��4
        %http://www.cnblogs.com/tornadomeet/archive/2012/03/20/2408086.html
        % 2015.6.11 ������imclose������imopen
%         I=imopen(I,se); % �ȸ�ʴ�����ͣ�����С���塢����ϸ���������塢ƽ��������߽磬�˴���Ҫ�Ǹɵ���β��
        I=imclose(I,se); % �����ͺ�ʴ������ڲ��ն��������ٽ����塢ƽ��������߽磬�˴�������
    end
    I=bwareaopen(I,minS,4); %��4�ھӲ�����ͨ���ţ�ȥ�����ص�����minS����
    
    % ����Եõ�������ͨ���Ž�����ʽ�ָ�
    [L,num] = bwlabel(I, 4); %����ͨ����
    S0 = regionprops(L,'Area'); %����ͨ�������
    S0 = cat(1,S0.Area);
    C0 = regionprops(L,'Centroid'); %����ͨ����λ��
    C0 = cat(1,C0.Centroid);
    len0 = regionprops(L,'MajorAxisLength'); % ��Բ����
    len0 = cat(1,len0.MajorAxisLength);
    wid0 = regionprops(L,'MinorAxisLength'); % ��Բ����
    wid0 = cat(1,wid0.MinorAxisLength);
    Ec0 = regionprops(L,'Eccentricity'); % ������
    Ec0 = cat(1,Ec0.Eccentricity);
    B0=regionprops(L,'BoundingBox'); % ����ͨ����߽�
    B0=cat(1,B0.BoundingBox);
global numfish;
numfish=num;
set(handles.edit3,'string',num2str(numfish));
set(handles.edit1,'string',num2str(minS));
set(handles.edit2,'string',num2str(maxS));
set(handles.edit5,'string',num2str(max(S0)));
set(handles.edit6,'string',num2str(min(S0)));
set(handles.edit7,'string',num2str(mean(S0)));
set(handles.edit8,'string',num2str(frame1));
imshow(I)
save('para.mat','minS','maxS');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con=1;
save('continue.mat','con');
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global fframe;
fframe=str2double(get(handles.edit8,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
