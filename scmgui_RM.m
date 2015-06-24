function varargout = scmgui_RM(varargin)
% this version is modified slighly to accomodate name changes to VRsets
% SCMGUI_RM MATLAB code for scmgui_RM.fig
%      SCMGUI_RM, by itself, creates a new SCMGUI_RM or raises the existing
%      singleton*.
%
%      H = SCMGUI_RM returns the handle to a new SCMGUI_RM or the handle to
%      the existing singleton*.
%
%      SCMGUI_RM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCMGUI_RM.M with the given input arguments.
%
%      SCMGUI_RM('Property','Value',...) creates a new SCMGUI_RM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scmgui_RM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scmgui_RM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scmgui_RM

% Last Modified by GUIDE v2.5 02-Mar-2015 14:15:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @scmgui_RM_OpeningFcn, ...
    'gui_OutputFcn',  @scmgui_RM_OutputFcn, ...
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


% --- Executes just before scmgui_RM is made visible.
function scmgui_RM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scmgui_RM (see VARARGIN)

% Choose default command line output for scmgui_RM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scmgui_RM wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = scmgui_RM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% get Image Sets (Load Image Set button)
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%gets folder containing set images
folder_name = uigetdir('C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\',...
    'Select Image folder');

%if file selection is cancelled, pathname should be zero and nothing should happen
if folder_name == 0
    return
end

userdata.folder_name = folder_name; %save to handle so can grab from other objects

%---get the name of images in the foldler---%
image_list = ls([folder_name '\','*.bmp']);
order = NaN(2,36);
for il = 1:size(image_list,1);
    num  = sscanf(image_list(il,6:end), '%d');
    %     if num < 10
    str  = sprintf(image_list(il,8:end), '%s');
    if strcmpi(str(1),'.')
        row = 1;%novel image
    else
        row = 2;%manipulated or repeat image
    end
    %     else
    %         str  = sprintf(image_list(il,3:end), '%s');
    %         if strcmpi(str(1),'.')
    %             row = 1;%novel image
    %         else
    %             row = 2;%manipulated or repeat image
    %         end
    %end
    order(row,num) = il;
end

image_names = cell(2,36);

if strcmp(image_list(1,11),' ') %forgot to add 0 for image nums < 10
    for img = 1:72
        period = strfind(image_list(img,:),'.');
        imgnam =image_list(img,6:period-1);
        if ~isempty(strfind(imgnam,'p')) || ~isempty(strfind(imgnam,'r')) || ~isempty(strfind(imgnam,'m'))
            imgnum =str2double(imgnam(1:end-1));
            image_names{2,imgnum} = strtrim(image_list(img,:));
        else
            imgnum =str2double(imgnam);
            image_names{1,imgnum} = strtrim(image_list(img,:));
        end
    end
else
    for img = 1:36;
        %strtrim removes extra spaces if they exist in the string names
        image_names{1,img} = strtrim(image_list(order(1,img),:)); %1st presentation
        image_names{2,img} = strtrim(image_list(order(2,img),:));%2nd presentation
    end
end
userdata.image_names = image_names;%save to can acess across objects
userdata.last_image = 0; %no images have been loaded yet

set(hObject,'UserData',userdata);%store structure array to pushbutton object

%% load the next image
% --- Executes on button press in pushbutton1.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
loadeddata  = get(pub1_h,'UserData');
last_image = loadeddata.last_image;
if last_image == 36
    return
else
    last_image = last_image+1;
    loadeddata.last_image = last_image;
    set(pub1_h,'UserData',loadeddata);%store
end


pathname = [loadeddata.folder_name '\']; %where the images are  stored
image_names = loadeddata.image_names; %the image names

if   ~isempty(strfind(image_names{2,last_image}(6:8),'p'))
    type = 2; %familiar image
    slashes = strfind(pathname,'\');
    %     repeated_dir = [pathname(1:slashes(end-1)) 'Repeat\ROI\'];
elseif ~isempty(strfind(image_names{2,last_image}(6:8),'r'))
    type = 3; %replaced
elseif   ~isempty(strfind(image_names{2,last_image}(6:8),'m'))
    type = 4; %moved image
else
    error('Image type not found')
end

img1 = imread([pathname image_names{1,last_image}]);
if type == 2%familiar image has ROI preselected which we can easily find
    %img2 = imread([repeated_dir image_names{2,last_image}(1:end-4) 'r.bmp']);
    img2 = imread([pathname image_names{2,last_image}]);
else
    img2 = imread([pathname image_names{2,last_image}]);
end

img_pixel_difference = abs(rgb2gray(img1)-rgb2gray(img2));
img_pixel_difference(img_pixel_difference < 5) = 0;

if type == 3 %replaced
    [r, c] = find(img_pixel_difference ~=0); %r and c are the row and column indices for the matrix replacedchange where the value=0.
    ROI = [min(c) max(c) 600-min(r) 600-max(r)]; %draw ROI where pixels changed
    if ROI(4) < ROI(3) %swap due to flipping below central axis
        ROI = [ROI(1:2) ROI(4) ROI(3)];
    end
elseif type == 4 %moved
    [r, c] = find(img_pixel_difference ~=0); %r and c are the row and column indices for the matrix replacedchange where the value=0.
    IDX = kmeans([r c], 2);%there should be 2 clusters, the old location and the new location
    %arbiratrily asign ROI1 and ROI2 as the 2 different clusters. Can Swap later
    ROI{1} = [min(c(IDX == 1)) max(c(IDX == 1)) 600-min(r(IDX == 1)) 600-max(r(IDX == 1))]; %draw ROI where pixels changed
    ROI{2} = [min(c(IDX == 2)) max(c(IDX == 2)) 600-min(r(IDX == 2)) 600-max(r(IDX == 2))]; %draw ROI where pixels changed
    if ROI{1}(4) < ROI{1}(3) %swap due to flipping below central axis
        ROI{1} = [ROI{1}(1:2) ROI{1}(4) ROI{1}(3)];
    end
    if ROI{2}(4) < ROI{2}(3) %swap due to flipping below central axis
        ROI{2} = [ROI{2}(1:2) ROI{2}(4) ROI{2}(3)];
    end
else %for repeat set ROI to blank
    ROI = [];
end


axes(handles.axes1);
cla
hold on
handles.h1=image(flipud(img1));
if type == 4 %moved has 2 ROIs
    plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
        [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
elseif type == 3 %replaced
    plot([ROI(1) ROI(2) ROI(2) ROI(1) ROI(1)],...
        [ROI(3) ROI(3) ROI(4) ROI(4) ROI(3)],'r');
end
hold off
box off
axis off
axis equal
set(handles.text3,'String',['Novel: ' image_names{1,last_image}])

axes(handles.axes2);
cla
hold on
handles.h2=image(flipud(img2));
if type == 4 %moved has 2 ROIs
    plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
        [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
elseif type == 3 %replaced ROI
    plot([ROI(1) ROI(2) ROI(2) ROI(1) ROI(1)],...
        [ROI(3) ROI(3) ROI(4) ROI(4) ROI(3)],'r');
end
hold off
box off
axis off
axis equal

if type == 2
    set(handles.text4,'String',['Familiar: ' image_names{2,last_image}])
elseif type == 3
    set(handles.text4,'String',['Replaced: ' image_names{2,last_image}])
elseif type == 4
    set(handles.text4,'String',['Moved: ' image_names{2,last_image}])
end

axes(handles.axes3);
cla
hold on
handles.h3=image(flipud(img_pixel_difference));
if type == 4 %moved has 2 ROIs
    plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
        [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
elseif type == 3 %replaced ROI
    plot([ROI(1) ROI(2) ROI(2) ROI(1) ROI(1)],...
        [ROI(3) ROI(3) ROI(4) ROI(4) ROI(3)],'r');
end
hold off
axis equal
set(handles.text3,'String','Pixel Difference')
box off
axis off

if ~iscell(ROI) %for replaced images
    ROI = {ROI};
end
%store ROIs to user data for this pushbutton
userdata = get(hObject,'UserData');
if isempty(userdata)
    set(hObject,'UserData',{ROI});
else
    set(hObject,'UserData',[userdata {ROI}]);
end


%% Save ROI Data to File
% --- Executes on button press in pushbutton3.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('ROIs.mat','Save file name');
pub3_h = findobj('Tag','pushbutton3');%get the hanlde for pushbutton3 that has ROI data
ROIs  = get(pub3_h,'UserData');
if length(ROIs) < 36
    disp('Warning Not all Images Have ROIs selected')
end
save([path file],'ROIs')

%% IF moved ROIs are labeled wrong switch so they are labeled correctly
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pub3_h = findobj('Tag','pushbutton3');%get the hanlde for pushbutton3 that has ROI data
ROIs  = get(pub3_h,'UserData');

if length(ROIs{end}) == 1%then it's not  a moved image
    return%we don't want to do anything
end

%if it is a moved image we want to swap the ROIs so the 1st ROI is the
%origional position in red and the 2nd ROI is the new location in green

ROI{1} = ROIs{end}{2};
ROI{2} = ROIs{end}{1};

%once swapped save to object
ROIs{end} = ROI;
set(pub3_h,'UserData',ROIs);

%then replot to visual confirm the change
pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
loadeddata  = get(pub1_h,'UserData');
last_image = loadeddata.last_image;

pathname = [loadeddata.folder_name '\']; %where the images are  stored
image_names = loadeddata.image_names; %the image names

type = 4; %moved image

img1 = imread([pathname image_names{1,last_image}]);
img2 = imread([pathname image_names{2,last_image}]);

img_pixel_difference = abs(rgb2gray(img1)-rgb2gray(img2));
img_pixel_difference(img_pixel_difference < 5) = 0;

axes(handles.axes1);
cla
hold on
handles.h1=image(flipud(img1));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
box off
axis off
axis equal

axes(handles.axes2);
cla
hold on
handles.h2=image(flipud(img2));

plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
box off
axis off
axis equal

axes(handles.axes3);
cla
hold on
handles.h3=image(flipud(img_pixel_difference));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
axis equal
box off
axis off

%% If Novel ROI is drawn wrong for moved or replaced images, or if need to draw repeat ROI
% this button allows you to draw new rectangle for the ROI
% --- Executes on button press in pushbutton4.
function pushbutton5_Callback(hObject, eventdata, handles)
axes(handles.axes1); %I don't think this actually limits you to the first axis
rect = getrect;

pub3_h = findobj('Tag','pushbutton3');%get the hanlde for pushbutton3 that has ROI data
ROIs  = get(pub3_h,'UserData');

if length(ROIs{end}) == 1 %replaced image
    ROI{1} = ROIs{end}{1};
elseif  length(ROIs{end}) == 2 %moved image
    ROI{1} = ROIs{end}{1};
    ROI{2} = ROIs{end}{2};
else %repeat image with no ROI
    ROI{1} = [];
end

xmin = rect(1);
xmax = rect(1)+rect(3);
ymin = rect(2);
ymax = rect(2)+rect(4);
ROI{1} = [xmin xmax ymin ymax]; %the new ROI has been chaned, regardless of type
%once swapped save to object
if iscell(ROI)
    ROIs{end} = ROI;
else
    ROIs{end} = {ROI};
end
set(pub3_h,'UserData',ROIs);

%then replot to visual confirm the change
pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
loadeddata  = get(pub1_h,'UserData');
last_image = loadeddata.last_image;

pathname = [loadeddata.folder_name '\']; %where the images are  stored
image_names = loadeddata.image_names; %the image names

img1 = imread([pathname image_names{1,last_image}]);
img2 = imread([pathname image_names{2,last_image}]);

img_pixel_difference = abs(rgb2gray(img1)-rgb2gray(img2));
img_pixel_difference(img_pixel_difference < 5) = 0;

axes(handles.axes1);
cla
hold on
handles.h1=image(flipud(img1));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
if  length(ROIs{end}) == 2 %moved image
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
end
hold off
box off
axis off
axis equal

axes(handles.axes2);
cla
hold on
handles.h2=image(flipud(img2));

plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
if  length(ROIs{end}) == 2 %moved image
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
end
hold off
box off
axis off
axis equal

axes(handles.axes3);
cla
hold on
handles.h3=image(flipud(img_pixel_difference));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
if  length(ROIs{end}) == 2 %moved image
    plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
        [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
end
hold off
axis equal
box off
axis off



%% If Moved ROI 2 is drawn wrong
% this button allows you to draw new rectangle for the ROI
% --- Executes on button press in pushbutton4.
function pushbutton6_Callback(hObject, eventdata, handles)

pub3_h = findobj('Tag','pushbutton3');%get the hanlde for pushbutton3 that has ROI data
ROIs  = get(pub3_h,'UserData');

if length(ROIs{end}) == 1%then it's not  a moved image, and ROI 2 does not apply
    return%we don't want to do anything
end

axes(handles.axes2)
rect = getrect; %I don't think this actually limits you to the second axis

ROI{1} = ROIs{end}{1};
ROI{2} = ROIs{end}{2};

xmin = rect(1);
xmax = rect(1)+rect(3);
ymin = rect(2);
ymax = rect(2)+rect(4);
ROI{2} = [xmin xmax ymin ymax]; %the new ROI has been chaned, regardless of type

%once swapped save to object
if iscell(ROI)
    ROIs{end} = ROI;
else
    ROIs{end} = {ROI};
end
set(pub3_h,'UserData',ROIs);

%then replot to visual confirm the change
pub1_h = findobj('Tag','pushbutton1');%get the hanlde for pushbutton1
loadeddata  = get(pub1_h,'UserData');
last_image = loadeddata.last_image;

pathname = [loadeddata.folder_name '\']; %where the images are  stored
image_names = loadeddata.image_names; %the image names

img1 = imread([pathname image_names{1,last_image}]);
img2 = imread([pathname image_names{2,last_image}]);

img_pixel_difference = abs(rgb2gray(img1)-rgb2gray(img2));
img_pixel_difference(img_pixel_difference < 5) = 0;

axes(handles.axes1);
cla
hold on
handles.h1=image(flipud(img1));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
box off
axis off
axis equal

axes(handles.axes2);
cla
hold on
handles.h2=image(flipud(img2));

plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
box off
axis off
axis equal

axes(handles.axes3);
cla
hold on
handles.h3=image(flipud(img_pixel_difference));
plot([ROI{1}(1) ROI{1}(2) ROI{1}(2) ROI{1}(1) ROI{1}(1)],...
    [ROI{1}(3) ROI{1}(3) ROI{1}(4) ROI{1}(4) ROI{1}(3)],'r');
plot([ROI{2}(1) ROI{2}(2) ROI{2}(2) ROI{2}(1) ROI{2}(1)],...
    [ROI{2}(3) ROI{2}(3) ROI{2}(4) ROI{2}(4) ROI{2}(3)],'g');
hold off
axis equal
box off
axis off
