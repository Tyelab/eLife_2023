function varargout = calibrationgui(varargin)
% CALIBRATIONGUI MATLAB code for calibrationgui.fig
%      CALIBRATIONGUI, by itself, creates a new CALIBRATIONGUI or raises the existing
%      singleton*.
%
%      H = CALIBRATIONGUI returns the handle to a new CALIBRATIONGUI or the handle to
%      the existing singleton*.
%
%      CALIBRATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATIONGUI.M with the given input arguments.
%
%      CALIBRATIONGUI('Property','Value',...) creates a new CALIBRATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibrationgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibrationgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibrationgui

% Last Modified by GUIDE v2.5 09-Nov-2016 17:19:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibrationgui_OpeningFcn, ...
                   'gui_OutputFcn',  @calibrationgui_OutputFcn, ...
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


% --- Executes just before calibrationgui is made visible.
function calibrationgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibrationgui (see VARARGIN)

% Choose default command line output for calibrationgui
handles.output.hObject = hObject;

% Parameters
radius_range = [111 200]; % If using imcircles, then rmax should be < 3*r_min and rmax - rmin < 100
handles.defaultRadius = 270;
handles.maxFibers = 8; % Original code = 7
num_bit = 16;

% Initialization
handles.ellipses = {};
% handles.cmap = hsv(handles.maxFibers);
% handles.cmap = handles.cmap(randperm(handles.maxFibers),:); % Don't randomize, EK 2019/07/17
handles.cmap = [0.9 0 0; 0.3 0.3 1; 0.8 0.4 0; 0 0.7 0.2; 1 0 1; 0.4 0.4 0.4; 0.2 0.6 1; 0.5 0 0.5]; % Predefine colors, EK 2019/07/17
% image(reshape(handles.cmap, [size(handles.cmap, 1), 1, 3])); % Plot colors

% Process the input image
handles.image = varargin{1};
max_value = max(handles.image(:));
if max_value >= 2^num_bit
    warning(sprintf('WARNING: Calibration image has saturated pixels for %d bits', num_bit));
end
% disp(sprintf('Max of calibration image uses ' num2str(max_value*100/(2^num_bit)))) '% of camera dynamic range (assuming %d bits).', num_bit));
fprintf('Max of calibration image uses %.1f%% of camera dynamic range (assuming %d bits and 1x1 binning).\n', max_value*100/(2^num_bit), num_bit);
if max_value*100/2^num_bit < 10
    warning('Consider increasing the light power or reducing the acquisition rate for better signal to noise.');
end
imagesc(handles.image, 'Parent', handles.img_ax);
handles.frameSize = size(handles.image);

% % Find initial circles: Defer
% [centers, radii, metric] = imfindcircles(handles.image, radius_range);
% radii = radii * 0.9; % Shrink slightly for automatic circles

% Added to use prior coordinates to simplify placements: EK 2019/07/17
old_topleft = [  ...
812.00	652.00	270.00	270.00	
1160.00	684.00	270.00	270.00	
1328.00	996.00	270.00	270.00	
908.00	992.00	270.00	270.00	
% 410.00	322.00	135.00	135.00	
% 586.00	340.00	135.00	135.00	
% 666.00	498.00	135.00	135.00	
% 454.00	494.00	135.00	135.00	
% 578.50	652.50	135.00	135.00	
% 398.50	660.50	135.00	135.00	
  ];
% Prior 2x2 binning = 1024 x 1024, keep as is; For 4x4 halve and for 1x1 double
radii = old_topleft(:, 3)/2; % Prior output is as diameters
centers = old_topleft(:, 1:2) + repmat(radii, 1, 2); % Prior is as top left corners

if centers
    % Initialize ellipses
    for i = 1:size(old_topleft, 1)
%         r = radii(i) * 0.9; % moved above
%         handles = placeEllipse([centers(i,1:2) - radii(i), 2*radii(i), 2*radii(i)], i, handles);
        handles = placeEllipse(old_topleft(i, :), i, handles); % converts from topleft to center to topleft...
    end
else
    % No ellipses found. Place one in the middle.
    r = handles.defaultRadius;
    handles = placeEllipse([handles.frameSize(1)/2 - r/2, handles.frameSize(2)/2 - r/2, r, r], 1, handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibrationgui wait for user response (see UIRESUME)
uiwait(handles.calibrationgui);

function handles = placeEllipse(pos, num, handles)
center = [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2];
axes(handles.img_ax);
% numh = text(center(1) - 12, center(2) - 12, num2str(num), 'FontSize', 24); % prior code with nudging
numh = text(center(1), center(2), num2str(num), 'FontSize', 24); % just go with middle/center
set(numh, 'HorizontalAlignment', 'center');
h = imellipse(handles.img_ax, pos);
handles.ellipses{num} = h;

color = handles.cmap(num,:);
h.setColor(color);
numh.Color = color;
h.addNewPositionCallback(@(pos) updateNumber(pos, numh));

function updateNumber(pos, h)
center = [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2];
h.Position = [center(1) - 12, center(2) - 12];

% --- Outputs from this function are returned to the command line.
function varargout = calibrationgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.calibrationgui);

% --- Executes on button press in add_btn.
function add_btn_Callback(hObject, eventdata, handles)
% hObject    handle to add_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nfibers = length(handles.ellipses);

r = handles.defaultRadius;
handles = placeEllipse([handles.frameSize(1)/2 - r/2, handles.frameSize(2)/2 - r/2, r, r], nfibers + 1, handles);

% Update handles structure
guidata(hObject, handles);

function labels_txt_Callback(hObject, eventdata, handles)
% hObject    handle to labels_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labels_txt as text
%        str2double(get(hObject,'String')) returns contents of labels_txt as a double


% --- Executes during object creation, after setting all properties.
function labels_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labels_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in done_btn.
function done_btn_Callback(hObject, eventdata, handles)
% hObject    handle to done_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pre-check for number of valid ellipses
nfibers = 0;
for i = 1:length(handles.ellipses)
    if isvalid(handles.ellipses{i})
        nfibers = nfibers + 1;
    end
end

% % Get the number of labels and warn if there's a different number of
% % labels and fibers
% raw_labels = get(handles.labels_txt, 'String');
% [rows, cols] = size(raw_labels);
% labels = cell(rows, 1);
% for row = 1:rows
%     str = deblank(raw_labels(row,:));
%     if length(str) == 0
%         errordlg('List of labels cannot include blank lines');
%         return;
%     end
%     labels{row} = str;
% end
% 
% if length(labels) ~= nfibers
%     errordlg('Number of labels must be equal to the number of fibers');
%     return;
% end

% Just give number names: EK 2019/07/17
labels = cell(nfibers, 1);
for i = 1:nfibers
    labels{i} = sprintf('%d', i);
end

frsz = handles.frameSize;
colors = zeros(nfibers, 3);
masks = zeros([frsz nfibers]);

% Store this labeled image for saving later
figImg = getframe(gcf);

i = 0;
fprintf('Circle/Ellipse data: Pos/Left Pos/Top X-Diam Y-Diam\n')
for j = 1:nfibers
    e = handles.ellipses{j};
    if isvalid(e)
        i = i +1;
        pos = e.getVertices();
        x = pos(:,1); y = pos(:,2);
        mask = poly2mask(x, y, frsz(1), frsz(2));
        masks(:,:,i) = mask;
        colors(i,:) = e.getColor();
        fprintf('%.2f\t', e.getPosition);
    end
    fprintf('\n');
end
fprintf('\n');
    
handles.output.figImg = figImg;
handles.output.colors = colors;
handles.output.masks = masks;
handles.output.labels = labels;

% Update handles structure
guidata(hObject, handles);

close(handles.calibrationgui);

% --- Executes when user attempts to close calibrationgui.
function calibrationgui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to calibrationgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
