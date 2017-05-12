function varargout = Kennlinienverfahren(varargin)
% KENNLINIENVERFAHREN MATLAB code for Kennlinienverfahren.fig
% KENNLINIENVERFAHREN implements the Kennlinienverfahren () for the
% maximum-power-point (MPP) search. It provides a GUI where input parameters
% can be edited and the MPP search is visualized.
%
%      KENNLINIENVERFAHREN(file) creates a new KENNLINIENVERFAHREN or raises the
%      existing singleton. The value of file is used as input file for
%      the characteristic curves at different intensities. The file needs
%      to be a mat file containing the fields U, I and intensity. U and
%      intensity are 1xn and 1xm vectors respectively, while I is a nxm
%      matrix. I(i,j) must contain the current for voltage U(i) and intensity
%      intensity(j).
%
%      KENNLINIENVERFAHREN, by itself, creates a new KENNLINIENVERFAHREN GUI or raises the existing
%      singleton. It uses the file 'kennlinie.mat' as input.
%
%      H = KENNLINIENVERFAHREN returns the handle to a new KENNLINIENVERFAHREN or the handle to
%      the existing singleton. It uses the file 'kennlinie.mat' as input.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Kennlinienverfahren_OpeningFcn, ...
    'gui_OutputFcn',  @Kennlinienverfahren_OutputFcn, ...
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
end


% --- Executes just before Kennlinienverfahren is made visible.
function Kennlinienverfahren_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Kennlinienverfahren (see VARARGIN)

% Choose default command line output for Kennlinienverfahren
handles.output = hObject;

% determine input source from input parameter
datafile = 'kennlinie.mat'; % default
if nargin > 4
    uiwait(warndlg('Zu viele Input-Parameter. Nur der erste Parameter wird verwendet.', 'Warnung', 'modal'))
    datafile = varargin{1};
elseif nargin == 4
    datafile = varargin{1};
end


%% load I, U, intensity from file
try
    handles.K = load(datafile);
catch
    uiwait(errordlg(['Datei ' datafile ' konnte nicht geöffnet werden.'], 'Fehler', 'modal'));
    close(hObject);
    return;
end

% check if necessary fields exist
if ~isfield(handles.K,'U') ...
        || ~isfield(handles.K,'I') ...
        || ~isfield(handles.K,'intensity')
    uiwait( ...
        errordlg( ...
        ['Die Datei ' datafile ...
        ' enthält nicht die benötigten Variablen intensity, U und I.'], ...
        'Fehler', ...
        'modal' ...
        ) ...
        );
    close(hObject);
    return;
end

% check if fields have consistent sizes
if size(handles.K.U,1) ~= 1 ...
        || size(handles.K.intensity, 1) ~= 1 ...
        || size(handles.K.U,2) ~= size(handles.K.I,1) ...
        || size(handles.K.intensity,2) ~= size(handles.K.I,2)
    uiwait( ...
        errordlg( ...
        ['Die Datei ' datafile ...
        ' enthält Daten, die nicht dem vorgegebenen Format ' ...
        'entsprechen. Benötigte Matrixgrößen: size(U) = [1 n], ' ...
        'size(I) = [n m], size(intensity) = [1 m].'], ...
        'Fehler', ...
        'modal'));
    close(hObject);
    return;
end

%% set default parameters
handles.intensity = 50; % intensity of incident radiation
handles.k = 0.8;        % constant in range (0, 1) that determines distance between reversal points of maximum-power-point search
handles.R = 280;        % resistance
handles.dR = 7;         % resistance change per time step
handles.s = -1;         % direction of resistance change (- -> U grows smaller)
handles.T = 1;          % duration of 1 time step / timer period
handles.dRmax = 100;    % maximum allowed value of dR
handles.Rmin = 0;       % minimum allowed value of R
handles.Rmax = 1000;    % maximum allowed value of R
handles.checkvalue = 0;
handles.minStepsBetweenFails = 5; % minimum number of time steps between non-positive resistance values (and corresponding error handling) before warning occurs

handles.maxU = max(handles.K.U);

%% calculate initial characteristic curve + current and voltage values
% characteristic curve for given intensity (handles.intensity)
I = interp2(handles.K.intensity, handles.K.U, handles.K.I, handles.intensity, handles.K.U);
% calculate current position on characteristic curve for given resistance (handles.R)
intersect = @(u) interp2(handles.K.intensity, handles.K.U, handles.K.I, handles.intensity, u) - u/handles.R;
handles.U_current = fzero(intersect, [0 handles.maxU]);
handles.I_current = handles.U_current/handles.R;
handles.currentUmax = handles.U_current;
handles.currentImax = handles.I_current;
% calculate ranges for static axes limits from maximum current and power values
handles.Irange = [0 ceil(max(max(handles.K.I)))];
handles.Prange = [0 ceil(0.1*max(max(handles.K.I' .* repmat(handles.K.U, size(handles.K.I, 2), 1))))*10];

%% set up plots
% plots on left y axis (current)
handles.ax1 = gca;
handles.ax1.YLim = handles.Irange;
handles.ax1.YColor = 'b';
xlabel('U (V)');
ylabel('I (A)');
hold on
% plot characteristic curve (current and voltage)
handles.plots.kennlinie = plot(handles.K.U, I);
% plot current point
handles.plots.point = plot(handles.U_current, handles.I_current, 'kx');
% plot I = U/R as a line, but only for values smaller than current voltage
cond = (handles.K.U <= handles.U_current);
handles.plots.rline = plot(handles.K.U(cond), handles.K.U(cond)/handles.R, 'k:');
hold off

% plots on right y axis (power)
% set up second y axis
handles.ax2 = axes('Position',get(handles.ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','r', ...
    'YLim', handles.Prange, ...
    'XTickLabel', []);
% link axes so that all plots share the same x axis
linkaxes([handles.ax1 handles.ax2],'x');
ylabel('P (W)');
hold on % we are now in a different plot so we need this again
% plot characteristic curve (power and voltage)
handles.plots.power = plot(handles.K.U, handles.K.U .* I', ...
    'Parent', handles.ax2, 'LineStyle', '-', 'Color', [1 0 0]);
% plot current point
handles.plots.powerpoint = plot(handles.U_current, handles.U_current .* handles.I_current, ...
    'Parent', handles.ax2, 'Marker', 'x', 'Color', 'k');
hold off


%% create timer
handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly.
    'Period', 1, ...                        % Initial period is 1 sec.
    'TimerFcn', {@next_step,hObject}); % Specify callback function.

%% display initial parameter values
set(handles.dRedit, 'String', num2str(handles.dR));
set(handles.resistanceedit, 'String', num2str(handles.R));
set(handles.intensityedit, 'String', num2str(handles.intensity));
set(handles.kedit, 'String', num2str(handles.k));
set(handles.speedslider, 'Value', handles.T);
set(handles.speededit, 'String', num2str(handles.T));
set(handles.intensityslider, 'Min', min(handles.K.intensity));
set(handles.intensityslider, 'Max', max(handles.K.intensity));
set(handles.intensityslider, 'Value', handles.intensity);
set(handles.scalemenu, 'Value', 1);
if handles.s == 1
    set(handles.sposbutton, 'Value', handles.sposbutton.Max);
else % handles.s == -1
    set(handles.snegbutton, 'Value', handles.snegbutton.Max);
end
set(handles.udisp, 'String', num2str(handles.U_current));
set(handles.idisp, 'String', num2str(handles.I_current));
set(handles.pdisp, 'String', num2str(handles.U_current*handles.I_current));

%% define error/info messages
handles.msg.nan = 'Bitte nur Zahlen eingeben.';
handles.msg.dRexceedsMax = ['dR darf maximal ' num2str(handles.dRmax) ' Ohm betragen.'];
handles.msg.intensityBelowMin = ['Intensität muss größer als ' num2str(min(handles.K.intensity)) ' W/m^2 sein.'];
handles.msg.intensityExceedsMax = ['Intensität muss kleiner als ' num2str(max(handles.K.intensity)) ' W/m^2 sein.'];
handles.msg.kOutsideRange = 'k muss zwischen 0 und 1 liegen.';
handles.msg.RBelowMin = ['Widerstand muss größer ' num2str(handles.Rmin) ' Ohm sein.'];
handles.msg.RExceedsMax = ['Widerstand muss kleiner ' num2str(handles.Rmax) ' Ohm sein.'];
handles.msg.speedBelowMin = ['Wert muss größer als ' num2str(handles.speedslider.Min) ' s sein.'];
handles.msg.speedExceedsMax = ['Wert muss kleiner als ' num2str(handles.speedslider.Max) ' s sein.'];
handles.msg.dRTooSmall = 'Der Wert von dR ist für die gewählte Intensität zu niedrig, um vernünftige Ergebnisse zu erzielen.';
handles.msg.help = {'Diese Anwendung implementiert das Kennlinienverfahren (Boehringer 1969) für die Maximum-Power-Point-(MPP-)Suche.', ...
    '', ...
    'Parameter:', ...
    'Intensität: Intensität der einfallenden Strahlung. Der Verlauf der Kennlinie hängt hiervon ab.', ...
    'R: Lastwiderstand.', ...
    'dR: Widerstandsänderung pro Zeitschritt.', ...
    'Richtung: Suchrichtung. Links = Widerstand wird erniedrigt. Rechts = Widerstand wird erhöht.', ...
    'k: Vergleichsfaktor mit 0 < k < 1, der den Abstand der Umkehrpunkte bei der MPP-Suche bestimmt.', ...
    '', ...
    'Anzeigeoptionen:', ...
    'Autoscaling: An = automatische Skalierung der Plots. Aus = fester y-Bereich zum besseren Vergleich bei verschiedenen Intensitäten.', ...
    'Dauer eines Zeitschritts: Angabe in Sekunden.'};


%% Update handles structure
guidata(hObject, handles);
end

function next_step(~,~,hObject)
handles = guidata(hObject);

% update current point
handles = update_point(handles);

% update maximum values for current cycle
if handles.U_current > handles.currentUmax
    handles.currentUmax = handles.U_current;
end
if handles.I_current > handles.currentImax
    handles.currentImax = handles.I_current;
end

% check for direction change conditions
if handles.s == -1
    if handles.U_current <= handles.k*handles.currentUmax
        % change direction
        handles.s = 1;
        set(handles.sposbutton, 'Value', handles.sposbutton.Max);
        % starting new cycle: restart maximum tracking
        handles.currentImax = handles.I_current;
        handles.currentUmax = handles.U_current;
    end
else % handles.s == 1
    if handles.I_current <= handles.k*handles.currentImax
        % change direction
        handles.s = -1;
        set(handles.snegbutton, 'Value', handles.snegbutton.Max);
        % starting new cycle: restart maximum tracking
        handles.currentImax = handles.I_current;
        handles.currentUmax = handles.U_current;
    end
end

% update resistance
handles.R = handles.R + handles.s*handles.dR;
if handles.R <= 0
    % can't deal with non-positive resistance
    % undo last step
    handles.R = handles.R - 2*handles.s*handles.dR;
    % change direction
    handles.s = handles.s*(-1);
    if handles.checkvalue > 0
        % this error has now happened once too often in a number of time steps
        % warn user
        warn(handles, handles.msg.dRTooSmall);
        % assume that user changes something and act as if fail never
        % happened
        handles.checkvalue = 0;
    else
        handles.checkvalue = handles.minStepsBetweenFails;
    end
else
    handles.checkvalue = handles.checkvalue - 1;
    if handles.checkvalue < 0
        handles.checkvalue = 0;
    end
end
set(handles.resistanceedit, 'String', num2str(handles.R));

% update handles
guidata(hObject,handles);
end



% --- Outputs from this function are returned to the command line.
function varargout = Kennlinienverfahren_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
catch exception
    %ignore
end
end

% --- Executes on button press in startbutton.
function startbutton_Callback(hObject, ~, handles)
% hObject    handle to startbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% toggle timer status and button text
if strcmp(get(handles.timer, 'Running'), 'off')
    start(handles.timer);
    set(hObject,'String', 'Pause');
else
    stop(handles.timer);
    set(hObject,'String', 'Start');
end
end



% --- Executes on slider movement.
function speedslider_Callback(hObject, ~, handles)
% hObject    handle to speedslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.T = round(get(hObject, 'Value')*10)/10;
set(handles.speededit, 'String', num2str(handles.T));
if strcmp(get(handles.timer, 'Running'), 'on')
    % timer needs to be stopped to change period
    stop(handles.timer);
    set(handles.timer, 'Period', handles.T);
    start(handles.timer);
else
    set(handles.timer, 'Period', handles.T);
end

guidata(hObject, handles);
end


function intensityedit_Callback(hObject, ~, handles)
% hObject    handle to intensityedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2double(get(hObject, 'String'));
% check for invalid input
if isnan(val)
    set(hObject, 'String', num2str(handles.intensity));
    warn(handles, handles.msg.nan)
elseif val < min(handles.K.intensity)
    set(hObject, 'String', num2str(handles.intensity));
    warn(handles, handles.msg.intensityBelowMin)
elseif val > max(handles.K.intensity)
    set(hObject, 'String', num2str(handles.intensity));
    warn(handles, handles.msg.intensityExceedsMax)
else
    % input is valid
    handles.intensity = val;
end
set(handles.intensityslider, 'Value', handles.intensity);

% calculate new characteristic curve and update plots
I = interp2(handles.K.intensity, handles.K.U, handles.K.I, handles.intensity, handles.K.U);
set(handles.plots.kennlinie, 'YData', I);
set(handles.plots.power, 'YData', handles.K.U .* I');

handles = update_point(handles);

guidata(hObject, handles);
end


function resistanceedit_Callback(hObject, ~, handles)
% hObject    handle to resistanceedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = str2double(get(hObject, 'String'));
% check for invalid input
if isnan(val)
    set(hObject, 'String', num2str(handles.R));
    warn(handles, handles.msg.nan)
elseif val <= handles.Rmin
    set(hObject, 'String', num2str(handles.R));
    warn(handles, handles.msg.RBelowMin)
elseif val >= handles.Rmax
    set(hObject, 'String', num2str(handles.R));
    warn(handles, handles.msg.RExceedsMax)
else
    % input is valid
    handles.R = val;
    handles = update_point(handles);
    % reset maxima of observed current and voltage
    handles.currentUmax = handles.U_current;
    handles.currentImax = handles.I_current;
end

guidata(hObject, handles);
end


function dRedit_Callback(hObject, ~, handles)
% hObject    handle to dRedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dR = str2double(get(hObject, 'String'));
% check if input is valid
if isnan(dR)
    set(hObject, 'String', num2str(handles.dR));
    warn(handles, handles.msg.nan)
elseif abs(dR) > handles.dRmax
    set(hObject, 'String', num2str(handles.dR));
    warn(handles, handles.msg.dRexceedsMax)
elseif dR < 0
    % input is valid but change of direction required
    handles.s = -handles.s;
    handles.dR = -dR;
    set(hObject, 'String', num2str(handles.dR));
else
    % input is valid
    handles.dR = dR;
end

guidata(hObject, handles);
end


function kedit_Callback(hObject, ~, handles)
% hObject    handle to kedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = str2double(get(hObject, 'String'));
% check if input is invalid
if isnan(k)
    set(hObject, 'String', num2str(handles.k));
    warn(handles, handles.msg.nan)
else
    % input is valid
    if k <= 0 || k >= 1
        set(hObject, 'String', num2str(handles.k));
        warn(handles, handles.msg.kOutsideRange);
    else
        handles.k = k;
    end
end

guidata(hObject, handles);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if strcmp(get(handles.timer, 'Running'), 'on')
        stop(handles.timer);
    end
    % Destroy timer
    delete(handles.timer)
catch exception
    % ignore, just means that no timer exists -> close anyway
end
% Hint: delete(hObject) closes the figure
delete(hObject);
end



function speededit_Callback(hObject, ~, handles)
% hObject    handle to speededit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = round(str2double(get(hObject, 'String'))*10)/10;
% check for invalid input
if isnan(val)
    set(hObject, 'String', num2str(handles.T));
    warn(handles, handles.msg.nan)
elseif val < handles.speedslider.Min
    set(hObject, 'String', num2str(handles.T));
    warn(handles, handles.msg.speedBelowMin)
elseif val > handles.speedslider.Max
    set(hObject, 'String', num2str(handles.T));
    warn(handles, handles.msg.speedExceedsMax)
else
    % input is valid
    handles.T = val;
    set(handles.speedslider, 'Value', handles.T);
    if strcmp(get(handles.timer, 'Running'), 'on')
        % timer needs to be stopped to change period
        stop(handles.timer);
        set(handles.timer, 'Period', handles.T);
        start(handles.timer);
    else
        set(handles.timer, 'Period', handles.T);
    end
end

guidata(hObject, handles);
end


% --- Executes on selection change in scalemenu.
function scalemenu_Callback(hObject, ~, handles)
% hObject    handle to scalemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
% apply scaling mode to the plots
switch val
    case 2 % an
        set(handles.ax1, 'YLimMode', 'auto');
        set(handles.ax2, 'YLimMode', 'auto');
    case 1 % aus
        set(handles.ax1, 'YLim', handles.Irange);
        set(handles.ax2, 'YLim', handles.Prange);
end
% Hints: contents = cellstr(get(hObject,'String')) returns scalemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scalemenu
guidata(hObject, handles);
end


% creates a warning dialog that pauses the main window
function warn(handles, message)
% handles   structure with handles and user data
% message   message to be displayed
stop(handles.timer);
set(handles.startbutton, 'String', 'Start');
uiwait(warndlg(message, 'Fehlerhafte Eingabe', 'modal'));
end


% --- Executes on slider movement.
function intensityslider_Callback(hObject, ~, handles)
% hObject    handle to intensityslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.intensity = round(get(hObject, 'Value'));
set(handles.intensityedit, 'String', num2str(handles.intensity));

I = interp2(handles.K.intensity, handles.K.U, handles.K.I, handles.intensity, handles.K.U);
set(handles.plots.kennlinie, 'YData', I);
set(handles.plots.power, 'YData', handles.K.U .* I');

handles = update_point(handles);

guidata(hObject, handles);
end


% --- Executes on button press in snegbutton.
function snegbutton_Callback(hObject, ~, handles)
% hObject    handle to snegbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.s = -1;
% reset maxima of observed current and voltage
handles.currentUmax = handles.U_current;
handles.currentImax = handles.I_current;

guidata(hObject, handles);
end

% --- Executes on button press in sposbutton.
function sposbutton_Callback(hObject, ~, handles)
% hObject    handle to sposbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.s = 1;
% reset maxima of observed current and voltage
handles.currentUmax = handles.U_current;
handles.currentImax = handles.I_current;

guidata(hObject, handles);
end


% --- Executes on button press in helpbutton.
function helpbutton_Callback(~, ~, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
title = 'Info';
if strcmp(get(handles.timer, 'Running'), 'off')
    uiwait(msgbox(handles.msg.help, title, 'modal'));
else
    stop(handles.timer);
    uiwait(msgbox(handles.msg.help, title, 'modal'));
    start(handles.timer);
end
end

% --- Executes during object creation, after setting all properties.
function default_CreateFcn(hObject, ~, ~)
% hObject    handle to udisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function handles = update_point(handles)
% find U and I for current resistance value
% at intersection of characteristic curve (interpolated from available
% data)and I(U) = U/R
intersect = @(u) interp2(handles.K.intensity, handles.K.U, handles.K.I, handles.intensity, u) - u/handles.R;
handles.U_current = fzero(intersect, [0 handles.maxU]);
handles.I_current = handles.U_current/handles.R;
% current power, just needed for plotting
P_current = handles.U_current*handles.I_current;

% update display
set(handles.udisp, 'String', num2str(handles.U_current));
set(handles.idisp, 'String', num2str(handles.I_current));
set(handles.pdisp, 'String', num2str(P_current));

% update plots of current point
set(handles.plots.point, 'XData', handles.U_current);
set(handles.plots.point, 'YData', handles.I_current);
set(handles.plots.powerpoint, 'XData', handles.U_current);
set(handles.plots.powerpoint, 'YData', handles.U_current*handles.I_current);

% update plot of I = U/R
U = 0:0.1:handles.U_current;
I = U/handles.R;
set(handles.plots.rline, 'XData', U);
set(handles.plots.rline, 'YData', I);
end
