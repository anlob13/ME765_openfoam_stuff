%% ==============================================================
% PERIODIC HILL VELOCITY PROFILES
% ==============================================================

clear;
clc;
close all;

%% USER SETTINGS

caseDir = '.';
timeName = '0.08';

% OpenFOAM coordinates are already normalized by h
Ub = 1.0;       % bulk velocity

% Requested x/h stations
xStations = [ ...
    0 ...
    0.5 ...
    1 ...
    1.5 ...
    2 ...
    2.5 ...
    3 ...
    4 ...
    5 ...
    6 ...
    7 ...
    8 ];

% Width of x/h sampling window
xWindow = 0.03;

%% ==============================================================
% READ CELL CENTRES
% ==============================================================

fprintf('Reading C...\n');

Cfile = fullfile(caseDir,timeName,'C');

C = readOpenFOAMVectorField(Cfile);

fprintf('C: %d vectors\n',size(C,1));

%% ==============================================================
% READ VELOCITY
% ==============================================================

fprintf('Reading U...\n');

Ufile = fullfile(caseDir,timeName,'U');

U = readOpenFOAMVectorField(Ufile);

fprintf('U: %d vectors\n',size(U,1));

%% ==============================================================
% CHECK
% ==============================================================

if size(C,1) ~= size(U,1)

    error( ...
        'C and U contain different numbers of cells: %d vs %d', ...
        size(C,1),size(U,1));

end

%% ==============================================================
% EXTRACT COORDINATES AND VELOCITY
% ==============================================================

x  = C(:,1);
y  = C(:,2);
z  = C(:,3);

Ux = U(:,1);
Uy = U(:,2);
Uz = U(:,3);

%% ==============================================================
% VELOCITY PROFILES OVER HILL
% ==============================================================

Nx = 120;
Ny = 40;
Nz = 60;

% Check mesh size
if size(C,1) ~= Nx*Ny*Nz
    error('Unexpected mesh size.');
end

% ---------------------------------------------------------------
% Reshape fields
%
% Based on your C file, x varies fastest, then y, then z.
% ---------------------------------------------------------------

X  = reshape(x,  [Nx Ny Nz]);
Y  = reshape(y,  [Nx Ny Nz]);

UX = reshape(Ux, [Nx Ny Nz]);

%% ==============================================================
% PLOT
% ==============================================================

figure('Color','w');
hold on;

% Hill
xhHill = linspace(0,9,1000);
ywHill = arrayfun(@hillHeight,xhHill);

plot( ...
    xhHill, ...
    ywHill, ...
    'k-', ...
    'LineWidth',2);

% Top wall
plot( ...
    [0 9], ...
    [3.035 3.035], ...
    'k--', ...
    'LineWidth',1);

%% ==============================================================
% PROFILE SCALE
% ==============================================================

profileScale = 0.5;

%% ==============================================================
% EXTRACT EACH PROFILE
% ==============================================================

for k = 1:length(xStations)

    xh = xStations(k);

    %% ----------------------------------------------------------
    % Hill height
    % ----------------------------------------------------------

    yw = hillHeight(xh);

    %% ----------------------------------------------------------
    % Find nearest x-index for EVERY y layer
    % ----------------------------------------------------------

    nLevels = Ny;

    yProfile = zeros(nLevels,1);
    uProfile = zeros(nLevels,1);

    for j = 1:Ny

        % x coordinates of this wall-normal layer
        xLayer = squeeze(X(:,j,:));

        % Average x over z
        xLayerMean = mean(xLayer,2);

        % Find closest x location
        [~,i] = min(abs(xLayerMean - xh));

        % Average velocity over all z
        uLayer = squeeze(UX(i,j,:));

        uProfile(j) = mean(uLayer);

        % Average y over z
        yLayer = squeeze(Y(i,j,:));

        yProfile(j) = mean(yLayer);

    end

    %% ----------------------------------------------------------
    % Sort profile from bottom to top
    % ----------------------------------------------------------

    [yProfile,idx] = sort(yProfile);
    uProfile = uProfile(idx);

    %% ----------------------------------------------------------
    % Convert to distance above hill
    % ----------------------------------------------------------

    nProfile = yProfile - yw;

    %% ----------------------------------------------------------
    % Keep only points above the wall
    % ----------------------------------------------------------

    valid = nProfile >= 0;

    nProfile = nProfile(valid);
    uProfile = uProfile(valid);

    %% ----------------------------------------------------------
    % Add wall point
    %
    % No-slip: U = 0 at the wall.
    % ----------------------------------------------------------

    nProfile = [0; nProfile];
    uProfile = [0; uProfile];

    %% ----------------------------------------------------------
    % Normalize
    % ----------------------------------------------------------

    Unorm = uProfile / Ub;

    %% ----------------------------------------------------------
    % Plot spatial representation
    % ----------------------------------------------------------

    xPlot = xh + profileScale*Unorm;
    yPlot = yw + nProfile;

    plot( ...
        xPlot, ...
        yPlot, ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('$x/h=%.1f$',xh));

    %% ----------------------------------------------------------
    % Vertical sampling line
    % ----------------------------------------------------------

    plot( ...
        [xh xh], ...
        [yw 3.035], ...
        ':', ...
        'Color',[0.6 0.6 0.6], ...
        'HandleVisibility','off');

    %% ----------------------------------------------------------
    % Sampling point on hill
    % ----------------------------------------------------------

    plot( ...
        xh, ...
        yw, ...
        'ko', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',4, ...
        'HandleVisibility','off');

end

%% ==============================================================
% FORMATTING
% ==============================================================

xlabel('$x/h$','Interpreter','latex');
ylabel('$y/h$','Interpreter','latex');

title( ...
    '$U/U_b$ velocity profiles over periodic hill', ...
    'Interpreter','latex');

xlim([0 9]);
ylim([0 3.035]);

grid on;
box on;

%legend( ...'Interpreter','latex', ...'Location','best');

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1);

%% ==============================================================
% FUNCTION: READ OPENFOAM VECTOR FIELD
% ==============================================================

function V = readOpenFOAMVectorField(filename)

    fid = fopen(filename,'r');

    if fid == -1
        error('Could not open %s',filename);
    end

    text = fread(fid,'*char')';
    fclose(fid);

    %% Find internalField

    idx = strfind(text,'internalField');

    if isempty(idx)
        error('Could not find internalField in %s',filename);
    end

    idx = idx(1);

    localText = text(idx:end);

    %% Check for nonuniform vector field

    token = 'nonuniform List<vector>';

    idx2 = strfind(localText,token);

    if isempty(idx2)

        error( ...
            'Expected nonuniform List<vector> in %s', ...
            filename);

    end

    idx2 = idx + idx2(1) - 1;

    %% Number of vectors

    afterType = text(idx2 + length(token):end);

    countMatch = regexp( ...
        afterType, ...
        '^\s*(\d+)', ...
        'tokens', ...
        'once');

    if isempty(countMatch)

        error( ...
            'Could not determine vector count in %s', ...
            filename);

    end

    N = str2double(countMatch{1});

    %% Opening parenthesis

    pRel = strfind(afterType,'(');

    if isempty(pRel)

        error('Could not find vector list in %s',filename);

    end

    p = idx2 + length(token) + pRel(1) - 1;

    %% Data after opening parenthesis

    dataText = text(p+1:end);

    %% Extract vectors

    expr = [ ...
        '\(\s*' ...
        '([-+0-9.eE]+)\s+' ...
        '([-+0-9.eE]+)\s+' ...
        '([-+0-9.eE]+)\s*\)' ];

    tokens = regexp( ...
        dataText, ...
        expr, ...
        'tokens');

    if length(tokens) < N

        error( ...
            'Expected %d vectors but only found %d in %s', ...
            N,length(tokens),filename);

    end

    %% Take ONLY internalField vectors

    tokens = tokens(1:N);

    V = zeros(N,3);

    for i = 1:N

        V(i,1) = str2double(tokens{i}{1});
        V(i,2) = str2double(tokens{i}{2});
        V(i,3) = str2double(tokens{i}{3});

    end

end

%% ==============================================================
% FUNCTION: HILL HEIGHT
% ==============================================================

function yw = hillHeight(xh)
% HILLHEIGHT
%
% Exact ERCOFTAC periodic-hill geometry supplied by the user.
%
% Input:
%     xh = x/h
%
% Output:
%     yw = y_wall/h

    %% Convert x/h -> original geometry coordinate in mm

    x_mm = mod(xh,9.0)*28.0;

    %% Mirror second half of the hill

    if x_mm > 198.0
        xs = 252.0 - x_mm;
    else
        xs = x_mm;
    end

    %% Piecewise hill definition

    if xs >= 0 && xs < 9

        y_mm = min( ...
            28.0, ...
            28.0 ...
          + 0.000000000000E+00*xs ...
          + 6.775070969851E-03*xs.^2 ...
          - 2.124527775800E-03*xs.^3);

    elseif xs >= 9 && xs < 14

        y_mm = ...
            2.507355893131E+01 ...
          + 9.754803562315E-01*xs ...
          - 1.016116352781E-01*xs.^2 ...
          + 1.889794677828E-03*xs.^3;

    elseif xs >= 14 && xs < 20

        y_mm = ...
            2.579601052357E+01 ...
          + 8.206693007457E-01*xs ...
          - 9.055370274339E-02*xs.^2 ...
          + 1.626510569859E-03*xs.^3;

    elseif xs >= 20 && xs < 30

        y_mm = ...
            4.046435022819E+01 ...
          - 1.379581654948E+00*xs ...
          + 1.945884504128E-02*xs.^2 ...
          - 2.070318932190E-04*xs.^3;

    elseif xs >= 30 && xs < 40

        y_mm = ...
            1.792461334664E+01 ...
          + 8.743920332081E-01*xs ...
          - 5.567361123058E-02*xs.^2 ...
          + 6.277731764683E-04*xs.^3;

    elseif xs >= 40 && xs < 54

        y_mm = max( ...
            0.0, ...
            5.639011190988E+01 ...
          - 2.010520359035E+00*xs ...
          + 1.644919857549E-02*xs.^2 ...
          + 2.674976141766E-05*xs.^3);

    else

        y_mm = 0;

    end

    %% Convert back to h-normalized coordinates

    yw = y_mm/28.0;

end