%% ==============================================================
% PERIODIC HILL: SPALDING vs ML
%
% Compares spanwise-averaged streamwise velocity profiles from
% multiple OpenFOAM simulations.
%
% The OpenFOAM coordinates are already normalized by hill height h.
%
% Produces:
%   1. Velocity profiles spatially represented over the hill
%   2. Conventional U/Ub vs n/h comparison plots at each x/h
%
% Current simulations:
%   - Spalding
%   - ML
%
% Mesh:
%   Nx = 120
%   Ny = 40
%   Nz = 60
%
% Bulk velocity:
%   Ub = 1
%
% ==============================================================

clear;
clc;
close all;

%% ==============================================================
% CASE SETTINGS
% ==============================================================
% endTime = 2000 (200 flow-through times of averaging + ~22 FTT spin-up, per
% controlDict) - was '0.1', an early, barely-past-the-initial-transient instant
% under the OLD (wrong) endTime=0.5 convention. Reads UMean (the fieldAverage'd
% field, not the instantaneous U - see readOpenFOAMVectorField(Ufile) below and
% system/functions' own comment on why an instantaneous LES snapshot isn't a fair
% comparison against Spalding's law, which describes the converged mean flow).

time = '2000';

% --------------------------------------------------------------
% Spalding simulation
% --------------------------------------------------------------

cases(1).name  = 'Spalding';
cases(1).dir   = 'spaldingSim';
cases(1).time  = time;
cases(1).color = [0.0000 0.4470 0.7410];

% --------------------------------------------------------------
% ML simulation
% --------------------------------------------------------------

cases(2).name  = 'ML';
cases(2).dir   = 'mlSim';
cases(2).time  = time;
cases(2).color = [0.8500 0.3250 0.0980];

%% ==============================================================
% MESH / FLOW SETTINGS
% ==============================================================

Nx = 120;
Ny = 40;
Nz = 60;

Ub = 1.0;

%% ==============================================================
% x/h LOCATIONS FOR VELOCITY PROFILES
% ==============================================================

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

%% ==============================================================
% PROFILE VISUALIZATION SETTINGS
% ==============================================================

% Horizontal distance used to represent U/Ub in the hill plot:
%
% x_plot = x/h + profileScale*(U/Ub)
%
profileScale = 0.5;

%% ==============================================================
% READ ALL SIMULATIONS
% ==============================================================

for c = 1:length(cases)

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('Reading %s\n',cases(c).name);
    fprintf('============================================\n');

    %% ----------------------------------------------------------
    % File names
    % ----------------------------------------------------------

    Cfile = fullfile( ...
        cases(c).dir, ...
        cases(c).time, ...
        'C');

    Ufile = fullfile( ...
        cases(c).dir, ...
        cases(c).time, ...
        'UMean');

    %% ----------------------------------------------------------
    % Read OpenFOAM fields
    % ----------------------------------------------------------

    cases(c).C = readOpenFOAMVectorField(Cfile);
    cases(c).U = readOpenFOAMVectorField(Ufile);

    fprintf( ...
        'C vectors = %d\n', ...
        size(cases(c).C,1));

    fprintf( ...
        'U vectors = %d\n', ...
        size(cases(c).U,1));

    %% ----------------------------------------------------------
    % Check mesh size
    % ----------------------------------------------------------

    expectedCells = Nx*Ny*Nz;

    if size(cases(c).C,1) ~= expectedCells

        error( ...
            '%s: C does not contain %d cells.', ...
            cases(c).name, ...
            expectedCells);

    end

    if size(cases(c).U,1) ~= expectedCells

        error( ...
            '%s: U does not contain %d values.', ...
            cases(c).name, ...
            expectedCells);

    end

    %% ----------------------------------------------------------
    % Reshape fields
    %
    % Based on your mesh ordering:
    %
    %   x varies first
    %   y varies second
    %   z varies third
    %
    % Therefore:
    %
    %   [Nx Ny Nz]
    %
    % ----------------------------------------------------------

    cases(c).X = reshape( ...
        cases(c).C(:,1), ...
        [Nx Ny Nz]);

    cases(c).Y = reshape( ...
        cases(c).C(:,2), ...
        [Nx Ny Nz]);

    cases(c).Z = reshape( ...
        cases(c).C(:,3), ...
        [Nx Ny Nz]);

    cases(c).UX = reshape( ...
        cases(c).U(:,1), ...
        [Nx Ny Nz]);

    cases(c).UY = reshape( ...
        cases(c).U(:,2), ...
        [Nx Ny Nz]);

    cases(c).UZ = reshape( ...
        cases(c).U(:,3), ...
        [Nx Ny Nz]);

end

%% ==============================================================
% FIGURE 1
%
% VELOCITY PROFILES SPATIALLY REPRESENTED OVER THE HILL
% ==============================================================

figure('Color','w');
hold on;

%% --------------------------------------------------------------
% Plot exact hill geometry
% --------------------------------------------------------------

xhHill = linspace(0,9,1000);

ywHill = arrayfun( ...
    @hillHeight, ...
    xhHill);

plot( ...
    xhHill, ...
    ywHill, ...
    'k-', ...
    'LineWidth',2.0, ...
    'DisplayName','Hill');

%% --------------------------------------------------------------
% Plot top wall
% --------------------------------------------------------------

plot( ...
    [0 9], ...
    [3.035 3.035], ...
    'k--', ...
    'LineWidth',1.0, ...
    'HandleVisibility','off');

%% --------------------------------------------------------------
% Plot velocity profiles
% --------------------------------------------------------------

legendHandles = gobjects(length(cases),1);

for k = 1:length(xStations)

    xh = xStations(k);

    %% ----------------------------------------------------------
    % Hill height at sampling station
    % ----------------------------------------------------------

    yw = hillHeight(xh);

    %% ----------------------------------------------------------
    % Vertical sampling line
    % ----------------------------------------------------------

    plot( ...
        [xh xh], ...
        [yw 3.035], ...
        ':', ...
        'Color',[0.7 0.7 0.7], ...
        'HandleVisibility','off');

    %% ----------------------------------------------------------
    % Extract each simulation profile
    % ----------------------------------------------------------

    for c = 1:length(cases)

        yProfile = zeros(Ny,1);
        uProfile = zeros(Ny,1);

        %% ------------------------------------------------------
        % Extract profile at x/h
        % -------------------------------------------------------

        for j = 1:Ny

            %% Coordinates in this y layer

            Xlayer = squeeze( ...
                cases(c).X(:,j,:));

            %% Average x over span

            Xmean = mean(Xlayer,2);

            %% Find nearest x index

            [~,i] = min(abs(Xmean-xh));

            %% Velocity at this x/y over all z

            ULayer = squeeze( ...
                cases(c).UX(i,j,:));

            %% Coordinates over z

            YLayer = squeeze( ...
                cases(c).Y(i,j,:));

            %% Spanwise average

            uProfile(j) = mean(ULayer);
            yProfile(j) = mean(YLayer);

        end

        %% ------------------------------------------------------
        % Sort from bottom to top
        % -------------------------------------------------------

        [yProfile,idx] = sort(yProfile);

        uProfile = uProfile(idx);

        %% ------------------------------------------------------
        % Distance above hill
        % -------------------------------------------------------

        n = yProfile - yw;

        %% ------------------------------------------------------
        % Remove points below the wall
        % -------------------------------------------------------

        valid = n >= 0;

        n = n(valid);
        uProfile = uProfile(valid);

        %% ------------------------------------------------------
        % Add no-slip wall point
        % -------------------------------------------------------

        n = [0; n];
        uProfile = [0; uProfile];

        %% ------------------------------------------------------
        % Normalize velocity
        % -------------------------------------------------------

        U_norm = uProfile / Ub;

        %% ------------------------------------------------------
        % Convert into spatial representation
        %
        % xPlot = station + velocity displacement
        % yPlot = wall + wall-normal distance
        % -------------------------------------------------------

        xPlot = xh + profileScale*U_norm;

        yPlot = yw + n;

        %% ------------------------------------------------------
        % Plot
        % -------------------------------------------------------

        if k == 1

            legendHandles(c) = plot( ...
                xPlot, ...
                yPlot, ...
                'Color',cases(c).color, ...
                'LineWidth',1.8, ...
                'DisplayName',cases(c).name);

        else

            plot( ...
                xPlot, ...
                yPlot, ...
                'Color',cases(c).color, ...
                'LineWidth',1.8, ...
                'HandleVisibility','off');

        end

    end

end

%% --------------------------------------------------------------
% Formatting
% --------------------------------------------------------------

xlabel('$x/h$', ...
    'Interpreter','latex');

ylabel('$y/h$', ...
    'Interpreter','latex');

title( ...
    '$U/U_b$ profiles over periodic hill', ...
    'Interpreter','latex');

legend( ...
    legendHandles, ...
    {cases.name}, ...
    'Location','best');

xlim([0 9]);
ylim([0 3.035]);

grid on;
box on;

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1);

%% ==============================================================
% FIGURES 2+
%
% CONVENTIONAL VELOCITY COMPARISON
%
% U/Ub vs n/h
% ==============================================================

for k = 1:length(xStations)

    xh = xStations(k);

    %% ----------------------------------------------------------
    % New figure
    % ----------------------------------------------------------

    figure('Color','w');
    hold on;

    %% ----------------------------------------------------------
    % Hill height
    % ----------------------------------------------------------

    yw = hillHeight(xh);

    %% ----------------------------------------------------------
    % Plot each simulation
    % ----------------------------------------------------------

    for c = 1:length(cases)

        yProfile = zeros(Ny,1);
        uProfile = zeros(Ny,1);

        %% ------------------------------------------------------
        % Extract profile
        % -------------------------------------------------------

        for j = 1:Ny

            Xlayer = squeeze( ...
                cases(c).X(:,j,:));

            Xmean = mean(Xlayer,2);

            [~,i] = min(abs(Xmean-xh));

            ULayer = squeeze( ...
                cases(c).UX(i,j,:));

            YLayer = squeeze( ...
                cases(c).Y(i,j,:));

            uProfile(j) = mean(ULayer);
            yProfile(j) = mean(YLayer);

        end

        %% ------------------------------------------------------
        % Sort
        % -------------------------------------------------------

        [yProfile,idx] = sort(yProfile);

        uProfile = uProfile(idx);

        %% ------------------------------------------------------
        % Wall-normal distance
        % -------------------------------------------------------

        n = yProfile - yw;

        %% ------------------------------------------------------
        % Remove points below wall
        % -------------------------------------------------------

        valid = n >= 0;

        n = n(valid);
        uProfile = uProfile(valid);

        %% ------------------------------------------------------
        % Add wall point
        % -------------------------------------------------------

        n = [0; n];
        uProfile = [0; uProfile];

        %% ------------------------------------------------------
        % Normalize
        % -------------------------------------------------------

        U_norm = uProfile / Ub;

        %% ------------------------------------------------------
        % Plot
        % -------------------------------------------------------

        plot( ...
            U_norm, ...
            n, ...
            'Color',cases(c).color, ...
            'LineWidth',1.8, ...
            'DisplayName',cases(c).name);

    end

    %% ----------------------------------------------------------
    % Formatting
    % ----------------------------------------------------------

    xlabel('$U/U_b$', ...
        'Interpreter','latex');

    ylabel('$n/h$', ...
        'Interpreter','latex');

    title( ...
        sprintf( ...
            'Velocity comparison at $x/h=%.1f$', ...
            xh), ...
        'Interpreter','latex');

    legend( ...
        'Location','best');

    grid on;
    box on;

    set(gca, ...
        'FontSize',11, ...
        'LineWidth',1);

end

%% ==============================================================
% FUNCTION: READ OPENFOAM VECTOR FIELD
% ==============================================================

function V = readOpenFOAMVectorField(filename)

    %% ----------------------------------------------------------
    % Open file
    % ----------------------------------------------------------

    fid = fopen(filename,'r');

    if fid == -1

        error( ...
            'Could not open %s', ...
            filename);

    end

    %% ----------------------------------------------------------
    % Read entire file
    % ----------------------------------------------------------

    text = fread(fid,'*char')';

    fclose(fid);

    %% ----------------------------------------------------------
    % Locate internalField
    % ----------------------------------------------------------

    idx = strfind( ...
        text, ...
        'internalField');

    if isempty(idx)

        error( ...
            'Could not find internalField in %s', ...
            filename);

    end

    idx = idx(1);

    localText = text(idx:end);

    %% ----------------------------------------------------------
    % Find nonuniform vector field
    % ----------------------------------------------------------

    token = 'nonuniform List<vector>';

    idx2 = strfind( ...
        localText, ...
        token);

    if isempty(idx2)

        error( ...
            'Could not find nonuniform List<vector> in %s', ...
            filename);

    end

    idx2 = idx + idx2(1) - 1;

    %% ----------------------------------------------------------
    % Find number of entries
    % ----------------------------------------------------------

    afterType = ...
        text(idx2 + length(token):end);

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

    N = str2double( ...
        countMatch{1});

    %% ----------------------------------------------------------
    % Find opening parenthesis
    % ----------------------------------------------------------

    pRel = strfind( ...
        afterType, ...
        '(');

    if isempty(pRel)

        error( ...
            'Could not find vector list in %s', ...
            filename);

    end

    p = ...
        idx2 + length(token) + pRel(1) - 1;

    %% ----------------------------------------------------------
    % Extract text after opening parenthesis
    % ----------------------------------------------------------

    dataText = text(p+1:end);

    %% ----------------------------------------------------------
    % Vector regex
    % ----------------------------------------------------------

    expr = [ ...
        '\(\s*' ...
        '([-+0-9.eE]+)\s+' ...
        '([-+0-9.eE]+)\s+' ...
        '([-+0-9.eE]+)\s*\)' ];

    tokens = regexp( ...
        dataText, ...
        expr, ...
        'tokens');

    %% ----------------------------------------------------------
    % Check number of vectors
    % ----------------------------------------------------------

    if length(tokens) < N

        error( ...
            ['Expected %d vectors but found %d ' ...
             'in %s'], ...
            N, ...
            length(tokens), ...
            filename);

    end

    %% ----------------------------------------------------------
    % Keep only internalField
    % ----------------------------------------------------------

    tokens = tokens(1:N);

    %% ----------------------------------------------------------
    % Convert strings to numbers
    % ----------------------------------------------------------

    V = zeros(N,3);

    for i = 1:N

        V(i,1) = ...
            str2double(tokens{i}{1});

        V(i,2) = ...
            str2double(tokens{i}{2});

        V(i,3) = ...
            str2double(tokens{i}{3});

    end

end

%% ==============================================================
% FUNCTION: EXACT PERIODIC HILL GEOMETRY
% ==============================================================

function yw = hillHeight(xh)
% HILLHEIGHT
%
% Exact hill geometry supplied in the OpenFOAM blockMeshDict.
%
% Input:
%     xh = x/h
%
% Output:
%     yw = y_wall/h
%
% Period:
%     9h
%
% Hill height:
%     h
%
% Original polynomial coordinates:
%     mm
% ==============================================================

    %% ----------------------------------------------------------
    % Convert x/h to original mm coordinate
    % ----------------------------------------------------------

    x_mm = mod(xh,9.0)*28.0;

    %% ----------------------------------------------------------
    % Mirror second half
    % ----------------------------------------------------------

    if x_mm > 198.0

        xs = 252.0 - x_mm;

    else

        xs = x_mm;

    end

    %% ----------------------------------------------------------
    % Piecewise hill polynomial
    % ----------------------------------------------------------

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

    %% ----------------------------------------------------------
    % Normalize by hill height
    % ----------------------------------------------------------

    yw = y_mm/28.0;

end