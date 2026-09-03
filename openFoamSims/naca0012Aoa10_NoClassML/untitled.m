%% ============================================================
% NACA 0012 LES: AoA = 10 deg
%
% Compare:
%   1. LES Cp against Ladson et al.
%   2. Cp difference LES - Ladson
%   3. Cumulative pressure drag
%   4. LES CL against Gregory
%
% Files:
%   naca0012LES_aoa10/comparison_results/Cp.csv
%   naca0012LES_aoa10/comparison_results/summary.csv
%   exp.dat
%   Gregory_CL_expdata.dat
%
% Cp.csv:
%   x/c, z/c, Cp
%
% ============================================================

clear;
close all;
clc;

%% ============================================================
% SETTINGS
% ============================================================

AoA = 10;

caseDir = ...
    'fine_naca0012LES_aoa10/comparison_results';

cpFile = fullfile(caseDir,'Cp.csv');

summaryFile = ...
    fullfile(caseDir,'summary.csv');

expFile = 'exp.dat';

gregoryFile = ...
    'Gregory_CL_expdata.dat';

%% ============================================================
% CHECK FILES
% ============================================================

requiredFiles = {
    cpFile
    summaryFile
    expFile
    gregoryFile
};

for i = 1:length(requiredFiles)

    if ~isfile(requiredFiles{i})

        error( ...
            'Could not find file:\n%s', ...
            requiredFiles{i});

    end

end

%% ============================================================
% READ LES Cp
% ============================================================

cpData = readmatrix(cpFile);

cpData = cpData(all(isfinite(cpData),2),:);

if size(cpData,2) < 3

    error( ...
        '%s must contain x/c, z/c, Cp.', ...
        cpFile);

end

LES_x  = cpData(:,1);
LES_z  = cpData(:,2);
LES_Cp = cpData(:,3);

%% ============================================================
% SPLIT LES UPPER / LOWER
% ============================================================

upperMask = LES_z >= 0;

lowerMask = LES_z < 0;

LES_xu = LES_x(upperMask);
LES_cpu = LES_Cp(upperMask);

LES_xl = LES_x(lowerMask);
LES_cpl = LES_Cp(lowerMask);

%% ============================================================
% SORT LES
% ============================================================

[LES_xu,idx] = sort(LES_xu);

LES_cpu = LES_cpu(idx);

[LES_xl,idx] = sort(LES_xl);

LES_cpl = LES_cpl(idx);

%% ============================================================
% REMOVE DUPLICATE x/c
% ============================================================

[LES_xu,idx] = unique( ...
    LES_xu, ...
    'stable');

LES_cpu = LES_cpu(idx);

[LES_xl,idx] = unique( ...
    LES_xl, ...
    'stable');

LES_cpl = LES_cpl(idx);

fprintf('\nLES Cp data:\n');

fprintf( ...
    '  Total points = %d\n', ...
    length(LES_x));

fprintf( ...
    '  Upper points = %d\n', ...
    length(LES_xu));

fprintf( ...
    '  Lower points = %d\n', ...
    length(LES_xl));

%% ============================================================
% READ LES SUMMARY
% ============================================================

summary = readtable(summaryFile);

idx = abs(summary.AoA_deg - AoA) < 1e-8;

if ~any(idx)

    error( ...
        'No result found in summary.csv for AoA = %.1f deg.', ...
        AoA);

end

row = summary(find(idx,1),:);

LES_CL  = row.Cl;

LES_CD  = row.Cd;

LES_CLp = row.Cl_pressure;

LES_CLf = row.Cl_viscous;

LES_CDp = row.Cd_pressure;

LES_CDf = row.Cd_viscous;

%% ============================================================
% READ LADSON Cp DATA
% ============================================================

fid = fopen(expFile,'rt');

if fid == -1

    error( ...
        'Could not open experimental file: %s', ...
        expFile);

end

expData = struct( ...
    'alpha',{}, ...
    'x',{}, ...
    'cp',{});

currentAlpha = NaN;

xData = [];

cpDataExp = [];

while true

    line = fgetl(fid);

    if isequal(line,-1)

        break;

    end

    line = char(line);

    line = strtrim(line);

    if isempty(line)

        continue;

    end

    % --------------------------------------------------------
    % Zone
    % --------------------------------------------------------

    if length(line) >= 4 && ...
            strncmpi(line,'zone',4)

        if ~isempty(xData) && ...
                ~isnan(currentAlpha)

            expData(end+1).alpha = ...
                currentAlpha;

            expData(end).x = ...
                xData;

            expData(end).cp = ...
                cpDataExp;

        end

        token = regexp( ...
            line, ...
            'alpha\s*=\s*([-+]?\d*\.?\d+)', ...
            'tokens', ...
            'once');

        currentAlpha = NaN;

        if ~isempty(token)

            currentAlpha = ...
                str2double(token{1});

        end

        xData = [];

        cpDataExp = [];

        continue;

    end

    % --------------------------------------------------------
    % Comments
    % --------------------------------------------------------

    if line(1) == '#'

        continue;

    end

    % --------------------------------------------------------
    % Variables line
    % --------------------------------------------------------

    if length(line) >= 9 && ...
            strncmpi(line,'variables',9)

        continue;

    end

    % --------------------------------------------------------
    % Numeric x/c, Cp
    % --------------------------------------------------------

    vals = sscanf( ...
        line, ...
        '%f %f');

    if numel(vals) >= 2

        xData(end+1,1) = vals(1);

        cpDataExp(end+1,1) = vals(2);

    end

end

fclose(fid);

%% ============================================================
% SAVE FINAL LADSON ZONE
% ============================================================

if ~isempty(xData) && ...
        ~isnan(currentAlpha)

    expData(end+1).alpha = ...
        currentAlpha;

    expData(end).x = ...
        xData;

    expData(end).cp = ...
        cpDataExp;

end

%% ============================================================
% SELECT LADSON AoA
% ============================================================

if isempty(expData)

    error( ...
        'No Ladson experimental zones were found.');

end

[minDiff,idx] = ...
    min(abs([expData.alpha] - AoA));

if minDiff > 0.5

    error( ...
        'No Ladson Cp data within 0.5 deg of AoA = %.1f deg.', ...
        AoA);

end

Ladson = expData(idx);

fprintf('\nLadson data:\n');

fprintf( ...
    '  Requested AoA = %.1f deg\n', ...
    AoA);

fprintf( ...
    '  Actual AoA    = %.4f deg\n', ...
    Ladson.alpha);

fprintf( ...
    '  Points        = %d\n', ...
    length(Ladson.x));

%% ============================================================
% SPLIT LADSON
%
% Data ordering:
%
%   upper -> leading edge -> lower
%
% ============================================================

%% ============================================================
% SPLIT LADSON INTO TWO SURFACE BRANCHES
% ============================================================

[~,iLE] = min(abs(Ladson.x));

branch1_x  = Ladson.x(1:iLE);
branch1_cp = Ladson.cp(1:iLE);

branch2_x  = Ladson.x(iLE:end);
branch2_cp = Ladson.cp(iLE:end);

% Sort both branches
[branch1_x,idx] = sort(branch1_x);
branch1_cp = branch1_cp(idx);

[branch2_x,idx] = sort(branch2_x);
branch2_cp = branch2_cp(idx);

% Remove duplicate x/c values
[branch1_x,idx] = unique(branch1_x,'stable');
branch1_cp = branch1_cp(idx);

[branch2_x,idx] = unique(branch2_x,'stable');
branch2_cp = branch2_cp(idx);

%% ============================================================
% IDENTIFY UPPER / LOWER AUTOMATICALLY
%
% At positive AoA, the upper surface should have the more
% negative Cp around mid-chord.
% ============================================================

xTest = 0.5;

cp1Test = interp1( ...
    branch1_x, ...
    branch1_cp, ...
    xTest, ...
    'linear', ...
    NaN);

cp2Test = interp1( ...
    branch2_x, ...
    branch2_cp, ...
    xTest, ...
    'linear', ...
    NaN);

fprintf('\nLadson branch check at x/c = %.2f:\n',xTest);
fprintf('  Branch 1 Cp = %.6f\n',cp1Test);
fprintf('  Branch 2 Cp = %.6f\n',cp2Test);

if cp1Test < cp2Test

    % Branch 1 is more negative -> upper
    Ladson_xu  = branch1_x;
    Ladson_cpu = branch1_cp;

    Ladson_xl  = branch2_x;
    Ladson_cpl = branch2_cp;

    fprintf('  Branch 1 identified as UPPER.\n');
    fprintf('  Branch 2 identified as LOWER.\n');

else

    % Branch 2 is more negative -> upper
    Ladson_xu  = branch2_x;
    Ladson_cpu = branch2_cp;

    Ladson_xl  = branch1_x;
    Ladson_cpl = branch1_cp;

    fprintf('  Branch 2 identified as UPPER.\n');
    fprintf('  Branch 1 identified as LOWER.\n');

end

%% ============================================================
% SORT LADSON
% ============================================================

[Ladson_xu,idx] = ...
    sort(Ladson_xu);

Ladson_cpu = ...
    Ladson_cpu(idx);

[Ladson_xl,idx] = ...
    sort(Ladson_xl);

Ladson_cpl = ...
    Ladson_cpl(idx);

%% ============================================================
% REMOVE DUPLICATES
% ============================================================

[Ladson_xu,idx] = ...
    unique(Ladson_xu,'stable');

Ladson_cpu = ...
    Ladson_cpu(idx);

[Ladson_xl,idx] = ...
    unique(Ladson_xl,'stable');

Ladson_cpl = ...
    Ladson_cpl(idx);

fprintf('\nLadson surfaces:\n');

fprintf( ...
    '  Upper = %d points\n', ...
    length(Ladson_xu));

fprintf( ...
    '  Lower = %d points\n', ...
    length(Ladson_xl));

%% ============================================================
% INTERPOLATE LES Cp ONTO LADSON
% ============================================================

LES_Cp_u_interp = interp1( ...
    LES_xu, ...
    LES_cpu, ...
    Ladson_xu, ...
    'linear', ...
    NaN);

LES_Cp_l_interp = interp1( ...
    LES_xl, ...
    LES_cpl, ...
    Ladson_xl, ...
    'linear', ...
    NaN);

%% ============================================================
% Cp ERRORS
% ============================================================

validU = ...
    isfinite(LES_Cp_u_interp) & ...
    isfinite(Ladson_cpu);

validL = ...
    isfinite(LES_Cp_l_interp) & ...
    isfinite(Ladson_cpl);

Cp_error_upper = ...
    LES_Cp_u_interp(validU) - ...
    Ladson_cpu(validU);

Cp_error_lower = ...
    LES_Cp_l_interp(validL) - ...
    Ladson_cpl(validL);

Cp_upper_RMSE = ...
    sqrt(mean(Cp_error_upper.^2));

Cp_lower_RMSE = ...
    sqrt(mean(Cp_error_lower.^2));

Cp_all_error = [
    Cp_error_upper
    Cp_error_lower
];

Cp_RMSE = ...
    sqrt(mean(Cp_all_error.^2));

Cp_MAE = ...
    mean(abs(Cp_all_error));

%% ============================================================
% READ GREGORY CL
% ============================================================

fid = fopen(gregoryFile,'rt');

if fid == -1

    error( ...
        'Could not open Gregory file: %s', ...
        gregoryFile);

end

gregoryAlpha = [];

gregoryCL = [];

while ~feof(fid)

    line = fgetl(fid);

    if ~ischar(line)

        continue;

    end

    line = strtrim(line);

    if isempty(line)

        continue;

    end

    if line(1) == '#'

        continue;

    end

    if length(line) >= 9 && ...
            strncmpi(line,'variables',9)

        continue;

    end

    if length(line) >= 4 && ...
            strncmpi(line,'zone',4)

        continue;

    end

    vals = sscanf( ...
        line, ...
        '%f %f');

    if numel(vals) >= 2

        gregoryAlpha(end+1,1) = ...
            vals(1);

        gregoryCL(end+1,1) = ...
            vals(2);

    end

end

fclose(fid);

valid = ...
    isfinite(gregoryAlpha) & ...
    isfinite(gregoryCL);

gregoryAlpha = ...
    gregoryAlpha(valid);

gregoryCL = ...
    gregoryCL(valid);

%% ============================================================
% GREGORY CL AT 10 DEG
% ============================================================

Gregory_CL = interp1( ...
    gregoryAlpha, ...
    gregoryCL, ...
    AoA, ...
    'linear', ...
    'extrap');

CL_error = ...
    LES_CL - Gregory_CL;

CL_abs_error = ...
    abs(CL_error);

CL_percent_error = ...
    100 * CL_error / abs(Gregory_CL);

%% ============================================================
% DIRECT 2D PRESSURE FORCE
%
% Reconstruct the LES contour and integrate:
%
%   Cx = - integral(Cp dz)
%   Cz =   integral(Cp dx)
%
% ============================================================

% ------------------------------------------------------------
% Build closed LES contour
% ------------------------------------------------------------

xc = [
    LES_xu
    flipud(LES_xl)
];

zc = [
    getSurfaceZ(LES_xu,LES_xu,LES_xu) % temporary
];

% Replace with actual upper/lower z data reconstructed
% from the original LES data.
xuZ = LES_z(upperMask);
xlZ = LES_z(lowerMask);

[~,iu] = sort(LES_x(upperMask));
xuZ = xuZ(iu);

[~,il] = sort(LES_x(lowerMask));
xlZ = xlZ(il);

% Remove duplicate x locations consistently
[~,iu] = unique(LES_x(upperMask),'stable');
xuZ = xuZ(iu);

[~,il] = unique(LES_x(lowerMask),'stable');
xlZ = xlZ(il);

zc = [
    xuZ
    flipud(xlZ)
];

cpc = [
    LES_cpu
    flipud(LES_cpl)
];

xc = [
    LES_xu
    flipud(LES_xl)
];

% Remove duplicate consecutive points
ds = sqrt( ...
    diff(xc).^2 + ...
    diff(zc).^2);

keep = [
    true
    ds > 1e-12
];

xc = xc(keep);

zc = zc(keep);

cpc = cpc(keep);

% Close contour
if xc(end) ~= xc(1) || ...
        zc(end) ~= zc(1)

    xc(end+1) = xc(1);

    zc(end+1) = zc(1);

    cpc(end+1) = cpc(1);

end

%% ============================================================
% ENFORCE COUNTER-CLOCKWISE CONTOUR
% ============================================================

areaSigned = ...
    0.5 * sum( ...
        xc(1:end-1).*zc(2:end) - ...
        xc(2:end).*zc(1:end-1));

if areaSigned < 0

    xc = flipud(xc);

    zc = flipud(zc);

    cpc = flipud(cpc);

end

%% ============================================================
% LOCAL PRESSURE FORCE
% ============================================================

dx = diff(xc);

dz = diff(zc);

CpMid = ...
    0.5 * ...
    (cpc(1:end-1) + ...
     cpc(2:end));

dCx = ...
    -CpMid .* dz;

dCz = ...
     CpMid .* dx;

alphaRad = ...
    deg2rad(AoA);

dCd = ...
    dCx*cos(alphaRad) + ...
    dCz*sin(alphaRad);

dCl = ...
    -dCx*sin(alphaRad) + ...
     dCz*cos(alphaRad);

%% ============================================================
% CUMULATIVE DRAG
%
% The contour starts at the leading edge and follows
% the airfoil surface. We calculate cumulative force
% along the contour.
% ============================================================

cumCd = cumsum(dCd);

cumCl = cumsum(dCl);

xMid = ...
    0.5 * ...
    (xc(1:end-1) + xc(2:end));

%% ============================================================
% SORT CUMULATIVE CONTRIBUTION BY STREAMWISE x/c
%
% The contour goes:
%
%   LE -> upper -> TE -> lower -> LE
%
% Therefore x/c alone is not monotonic.
%
% We separately calculate upper and lower cumulative
% drag from LE to TE.
% ============================================================

% ------------------------------------------------------------
% Upper surface
% ------------------------------------------------------------

upperContour = ...
    zc(1:end-1) >= 0;

xuSeg = xMid(upperContour);

dCdUpper = dCd(upperContour);

[~,order] = sort(xuSeg);

xuSeg = xuSeg(order);

dCdUpper = dCdUpper(order);

CdUpperCum = cumsum(dCdUpper);

% ------------------------------------------------------------
% Lower surface
% ------------------------------------------------------------

lowerContour = ...
    zc(1:end-1) < 0;

xlSeg = xMid(lowerContour);

dCdLower = dCd(lowerContour);

[~,order] = sort(xlSeg);

xlSeg = xlSeg(order);

dCdLower = dCdLower(order);

CdLowerCum = cumsum(dCdLower);

%% ============================================================
% INTERPOLATE CUMULATIVE CONTRIBUTIONS TO COMMON x GRID
% ============================================================

xGrid = linspace( ...
    0, ...
    1, ...
    500)';

CdUpperGrid = interp1( ...
    xuSeg, ...
    CdUpperCum, ...
    xGrid, ...
    'linear', ...
    'extrap');

CdLowerGrid = interp1( ...
    xlSeg, ...
    CdLowerCum, ...
    xGrid, ...
    'linear', ...
    'extrap');

CdTotalGrid = ...
    CdUpperGrid + ...
    CdLowerGrid;

%% ============================================================
% PRINT DIRECT INTEGRATION
% ============================================================

CdDirect = ...
    sum(dCd);

ClDirect = ...
    sum(dCl);

fprintf('\n');
fprintf('============================================================\n');
fprintf('PRESSURE FORCE CHECK\n');
fprintf('============================================================\n');

fprintf( ...
    'Direct 2D Cd,p = %.8f\n', ...
    CdDirect);

fprintf( ...
    'Direct 2D Cl,p = %.8f\n', ...
    ClDirect);

fprintf( ...
    '3D Cd,p        = %.8f\n', ...
    LES_CDp);

fprintf( ...
    '3D Cl,p        = %.8f\n', ...
    LES_CLp);

fprintf( ...
    'Cd difference  = %.8f\n', ...
    CdDirect - LES_CDp);

%% ============================================================
% ONE WINDOW WITH TABS
% ============================================================

plotWindow = figure( ...
    'Name','NACA 0012 LES AoA = 10 deg', ...
    'NumberTitle','off', ...
    'Color','w', ...
    'Position',[100 100 1100 750]);

tabGroup = uitabgroup(plotWindow);

fontName = 'Times New Roman';
axisFontSize = 9;
labelFontSize = 10;
titleFontSize = 11;

%% ============================================================
% TAB 1: Cp
% ============================================================

tab = uitab( ...
    tabGroup, ...
    'Title','Cp');

ax = axes('Parent',tab);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

plot( ...
    ax, ...
    Ladson_xu, ...
    Ladson_cpu, ...
    'ko', ...
    'MarkerSize',4, ...
    'LineStyle','none', ...
    'DisplayName','Ladson et al.');

plot( ...
    ax, ...
    Ladson_xl, ...
    Ladson_cpl, ...
    'ko', ...
    'MarkerSize',4, ...
    'LineStyle','none', ...
    'HandleVisibility','off');

plot( ...
    ax, ...
    LES_xu, ...
    LES_cpu, ...
    'b-', ...
    'LineWidth',1.5, ...
    'DisplayName','LES');

plot( ...
    ax, ...
    LES_xl, ...
    LES_cpl, ...
    'b-', ...
    'LineWidth',1.5, ...
    'HandleVisibility','off');

set(ax,'YDir','reverse');

xlim(ax,[0 1]);

xlabel(ax,'$x/c$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

ylabel(ax,'$C_p$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

title(ax, ...
    sprintf('NACA 0012: $C_p$, $\\alpha=%g^\\circ$',AoA), ...
    'Interpreter','latex', ...
    'FontSize',titleFontSize);

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'LineWidth',0.7, ...
    'TickDir','out');

%% ============================================================
% TAB 2: Cp ERROR
% ============================================================

tab = uitab( ...
    tabGroup, ...
    'Title','Cp Error');

ax = axes('Parent',tab);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

plot( ...
    ax, ...
    Ladson_xu(validU), ...
    Cp_error_upper, ...
    'b-', ...
    'LineWidth',1.5, ...
    'DisplayName','Upper');

plot( ...
    ax, ...
    Ladson_xl(validL), ...
    Cp_error_lower, ...
    'r-', ...
    'LineWidth',1.5, ...
    'DisplayName','Lower');

yline( ...
    ax, ...
    0, ...
    ':k', ...
    'HandleVisibility','off');

xlim(ax,[0 1]);

xlabel(ax,'$x/c$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

ylabel(ax,'$\Delta C_p$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

title(ax, ...
    '$C_p$ Error: LES - Ladson', ...
    'Interpreter','latex', ...
    'FontSize',titleFontSize);

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'LineWidth',0.7, ...
    'TickDir','out');

%% ============================================================
% TAB 3: CUMULATIVE Cd,p
% ============================================================

tab = uitab( ...
    tabGroup, ...
    'Title','Cumulative Cd,p');

ax = axes('Parent',tab);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

plot( ...
    ax, ...
    xGrid, ...
    CdUpperGrid, ...
    'b-', ...
    'LineWidth',1.5, ...
    'DisplayName','Upper');

plot( ...
    ax, ...
    xGrid, ...
    CdLowerGrid, ...
    'r-', ...
    'LineWidth',1.5, ...
    'DisplayName','Lower');

plot( ...
    ax, ...
    xGrid, ...
    CdTotalGrid, ...
    'k-', ...
    'LineWidth',2.0, ...
    'DisplayName','Total');

yline( ...
    ax, ...
    LES_CDp, ...
    '--k', ...
    'LineWidth',1.0, ...
    'DisplayName','$C_{D,p}$');

yline( ...
    ax, ...
    0, ...
    ':k', ...
    'HandleVisibility','off');

xlim(ax,[0 1]);

xlabel(ax,'$x/c$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

ylabel(ax,'Cumulative $C_{D,p}$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

title(ax, ...
    sprintf('Cumulative Pressure Drag, $\\alpha=%g^\\circ$',AoA), ...
    'Interpreter','latex', ...
    'FontSize',titleFontSize);

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'LineWidth',0.7, ...
    'TickDir','out');

%% ============================================================
% TAB 4: CL
% ============================================================

tab = uitab( ...
    tabGroup, ...
    'Title','CL');

ax = axes('Parent',tab);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

plot( ...
    ax, ...
    gregoryAlpha, ...
    gregoryCL, ...
    'ko-', ...
    'LineWidth',0.8, ...
    'MarkerSize',4, ...
    'DisplayName','Gregory');

plot( ...
    ax, ...
    AoA, ...
    LES_CL, ...
    'bo', ...
    'LineWidth',1.4, ...
    'MarkerSize',7, ...
    'DisplayName','LES');

plot( ...
    ax, ...
    AoA, ...
    Gregory_CL, ...
    'ks', ...
    'LineWidth',1.0, ...
    'MarkerSize',6, ...
    'DisplayName','Gregory at 10$^\circ$');

xline( ...
    ax, ...
    AoA, ...
    ':k', ...
    'HandleVisibility','off');

xlabel(ax,'$\alpha$ [deg]', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

ylabel(ax,'$C_L$', ...
    'Interpreter','latex', ...
    'FontSize',labelFontSize);

title(ax, ...
    '$C_L$ Comparison', ...
    'Interpreter','latex', ...
    'FontSize',titleFontSize);

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName',fontName, ...
    'FontSize',axisFontSize, ...
    'LineWidth',0.7, ...
    'TickDir','out');

%% ============================================================
% PRINT SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('TABBED COMPARISON WINDOW CREATED\n');
fprintf('============================================================\n');

fprintf('\nTabs:\n');
fprintf('  1. Cp\n');
fprintf('  2. Cp Error\n');
fprintf('  3. Cumulative Cd,p\n');
fprintf('  4. CL\n');
fprintf('\n');
%% ============================================================
% FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nCL:\n');

fprintf( ...
    '  Gregory = %.6f\n', ...
    Gregory_CL);

fprintf( ...
    '  LES     = %.6f\n', ...
    LES_CL);

fprintf( ...
    '  Error   = %.6f\n', ...
    CL_error);

fprintf( ...
    '  %% Error = %.3f %%\n', ...
    CL_percent_error);

fprintf('\nPressure drag:\n');

fprintf( ...
    '  LES 3D integration = %.8f\n', ...
    LES_CDp);

fprintf( ...
    '  Direct 2D Cp       = %.8f\n', ...
    CdDirect);

fprintf( ...
    '  Difference         = %.8f\n', ...
    CdDirect - LES_CDp);

fprintf('\nCp errors:\n');

fprintf( ...
    '  Upper RMSE = %.6f\n', ...
    Cp_upper_RMSE);

fprintf( ...
    '  Lower RMSE = %.6f\n', ...
    Cp_lower_RMSE);

fprintf( ...
    '  Total RMSE = %.6f\n', ...
    Cp_RMSE);

fprintf('\nLES drag:\n');

fprintf( ...
    '  Cd total    = %.8f\n', ...
    LES_CD);

fprintf( ...
    '  Cd pressure = %.8f\n', ...
    LES_CDp);

fprintf( ...
    '  Cd viscous  = %.8f\n', ...
    LES_CDf);

fprintf('\nGenerated figures:\n');

fprintf( ...
    '  %s\n', ...
    cpPdf);

fprintf( ...
    '  %s\n', ...
    cpErrorPdf);

fprintf( ...
    '  %s\n', ...
    dragPdf);

fprintf( ...
    '  %s\n', ...
    clPdf);

fprintf('\nDone.\n');

%% ============================================================
% LOCAL FUNCTION
% ============================================================

function zOut = getSurfaceZ(~,~,~)

    % Placeholder used only to initialize the contour array.
    zOut = [];

end