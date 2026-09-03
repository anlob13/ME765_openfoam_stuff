%% ============================================================
% NACA 0012 LES: AoA = 10 deg
%
% Compare:
%   1. LES Cp against Ladson et al.
%   2. LES CL against Gregory
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
% summary.csv:
%   time, AoA_deg, Cl, Cd, Cl_pressure, ...
%
% ============================================================

clear;
close all;
clc;

%% ============================================================
% SETTINGS
% ============================================================

AoA = 10;

caseDir = 'fine_naca0012LES_aoa10/comparison_results';

cpFile = fullfile(caseDir,'Cp.csv');
summaryFile = fullfile(caseDir,'summary.csv');

expFile = 'exp.dat';
gregoryFile = 'Gregory_CL_expdata.dat';

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
        error('Could not find file:\n%s',requiredFiles{i});
    end
end

%% ============================================================
% READ LES Cp
% ============================================================

cpData = readmatrix(cpFile);

% Remove header/non-numeric rows
cpData = cpData(all(isfinite(cpData),2),:);

if size(cpData,2) < 3
    error('%s must contain x/c, z/c, Cp.',cpFile);
end

LES_x  = cpData(:,1);
LES_z  = cpData(:,2);
LES_Cp = cpData(:,3);

% ------------------------------------------------------------
% Split upper/lower
% ------------------------------------------------------------

upperMask = LES_z >= 0;
lowerMask = LES_z < 0;

LES_xu  = LES_x(upperMask);
LES_cpu = LES_Cp(upperMask);

LES_xl  = LES_x(lowerMask);
LES_cpl = LES_Cp(lowerMask);

% Sort by x/c
[LES_xu,idx] = sort(LES_xu);
LES_cpu = LES_cpu(idx);

[LES_xl,idx] = sort(LES_xl);
LES_cpl = LES_cpl(idx);

fprintf('\nLES Cp data:\n');
fprintf('  Total points = %d\n',length(LES_x));
fprintf('  Upper points = %d\n',length(LES_xu));
fprintf('  Lower points = %d\n',length(LES_xl));

%% ============================================================
% READ LES SUMMARY
% ============================================================

summary = readtable(summaryFile);

% Find AoA = 10 row
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
    error('Could not open experimental file: %s',expFile);
end

expData = struct( ...
    'alpha',{}, ...
    'x',{}, ...
    'cp',{});

currentAlpha = NaN;
xData = [];
cpData = [];

while true

    line = fgetl(fid);

    % End of file
    if isequal(line,-1)
        break;
    end

    % fgetl should return char, but force row form
    line = char(line);

    % Remove whitespace
    line = strtrim(line);

    if isempty(line)
        continue;
    end

    % --------------------------------------------------------
    % Zone line
    %
    % Example:
    % zone, t="Re=3 million, alpha=10.0130, fixed transition"
    % --------------------------------------------------------

    if length(line) >= 4 && strncmpi(line,'zone',4)

        % Save previous zone
        if ~isempty(xData) && ~isnan(currentAlpha)

            expData(end+1).alpha = currentAlpha;
            expData(end).x = xData;
            expData(end).cp = cpData;

        end

        % Extract alpha
        token = regexp( ...
            line, ...
            'alpha\s*=\s*([-+]?\d*\.?\d+)', ...
            'tokens', ...
            'once');

        currentAlpha = NaN;

        if ~isempty(token)
            currentAlpha = str2double(token{1});
        end

        % Reset current zone
        xData = [];
        cpData = [];

        continue;
    end

    % Ignore comments
    if line(1) == '#'
        continue;
    end

    % Ignore variables line
    if length(line) >= 9 && strncmpi(line,'variables',9)
        continue;
    end

    % --------------------------------------------------------
    % Numeric x/c and Cp
    % --------------------------------------------------------

    vals = sscanf(line,'%f %f');

    if numel(vals) >= 2

        xData(end+1,1) = vals(1);
        cpData(end+1,1) = vals(2);

    end
end

fclose(fid);

% Save final zone
if ~isempty(xData) && ~isnan(currentAlpha)

    expData(end+1).alpha = currentAlpha;
    expData(end).x = xData;
    expData(end).cp = cpData;

end

%% ============================================================
% CHECK LADSON DATA
% ============================================================

if isempty(expData)
    error('No Ladson experimental zones were found.');
end

fprintf('\nLadson data:\n');

for i = 1:length(expData)

    fprintf( ...
        '  alpha = %.4f deg, %d points\n', ...
        expData(i).alpha, ...
        length(expData(i).x));

end

%% ============================================================
% SELECT LADSON DATA NEAREST AoA = 10
% ============================================================

[minDiff,idx] = min(abs([expData.alpha] - AoA));

if minDiff > 0.5

    error( ...
        'No Ladson Cp data within 0.5 deg of AoA = %.1f deg.', ...
        AoA);

end

Ladson = expData(idx);

fprintf('\nSelected Ladson data:\n');
fprintf('  Requested AoA = %.1f deg\n',AoA);
fprintf('  Actual AoA    = %.4f deg\n',Ladson.alpha);
fprintf('  Points        = %d\n',length(Ladson.x));

%% ============================================================
% SPLIT LADSON INTO UPPER / LOWER
%
% Your Ladson file is ordered:
%
%   upper surface -> leading edge -> lower surface
%
% with x/c = 0 at the leading edge.
% ============================================================

[~,iLE] = min(abs(Ladson.x));

% First branch = upper
Ladson_xu = Ladson.x(1:iLE);
Ladson_cpu = Ladson.cp(1:iLE);

% Second branch = lower
Ladson_xl = Ladson.x(iLE:end);
Ladson_cpl = Ladson.cp(iLE:end);

% Sort by increasing x/c
[Ladson_xu,idx] = sort(Ladson_xu);
Ladson_cpu = Ladson_cpu(idx);

[Ladson_xl,idx] = sort(Ladson_xl);
Ladson_cpl = Ladson_cpl(idx);

% Remove duplicates
[Ladson_xu,idx] = unique(Ladson_xu,'stable');
Ladson_cpu = Ladson_cpu(idx);

[Ladson_xl,idx] = unique(Ladson_xl,'stable');
Ladson_cpl = Ladson_cpl(idx);

fprintf('\nLadson surfaces:\n');
fprintf('  Upper points = %d\n',length(Ladson_xu));
fprintf('  Lower points = %d\n',length(Ladson_xl));

%% ============================================================
% INTERPOLATE LES Cp ONTO LADSON UPPER SURFACE
% ============================================================

LES_Cp_u_interp = interp1( ...
    LES_xu, ...
    LES_cpu, ...
    Ladson_xu, ...
    'linear', ...
    NaN);

validU = ...
    isfinite(LES_Cp_u_interp) & ...
    isfinite(Ladson_cpu);

Cp_error_upper = ...
    LES_Cp_u_interp(validU) - ...
    Ladson_cpu(validU);

%% ============================================================
% INTERPOLATE LES Cp ONTO LADSON LOWER SURFACE
% ============================================================

LES_Cp_l_interp = interp1( ...
    LES_xl, ...
    LES_cpl, ...
    Ladson_xl, ...
    'linear', ...
    NaN);

validL = ...
    isfinite(LES_Cp_l_interp) & ...
    isfinite(Ladson_cpl);

Cp_error_lower = ...
    LES_Cp_l_interp(validL) - ...
    Ladson_cpl(validL);

%% ============================================================
% Cp ERRORS
% ============================================================

Cp_upper_RMSE = sqrt( ...
    mean(Cp_error_upper.^2));

Cp_lower_RMSE = sqrt( ...
    mean(Cp_error_lower.^2));

Cp_upper_MAE = mean( ...
    abs(Cp_error_upper));

Cp_lower_MAE = mean( ...
    abs(Cp_error_lower));

Cp_all_error = [
    Cp_error_upper
    Cp_error_lower
];

Cp_RMSE = sqrt( ...
    mean(Cp_all_error.^2));

Cp_MAE = mean( ...
    abs(Cp_all_error));

%% ============================================================
% READ GREGORY CL DATA
% ============================================================

fid = fopen(gregoryFile,'rt');

if fid == -1
    error('Could not open Gregory file: %s',gregoryFile);
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

    % Comments/header
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

    % Numeric alpha, CL
    vals = sscanf(line,'%f %f');

    if numel(vals) >= 2

        gregoryAlpha(end+1,1) = vals(1);
        gregoryCL(end+1,1) = vals(2);

    end
end

fclose(fid);

valid = ...
    isfinite(gregoryAlpha) & ...
    isfinite(gregoryCL);

gregoryAlpha = gregoryAlpha(valid);
gregoryCL = gregoryCL(valid);

%% ============================================================
% GREGORY CL AT AoA = 10
% ============================================================

Gregory_CL = interp1( ...
    gregoryAlpha, ...
    gregoryCL, ...
    AoA, ...
    'linear', ...
    'extrap');

%% ============================================================
% CL ERROR
% ============================================================

CL_error = LES_CL - Gregory_CL;

CL_abs_error = abs(CL_error);

if abs(Gregory_CL) > 1e-12

    CL_percent_error = ...
        100 * CL_error / abs(Gregory_CL);

else

    CL_percent_error = NaN;

end

%% ============================================================
% RESULTS TABLE
% ============================================================

results = table( ...
    AoA, ...
    Ladson.alpha, ...
    Gregory_CL, ...
    LES_CL, ...
    CL_error, ...
    CL_abs_error, ...
    CL_percent_error, ...
    LES_CLp, ...
    LES_CLf, ...
    LES_CD, ...
    LES_CDp, ...
    LES_CDf, ...
    Cp_upper_RMSE, ...
    Cp_lower_RMSE, ...
    Cp_RMSE, ...
    Cp_upper_MAE, ...
    Cp_lower_MAE, ...
    Cp_MAE, ...
    'VariableNames',{ ...
    'AoA_deg', ...
    'Ladson_AoA_deg', ...
    'Gregory_CL', ...
    'LES_CL', ...
    'CL_error', ...
    'CL_abs_error', ...
    'CL_percent_error', ...
    'LES_CL_pressure', ...
    'LES_CL_viscous', ...
    'LES_CD', ...
    'LES_CD_pressure', ...
    'LES_CD_viscous', ...
    'Cp_upper_RMSE', ...
    'Cp_lower_RMSE', ...
    'Cp_RMSE', ...
    'Cp_upper_MAE', ...
    'Cp_lower_MAE', ...
    'Cp_MAE'});

%% ============================================================
% DISPLAY RESULTS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('NACA 0012 LES COMPARISON - AoA = %.1f deg\n',AoA);
fprintf('============================================================\n');

fprintf('\nCL comparison:\n');
fprintf('  Gregory CL       = %.6f\n',Gregory_CL);
fprintf('  LES CL            = %.6f\n',LES_CL);
fprintf('  Signed error      = %.6f\n',CL_error);
fprintf('  Absolute error    = %.6f\n',CL_abs_error);
fprintf('  Percent error     = %.3f %%\n',CL_percent_error);

fprintf('\nLES force breakdown:\n');
fprintf('  CL pressure       = %.6f\n',LES_CLp);
fprintf('  CL viscous        = %.6f\n',LES_CLf);
fprintf('  CL total          = %.6f\n',LES_CL);

fprintf('  CD pressure       = %.6f\n',LES_CDp);
fprintf('  CD viscous        = %.6f\n',LES_CDf);
fprintf('  CD total          = %.6f\n',LES_CD);

fprintf('\nCp comparison against Ladson:\n');
fprintf('  Upper RMSE        = %.6f\n',Cp_upper_RMSE);
fprintf('  Lower RMSE        = %.6f\n',Cp_lower_RMSE);
fprintf('  Total RMSE        = %.6f\n',Cp_RMSE);
fprintf('  Upper MAE         = %.6f\n',Cp_upper_MAE);
fprintf('  Lower MAE         = %.6f\n',Cp_lower_MAE);
fprintf('  Total MAE         = %.6f\n',Cp_MAE);

fprintf('\n');

disp(results);

%% ============================================================
% SAVE RESULTS
% ============================================================

resultsFile = fullfile( ...
    caseDir, ...
    'AoA10_Ladson_Gregory_comparison.csv');

writetable(results,resultsFile);

fprintf('Saved:\n  %s\n',resultsFile);

%% ============================================================
% SAVE Cp ERROR DATA
% ============================================================

CpComparison = table( ...
    [Ladson_xu(validU);Ladson_xl(validL)], ...
    [ones(sum(validU),1);-ones(sum(validL),1)], ...
    [Ladson_cpu(validU);Ladson_cpl(validL)], ...
    [LES_Cp_u_interp(validU);LES_Cp_l_interp(validL)], ...
    Cp_all_error, ...
    'VariableNames',{ ...
    'x_c', ...
    'surface', ...
    'Cp_Ladson', ...
    'Cp_LES', ...
    'Cp_error'});

cpComparisonFile = fullfile( ...
    caseDir, ...
    'AoA10_Cp_error.csv');

writetable(CpComparison,cpComparisonFile);

fprintf('Saved:\n  %s\n',cpComparisonFile);

%% ============================================================
% Cp PLOT
% ============================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[1 1 8 6]);

ax = axes('Parent',fig);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

% Ladson upper
plot( ...
    ax, ...
    Ladson_xu, ...
    Ladson_cpu, ...
    'ko', ...
    'MarkerSize',3, ...
    'LineStyle','none', ...
    'DisplayName','Ladson et al.');

% Ladson lower
plot( ...
    ax, ...
    Ladson_xl, ...
    Ladson_cpl, ...
    'ko', ...
    'MarkerSize',3, ...
    'LineStyle','none', ...
    'HandleVisibility','off');

% LES upper
plot( ...
    ax, ...
    LES_xu, ...
    LES_cpu, ...
    'b-', ...
    'LineWidth',1.2, ...
    'DisplayName','LES');

% LES lower
plot( ...
    ax, ...
    LES_xl, ...
    LES_cpl, ...
    'b-', ...
    'LineWidth',1.2, ...
    'HandleVisibility','off');

set(ax,'YDir','reverse');

xlim(ax,[0 1]);

xlabel(ax,'$x/c$', ...
    'Interpreter','latex');

ylabel(ax,'$C_p$', ...
    'Interpreter','latex');

title(ax, ...
    sprintf('NACA 0012, $\\alpha=%.0f^\\circ$',AoA), ...
    'Interpreter','latex');

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.7, ...
    'TickDir','out', ...
    'Layer','top');

%% ============================================================
% SAVE Cp PDF
% ============================================================

cpPdf = fullfile( ...
    caseDir, ...
    'Cp_AoA10_Ladson_comparison.pdf');

exportgraphics( ...
    ax, ...
    cpPdf, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

close(fig);

fprintf('Saved:\n  %s\n',cpPdf);

%% ============================================================
% CL PLOT
% ============================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[1 1 8 6]);

ax = axes('Parent',fig);

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

% Gregory curve
plot( ...
    ax, ...
    gregoryAlpha, ...
    gregoryCL, ...
    'ko-', ...
    'LineWidth',0.8, ...
    'MarkerSize',3, ...
    'DisplayName','Gregory');

% LES result
plot( ...
    ax, ...
    AoA, ...
    LES_CL, ...
    'bo', ...
    'LineWidth',1.2, ...
    'MarkerSize',5, ...
    'DisplayName','LES');

% Highlight comparison AoA
xline( ...
    ax, ...
    AoA, ...
    ':k', ...
    'HandleVisibility','off');

xlabel(ax,'$\alpha$ [deg]', ...
    'Interpreter','latex');

ylabel(ax,'$C_L$', ...
    'Interpreter','latex');

title(ax, ...
    'NACA 0012 Lift Comparison', ...
    'Interpreter','latex');

legend( ...
    ax, ...
    'Location','best', ...
    'Box','off');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'LineWidth',0.7, ...
    'TickDir','out', ...
    'Layer','top');

%% ============================================================
% SAVE CL PDF
% ============================================================

clPdf = fullfile( ...
    caseDir, ...
    'CL_AoA10_Gregory_comparison.pdf');

exportgraphics( ...
    ax, ...
    clPdf, ...
    'ContentType','vector', ...
    'BackgroundColor','white');

close(fig);

fprintf('Saved:\n  %s\n',clPdf);

%% ============================================================
% FINAL SUMMARY
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nAoA = %.1f deg\n',AoA);

fprintf('\nCL:\n');
fprintf('  Gregory = %.6f\n',Gregory_CL);
fprintf('  LES     = %.6f\n',LES_CL);
fprintf('  Error   = %.6f\n',CL_error);
fprintf('  %% Error = %.3f %%\n',CL_percent_error);

fprintf('\nCp RMSE:\n');
fprintf('  Upper = %.6f\n',Cp_upper_RMSE);
fprintf('  Lower = %.6f\n',Cp_lower_RMSE);
fprintf('  Total = %.6f\n',Cp_RMSE);

fprintf('\nFiles written:\n');
fprintf('  %s\n',resultsFile);
fprintf('  %s\n',cpComparisonFile);
fprintf('  %s\n',cpPdf);
fprintf('  %s\n',clPdf);

fprintf('\nDone.\n');

q = 0.5*1.225*43.824^2;
Aref = 0.2;

Cd_check = row.Fx/(q*Aref);
Cl_check = row.Fz/(q*Aref);

fprintf('Cd from Fx = %.6f\n',Cd_check);
fprintf('Cl from Fz = %.6f\n',Cl_check);
fprintf('Reported Cd = %.6f\n',row.Cd);
fprintf('Reported Cl = %.6f\n',row.Cl);