%% ============================================================
% DIRECT PRESSURE DRAG FROM Cp.csv
%
% Cp.csv:
%   x/c, z/c, Cp
%
% Computes:
%   Cx,p
%   Cz,p
%   Cd,p
%   Cl,p
%
% Uses the 2D airfoil contour reconstructed from
% upper and lower surfaces.
% ============================================================

clear;
close all;
clc;

%% ============================================================
% SETTINGS
% ============================================================

AoA = 10;

cpFile = ...
    'naca0012LES_aoa10/comparison_results/Cp.csv';

%% ============================================================
% READ Cp
% ============================================================

data = readmatrix(cpFile);

data = data(all(isfinite(data),2),:);

if size(data,2) < 3
    error('Cp.csv must contain x/c, z/c, Cp.');
end

x = data(:,1);
z = data(:,2);
Cp = data(:,3);

%% ============================================================
% SPLIT UPPER / LOWER
% ============================================================

upper = z >= 0;
lower = z < 0;

xu = x(upper);
zu = z(upper);
cpu = Cp(upper);

xl = x(lower);
zl = z(lower);
cpl = Cp(lower);

%% ============================================================
% SORT BY x/c
% ============================================================

[xu,iu] = sort(xu);
zu = zu(iu);
cpu = cpu(iu);

[xl,il] = sort(xl);
zl = zl(il);
cpl = cpl(il);

%% ============================================================
% REMOVE DUPLICATE x/c VALUES
% ============================================================

[xu,iu] = unique(xu,'stable');
zu = zu(iu);
cpu = cpu(iu);

[xl,il] = unique(xl,'stable');
zl = zl(il);
cpl = cpl(il);

%% ============================================================
% BUILD CLOSED AIRFOIL CONTOUR
%
% Upper:
%   LE -> TE
%
% Lower:
%   TE -> LE
%
% Therefore:
%   upper followed by reversed lower
% ============================================================

xc = [
    xu
    flipud(xl)
];

zc = [
    zu
    flipud(zl)
];

cpc = [
    cpu
    flipud(cpl)
];

%% ============================================================
% REMOVE DUPLICATE/VERY CLOSE POINTS
% ============================================================

ds = sqrt( ...
    diff(xc).^2 + ...
    diff(zc).^2);

keep = [true; ds > 1e-12];

xc = xc(keep);
zc = zc(keep);
cpc = cpc(keep);

%% ============================================================
% CLOSE CONTOUR
% ============================================================

if xc(end) ~= xc(1) || zc(end) ~= zc(1)

    xc(end+1) = xc(1);
    zc(end+1) = zc(1);
    cpc(end+1) = cpc(1);

end

%% ============================================================
% CHECK CONTOUR ORIENTATION
% ============================================================

areaSigned = 0.5 * sum( ...
    xc(1:end-1).*zc(2:end) - ...
    xc(2:end).*zc(1:end-1));

fprintf('\nSigned contour area = %.8e\n',areaSigned);

% Reverse contour if necessary
if areaSigned < 0

    xc = flipud(xc);
    zc = flipud(zc);
    cpc = flipud(cpc);

end

%% ============================================================
% INTEGRATE PRESSURE FORCE
%
% For unit chord:
%
%   Cx = - integral(Cp dz)
%   Cz =   integral(Cp dx)
%
% Midpoint integration is used.
% ============================================================

Cx = 0;
Cz = 0;

for i = 1:length(xc)-1

    dx = xc(i+1) - xc(i);
    dz = zc(i+1) - zc(i);

    CpMid = 0.5 * ...
        (cpc(i) + cpc(i+1));

    Cx = Cx - CpMid * dz;

    Cz = Cz + CpMid * dx;

end

%% ============================================================
% CONVERT TO DRAG / LIFT
% ============================================================

alpha = deg2rad(AoA);

CDp = ...
    Cx*cos(alpha) + ...
    Cz*sin(alpha);

CLp = ...
    -Cx*sin(alpha) + ...
    Cz*cos(alpha);

%% ============================================================
% PRINT RESULTS
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('DIRECT PRESSURE FORCE FROM Cp.csv\n');
fprintf('============================================================\n');

fprintf('\nAoA = %.2f deg\n',AoA);

fprintf('\nPressure force coefficients:\n');

fprintf('  Cx,p = %.8f\n',Cx);
fprintf('  Cz,p = %.8f\n',Cz);

fprintf('\nLift / drag:\n');

fprintf('  Cd,p = %.8f\n',CDp);
fprintf('  Cl,p = %.8f\n',CLp);

fprintf('\n');

%% ============================================================
% PLOT RECONSTRUCTED CONTOUR
% ============================================================

figure( ...
    'Color','w');

hold on;
box on;
grid on;

plot( ...
    xc, ...
    zc, ...
    'b-', ...
    'LineWidth',1.2);

axis equal;

xlabel('x/c');
ylabel('z/c');

title('Reconstructed NACA 0012 contour');

%% ============================================================
% CHECK AGAINST YOUR SUMMARY.CSV
% ============================================================

summaryFile = ...
    'naca0012LES_aoa10/comparison_results/summary.csv';

if isfile(summaryFile)

    summary = readtable(summaryFile);

    idx = abs(summary.AoA_deg - AoA) < 1e-8;

    if any(idx)

        row = summary(find(idx,1),:);

        fprintf('\nComparison with 3D force integration:\n');

        fprintf( ...
            '  3D Cd,p = %.8f\n', ...
            row.Cd_pressure);

        fprintf( ...
            '  2D Cd,p = %.8f\n', ...
            CDp);

        fprintf( ...
            '  Difference = %.8f\n', ...
            CDp - row.Cd_pressure);

        fprintf( ...
            '  Relative difference = %.3f %%\n', ...
            100 * ...
            (CDp-row.Cd_pressure) / ...
            abs(row.Cd_pressure));

        fprintf('\n');

        fprintf( ...
            '  3D Cl,p = %.8f\n', ...
            row.Cl_pressure);

        fprintf( ...
            '  2D Cl,p = %.8f\n', ...
            CLp);

    end

end