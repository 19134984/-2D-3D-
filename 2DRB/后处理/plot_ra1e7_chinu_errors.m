%% Ra = 1e7：原始与修正魔法参数分支的 chi_nu 绝对相对误差
% 网格：257 x 257；最终统计窗口：1050--1400 t_ff。
% Zhang Table 1 基准：Nu = 11.37，Re = 968。
% 纵轴采用绝对相对误差：100*abs(computed/reference - 1)，单位为 %。
% 原始分支 chi_nu = 0.5、0.7 的 Re 属于另一条统计解分支，不参与误差趋势比较。

clear; clc; close all;

chiNu = 0.0:0.1:0.9;
nuZhang = 11.37;
reZhang = 968.0;

% 行顺序：原始魔法参数、修正魔法参数。
nu = [ ...
    11.381782658, 11.361142109, 11.374216006, 11.365455288, 11.354660372, ...
    11.362469905, 11.351373550, 11.355030087, 11.353900815, 11.352738625; ...
    11.373182064, 11.368010248, 11.371935218, 11.366264777, 11.366025151, ...
    11.375363448, 11.385035567, 11.373957564, 11.384681271, 11.377382767];

re = [ ...
    965.714160, 964.753000, 965.227655, 964.778410, 964.801215, ...
    1045.729901, 965.741344, 1044.497893, 966.757592, 966.029251; ...
    964.616701, 964.333564, 964.424391, 964.306433, 964.394978, ...
    964.606121, 966.959972, 964.483347, 967.077478, 964.966651];

branchNames = { ...
    'Original magic TRT', ...
    'Corrected magic TRT'};

colors = [ ...
    0.8500, 0.3250, 0.0980; ...
    0.0000, 0.4470, 0.7410];
markers = {'o', 's'};

nuAbsErrorPct = 100.0 * abs(nu ./ nuZhang - 1.0);
reAbsErrorPct = 100.0 * abs(re ./ reZhang - 1.0);

% 仅从原始 Re 曲线中排除另一条统计解分支的两个点；修正分支保持完整。
rePointIncluded = true(size(reAbsErrorPct));
rePointIncluded(1, chiNu == 0.5 | chiNu == 0.7) = false;

fig = figure('Color', 'w', 'Position', [120, 100, 860, 700]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; box on; grid on;
for k = 1:numel(branchNames)
    plot(chiNu, nuAbsErrorPct(k, :), ...
        'Color', colors(k, :), 'LineStyle', '-', ...
        'Marker', markers{k}, 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'DisplayName', branchNames{k});
end
xlim([0.0, 0.9]);
xticks(chiNu);
ylabel('Absolute relative error of Nu (%)', 'Interpreter', 'none');
title('Nu relative to Zhang', 'FontWeight', 'normal');
legend('Location', 'best', 'Interpreter', 'tex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

nexttile;
hold on; box on; grid on;
for k = 1:numel(branchNames)
    included = rePointIncluded(k, :);
    plot(chiNu(included), reAbsErrorPct(k, included), ...
        'Color', colors(k, :), 'LineStyle', '-', ...
        'Marker', markers{k}, 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'DisplayName', branchNames{k});
end
xlim([0.0, 0.9]);
xticks(chiNu);
xlabel('\chi_\nu', 'Interpreter', 'tex');
ylabel('Absolute relative error of Re (%)', 'Interpreter', 'none');
title('Re relative to Zhang (original \chi_\nu = 0.5, 0.7 excluded)', ...
    'FontWeight', 'normal', 'Interpreter', 'tex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

title(layout, 'Ra = 10^7', ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal');

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    outputDir = pwd;
else
    outputDir = fileparts(scriptPath);
end

exportgraphics(fig, fullfile(outputDir, 'ra1e7_chinu_absolute_error_trends.png'), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(outputDir, 'ra1e7_chinu_absolute_error_trends.pdf'), ...
    'ContentType', 'vector');
