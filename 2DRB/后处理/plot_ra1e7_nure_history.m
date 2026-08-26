%% Ra = 1e7 原始魔法参数分支：Nu(t) 与 Re(t) 历史
% 数据来自服务器各工况自己的 results/ 或 merged_results/。
% 旧版本只从统计窗口起点开始保存历史，因此：
%   chi_nu = 0.4, 0.7, 0.8：现有数据为 700--1400 t_ff；
%   chi_nu = 0.5, 0.6：合并后的连续数据为 700--2800 t_ff。
% 图中不跨越缺失区间连线，也不把缺失数据插值成 0--700 t_ff。

clear; clc; close all;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    outputDir = pwd;
else
    outputDir = fileparts(scriptPath);
end
dataDir = fullfile(outputDir, 'ra1e7_history_data');

chiNu = [0.4, 0.5, 0.6, 0.7, 0.8];
fileTags = {'chi0p4', 'chi0p5', 'chi0p6', 'chi0p7', 'chi0p8'};
colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880; ...
    0.4940, 0.1840, 0.5560; ...
    0.3010, 0.7450, 0.9330];

nuHistory = cell(size(chiNu));
reHistory = cell(size(chiNu));
for k = 1:numel(chiNu)
    nuFile = fullfile(dataDir, [fileTags{k}, '_Nu.dat']);
    reFile = fullfile(dataDir, [fileTags{k}, '_Re.dat']);
    assert(isfile(nuFile), 'Missing Nu history: %s', nuFile);
    assert(isfile(reFile), 'Missing Re history: %s', reFile);

    nuHistory{k} = readmatrix(nuFile, 'FileType', 'text');
    reHistory{k} = readmatrix(reFile, 'FileType', 'text');
    assert(size(nuHistory{k}, 2) >= 2 && size(reHistory{k}, 2) >= 2, ...
        'History file must contain time and value columns for chi_nu = %.1f.', chiNu(k));
    assert(isequal(nuHistory{k}(:, 1), reHistory{k}(:, 1)), ...
        'Nu/Re time axes do not match for chi_nu = %.1f.', chiNu(k));
    assert(all(diff(nuHistory{k}(:, 1)) > 0), ...
        'History time is not strictly increasing for chi_nu = %.1f.', chiNu(k));
end

fig = figure('Color', 'w', 'Position', [100, 80, 980, 760]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axNu = nexttile;
hold(axNu, 'on'); box(axNu, 'on'); grid(axNu, 'on');
for k = 1:numel(chiNu)
    t = nuHistory{k}(:, 1);
    y = nuHistory{k}(:, 2);
    lightColor = 0.68 .* colors(k, :) + 0.32;
    plot(axNu, t, y, '-', 'Color', lightColor, 'LineWidth', 0.55, ...
        'HandleVisibility', 'off');
    plot(axNu, t, movmean(y, 51, 'Endpoints', 'shrink'), '-', ...
        'Color', colors(k, :), 'LineWidth', 1.65, ...
        'DisplayName', sprintf('\\chi_\\nu = %.1f', chiNu(k)));
end
yline(axNu, 11.37, '--', 'Zhang: Nu = 11.37', ...
    'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.0, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xline(axNu, 700, ':', 'history starts', 'HandleVisibility', 'off');
xline(axNu, 1400, ':', 'short runs end', 'HandleVisibility', 'off');
xlim(axNu, [0, 2800]);
ylabel(axNu, 'Nu(t)');
title(axNu, 'Nusselt-number history (thin: instantaneous; thick: 51-t_{ff} moving mean)', ...
    'FontWeight', 'normal', 'Interpreter', 'tex');
legend(axNu, 'Location', 'southoutside', 'NumColumns', 5, 'Interpreter', 'tex');
set(axNu, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

axRe = nexttile;
hold(axRe, 'on'); box(axRe, 'on'); grid(axRe, 'on');
for k = 1:numel(chiNu)
    t = reHistory{k}(:, 1);
    y = reHistory{k}(:, 2);
    lightColor = 0.68 .* colors(k, :) + 0.32;
    plot(axRe, t, y, '-', 'Color', lightColor, 'LineWidth', 0.55, ...
        'HandleVisibility', 'off');
    plot(axRe, t, movmean(y, 51, 'Endpoints', 'shrink'), '-', ...
        'Color', colors(k, :), 'LineWidth', 1.65, ...
        'DisplayName', sprintf('\\chi_\\nu = %.1f', chiNu(k)));
end
yline(axRe, 968.0, '--', 'Zhang: Re = 968', ...
    'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.0, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xline(axRe, 700, ':', 'history starts', 'HandleVisibility', 'off');
xline(axRe, 1400, ':', 'short runs end', 'HandleVisibility', 'off');
xlim(axRe, [0, 2800]);
xlabel(axRe, 't/t_{ff}', 'Interpreter', 'tex');
ylabel(axRe, 'Re(t)');
title(axRe, 'Reynolds-number history', 'FontWeight', 'normal');
legend(axRe, 'Location', 'southoutside', 'NumColumns', 5, 'Interpreter', 'tex');
set(axRe, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

title(layout, 'Ra = 10^7, original magic parameter branch', ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal');

exportgraphics(fig, fullfile(outputDir, 'ra1e7_original_nure_history.png'), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(outputDir, 'ra1e7_original_nure_history.pdf'), ...
    'ContentType', 'vector');

fprintf('Saved %s\n', fullfile(outputDir, 'ra1e7_original_nure_history.png'));
fprintf('Saved %s\n', fullfile(outputDir, 'ra1e7_original_nure_history.pdf'));
for k = 1:numel(chiNu)
    fprintf('chi_nu=%.1f: %d rows, %.0f--%.0f t_ff\n', chiNu(k), ...
        size(reHistory{k}, 1), reHistory{k}(1, 1), reHistory{k}(end, 1));
end
