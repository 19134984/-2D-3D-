%% Ra = 1e6：四个流场分支的 chi_nu 误差趋势
% 数据来源：P100 当前 FLOW-TEST/Ra1e6 完整统计窗口（501 个样本）。
% Zhang Table 1 基准：Nu = 6.3，Re = 279。
% 纵轴为有符号相对误差 100*(computed/reference - 1)，不是绝对值。

clear; clc; close all;

chiNu = 0.0:0.1:0.9;
nuZhang = 6.3;
reZhang = 279.0;

% 行顺序：原始魔法参数、严格修正魔法参数、固定 Sq=1、速度场 BGK。
nu = [ ...
    6.317437517900302, 6.318454409973400, 6.319471483436159, 6.320490102420832, 6.321512579776750, 6.322543124613951, 6.323590213328860, 6.324673711626106, 6.325853416946512, 6.327429554828948; ...
    6.317437517900301, 6.317441603059550, 6.317445380593760, 6.317448192725939, 6.317448609353887, 6.317443334281808, 6.317424006761923, 6.317365960961512, 6.317175571555908, 6.316266829696051; ...
    6.328011531500073, 6.328034230793171, 6.328061029130828, 6.328093276967998, 6.328133040329153, 6.328183677626624, 6.328251113061131, 6.328347090475629, 6.328499542841072, 6.328801187962280; ...
    6.331396375231226, 6.331407081898051, 6.331413983435105, 6.331415588260574, 6.331409507391085, 6.331391715646549, 6.331354808456999, 6.331283013929860, 6.331132900923720, 6.330697983271782];

re = [ ...
    279.2861509890117, 279.3417901597084, 279.3975299102407, 279.4533963605305, 279.5094341748340, 279.5657250774222, 279.6224342071531, 279.6799494589842, 279.7394406068576, 279.8067975888020; ...
    279.2861509890117, 279.2863320274637, 279.2865171836087, 279.2866922319185, 279.2868234501162, 279.2868285836082, 279.2864906147070, 279.2851509777359, 279.2802591965587, 279.2556489548828; ...
    279.8424570592343, 279.8428895893053, 279.8433968360359, 279.8440022495212, 279.8447410359571, 279.8456690178600, 279.8468812886531, 279.8485564998164, 279.8510799433454, 279.8554242105051; ...
    279.9386996507635, 279.9385870192313, 279.9383221764408, 279.9378525540281, 279.9370949936199, 279.9359102992140, 279.9340428347642, 279.9309466348902, 279.9251090493903, 279.9092982357992];

branchNames = { ...
    'Original magic TRT', ...
    'Effective magic TRT', ...
    'Fixed S_q = 1', ...
    'Flow BGK (S_q = S_\nu)'};

colors = [ ...
    0.8500, 0.3250, 0.0980; ...
    0.0000, 0.4470, 0.7410; ...
    0.4660, 0.6740, 0.1880; ...
    0.4940, 0.1840, 0.5560];
markers = {'o', 's', '^', 'd'};
lineStyles = {'-', '-', '--', '-.'};

nuErrorPct = 100.0 * (nu ./ nuZhang - 1.0);
reErrorPct = 100.0 * (re ./ reZhang - 1.0);

fig = figure('Color', 'w', 'Position', [120, 100, 860, 700]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; box on; grid on;
for k = 1:numel(branchNames)
    plot(chiNu, nuErrorPct(k, :), ...
        'Color', colors(k, :), 'LineStyle', lineStyles{k}, ...
        'Marker', markers{k}, 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'DisplayName', branchNames{k});
end
yline(0.0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlim([0.0, 0.9]); xticks(chiNu);
ylabel('Relative error of Nu (%)', 'Interpreter', 'none');
title('Nu relative to Zhang', 'FontWeight', 'normal');
legend('Location', 'best', 'Interpreter', 'tex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

nexttile;
hold on; box on; grid on;
for k = 1:numel(branchNames)
    plot(chiNu, reErrorPct(k, :), ...
        'Color', colors(k, :), 'LineStyle', lineStyles{k}, ...
        'Marker', markers{k}, 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'DisplayName', branchNames{k});
end
yline(0.0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlim([0.0, 0.9]); xticks(chiNu);
xlabel('\chi_\nu', 'Interpreter', 'tex');
ylabel('Relative error of Re (%)', 'Interpreter', 'none');
title('Re relative to Zhang', 'FontWeight', 'normal');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

title(layout, 'Ra = 10^6', ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal');

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    outputDir = pwd;
else
    outputDir = fileparts(scriptPath);
end

exportgraphics(fig, fullfile(outputDir, 'ra1e6_chinu_error_trends.png'), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(outputDir, 'ra1e6_chinu_error_trends.pdf'), ...
    'ContentType', 'vector');
