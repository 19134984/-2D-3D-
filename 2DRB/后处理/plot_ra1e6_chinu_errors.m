%% Ra = 1e6：四个流场分支的 chi_nu 相对误差趋势
% 数据来源：P100 当前 FLOW-TEST/Ra1e6 扫描。
% 所有算例均满足 run.status=0、完整历史 501 个样本、最终窗口 251 个样本，
% 且 Nu/Re 前后半窗相对误差均小于 1%。
% Zhang Table 1 基准：Nu = 6.30，Re = 279。
% 纵轴为有符号相对误差：100*(computed/reference - 1)，单位为 %。

clear; clc; close all;

chiNu = 0.0:0.1:0.9;
nuZhang = 6.30;
reZhang = 279.0;

% 行顺序：原始魔法参数、修正魔法参数、固定 Sq=1、流场 BGK。
nu = [ ...
    6.317437382980921, 6.318454273882961, 6.319471346437386, 6.320489966164613, 6.321512441301700, 6.322542986934131, 6.323590073796610, 6.324673573468429, 6.325853276955203, 6.327429414915536; ...
    6.317437384223000, 6.317441469475725, 6.317445246293720, 6.317448057639577, 6.317448474159084, 6.317443200113152, 6.317423870640674, 6.317365824067601, 6.317175434781114, 6.316266690746391; ...
    6.328011390294407, 6.328034088815191, 6.328060886356141, 6.328093134188153, 6.328132898800424, 6.328183536101220, 6.328250971132377, 6.328346948331016, 6.328499400575103, 6.328801045625344; ...
    6.331396229492414, 6.331406936853456, 6.331413838012160, 6.331415443835007, 6.331409361279170, 6.331391569769720, 6.331354663714464, 6.331282868146812, 6.331132755744535, 6.330697838469139];

re = [ ...
    279.2861509846898, 279.3417901103322, 279.3975298244149, 279.4533962311389, 279.5094340035673, 279.5657248699428, 279.6224339654435, 279.6799491947829, 279.7394403129910, 279.8067972831411; ...
    279.2861509942099, 279.2863320193791, 279.2865171736821, 279.2866922346116, 279.2868234424132, 279.2868285937921, 279.2864906117545, 279.2851509814524, 279.2802591998765, 279.2556489661146; ...
    279.8424567234592, 279.8428892495791, 279.8433964996934, 279.8440019145983, 279.8447407050305, 279.8456686905238, 279.8468809627661, 279.8485561792194, 279.8510796324987, 279.8554239212353; ...
    279.9386994089824, 279.9385867722600, 279.9383219318762, 279.9378523183869, 279.9370947594091, 279.9359100645270, 279.9340426076763, 279.9309464065867, 279.9251088248490, 279.9092980086878];

branchNames = { ...
    'Original magic TRT', ...
    'Corrected magic TRT', ...
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
xlim([0.0, 0.9]);
xticks(chiNu);
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
xlim([0.0, 0.9]);
xticks(chiNu);
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
