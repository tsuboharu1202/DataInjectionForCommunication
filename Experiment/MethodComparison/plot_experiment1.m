function plot_experiment1(results)
% plot_experiment1: 実験1（保守性）の結果をプロット
%
% 入力:
%   results: results_exp1 構造体
%       - param_list: パラメータリスト
%       - delta_ori_mean: original_thesisのdelta平均
%       - delta_reg_mean: regularization_sdpのdelta平均
%       - delta_ori_std: original_thesisのdelta標準偏差
%       - delta_reg_std: regularization_sdpのdelta標準偏差
%       - rho_rough: 理論的下界

param_list = results.param_list;
delta_ori_mean = results.delta_ori_mean;
delta_reg_mean = results.delta_reg_mean;
delta_ori_std = results.delta_ori_std;
delta_reg_std = results.delta_reg_std;
rho_rough = results.rho_rough;

% グラフ作成
figure('Name', '実験1: 保守性の評価', 'Position', [100, 100, 800, 500]);

% X軸位置（カテゴリカル）
x = 1:length(param_list);
width = 0.35;

% バーグラフ
hold on;

% original_thesis
bar_ori = bar(x - width/2, delta_ori_mean, width, 'FaceColor', [0.2, 0.4, 0.8]);
errorbar(x - width/2, delta_ori_mean, delta_ori_std, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% regularization_sdp
bar_reg = bar(x + width/2, delta_reg_mean, width, 'FaceColor', [0.8, 0.3, 0.3]);
errorbar(x + width/2, delta_reg_mean, delta_reg_std, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% rho_roughの参照線
yline(rho_rough, '--', 'Color', [0.3, 0.7, 0.3], 'LineWidth', 2, ...
    'Label', sprintf('rho\\_rough = %.3e', rho_rough), 'LabelHorizontalAlignment', 'left');

hold off;

% 装飾
xlabel('パラメータ (Phi11係数 / epsilon)', 'FontSize', 12);
ylabel('delta', 'FontSize', 12);
title('実験1: 保守性の評価（攻撃なし）', 'FontSize', 14);
legend([bar_ori, bar_reg], {'original\_thesis', 'regularization\_sdp'}, ...
    'Location', 'best', 'FontSize', 11);

% X軸ラベル
xticks(x);
xticklabels(arrayfun(@(p) sprintf('%.0e', p), param_list, 'UniformOutput', false));

% Y軸を対数スケールに
set(gca, 'YScale', 'log');

grid on;
set(gca, 'FontSize', 11);

% 保存
saveas(gcf, 'experiment1_conservatism.png');
saveas(gcf, 'experiment1_conservatism.fig');
fprintf('実験1のグラフを保存しました: experiment1_conservatism.png\n');

end

