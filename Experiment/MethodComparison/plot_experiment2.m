function plot_experiment2(results)
% plot_experiment2: 実験2（ロバスト性）の結果をプロット
%
% 4つのグラフを生成:
% 1. original_thesis + 判定1 (delta > rho_rough)
% 2. original_thesis + 判定2 (hinf > 1/delta)
% 3. regularization_sdp + 判定1 (delta > rho_rough)
% 4. regularization_sdp + 判定2 (hinf > 1/delta)
%
% 各グラフ:
%   横軸: 攻撃の大きさ
%   縦軸: 攻撃成功率 (%)
%   3本のライン: パラメータ別

param_list = results.param_list;
attack_sizes = results.attack_sizes;

% カラーマップ
colors = [
    0.2, 0.4, 0.8;   % 青
    0.8, 0.3, 0.3;   % 赤
    0.3, 0.7, 0.3;   % 緑
];

markers = {'o', 's', 'd'};
line_styles = {'-', '--', ':'};

% Figure作成（2x2のサブプロット）
figure('Name', '実験2: ロバスト性の評価', 'Position', [100, 100, 1200, 900]);

% === グラフ1: original_thesis + 判定1 ===
subplot(2, 2, 1);
hold on;
for p_idx = 1:length(param_list)
    success_rate = results.ori_success_rate_delta(p_idx, :);
    plot(1:length(attack_sizes), success_rate, ...
        'Color', colors(p_idx, :), ...
        'LineWidth', 2, ...
        'Marker', markers{p_idx}, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', colors(p_idx, :), ...
        'DisplayName', sprintf('param = %.0e', param_list(p_idx)));
end
hold off;
xlabel('攻撃の大きさ', 'FontSize', 11);
ylabel('攻撃成功率 (%)', 'FontSize', 11);
title('original\_thesis: delta > rho\_rough × 1.01', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
xticks(1:length(attack_sizes));
xticklabels(arrayfun(@(a) sprintf('%.0e', a), attack_sizes, 'UniformOutput', false));
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 10);

% === グラフ2: original_thesis + 判定2 ===
subplot(2, 2, 2);
hold on;
for p_idx = 1:length(param_list)
    success_rate = results.ori_success_rate_hinf(p_idx, :);
    plot(1:length(attack_sizes), success_rate, ...
        'Color', colors(p_idx, :), ...
        'LineWidth', 2, ...
        'Marker', markers{p_idx}, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', colors(p_idx, :), ...
        'DisplayName', sprintf('param = %.0e', param_list(p_idx)));
end
hold off;
xlabel('攻撃の大きさ', 'FontSize', 11);
ylabel('攻撃成功率 (%)', 'FontSize', 11);
title('original\_thesis: hinf > (1/delta) × 1.01', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
xticks(1:length(attack_sizes));
xticklabels(arrayfun(@(a) sprintf('%.0e', a), attack_sizes, 'UniformOutput', false));
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 10);

% === グラフ3: regularization_sdp + 判定1 ===
subplot(2, 2, 3);
hold on;
for p_idx = 1:length(param_list)
    success_rate = results.reg_success_rate_delta(p_idx, :);
    plot(1:length(attack_sizes), success_rate, ...
        'Color', colors(p_idx, :), ...
        'LineWidth', 2, ...
        'Marker', markers{p_idx}, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', colors(p_idx, :), ...
        'DisplayName', sprintf('epsilon = %.0e', param_list(p_idx)));
end
hold off;
xlabel('攻撃の大きさ', 'FontSize', 11);
ylabel('攻撃成功率 (%)', 'FontSize', 11);
title('regularization\_sdp: delta > rho\_rough × 1.01', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
xticks(1:length(attack_sizes));
xticklabels(arrayfun(@(a) sprintf('%.0e', a), attack_sizes, 'UniformOutput', false));
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 10);

% === グラフ4: regularization_sdp + 判定2 ===
subplot(2, 2, 4);
hold on;
for p_idx = 1:length(param_list)
    success_rate = results.reg_success_rate_hinf(p_idx, :);
    plot(1:length(attack_sizes), success_rate, ...
        'Color', colors(p_idx, :), ...
        'LineWidth', 2, ...
        'Marker', markers{p_idx}, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', colors(p_idx, :), ...
        'DisplayName', sprintf('epsilon = %.0e', param_list(p_idx)));
end
hold off;
xlabel('攻撃の大きさ', 'FontSize', 11);
ylabel('攻撃成功率 (%)', 'FontSize', 11);
title('regularization\_sdp: hinf > (1/delta) × 1.01', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
xticks(1:length(attack_sizes));
xticklabels(arrayfun(@(a) sprintf('%.0e', a), attack_sizes, 'UniformOutput', false));
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 10);

% 全体タイトル
sgtitle('実験2: ロバスト性の評価（攻撃あり）', 'FontSize', 14, 'FontWeight', 'bold');

% 保存
saveas(gcf, 'experiment2_robustness.png');
saveas(gcf, 'experiment2_robustness.fig');
fprintf('実験2のグラフを保存しました: experiment2_robustness.png\n');

end

