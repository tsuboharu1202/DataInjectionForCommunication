function plot_delta_hinf_ratio(csv_file)
% plot_delta_hinf_ratio - (1/delta)/hinf_normとgammaの関係を描画
%
% 使用方法:
%   plot_delta_hinf_ratio()  % デフォルトでresults.csvを読み込む
%   plot_delta_hinf_ratio('results.csv')  % ファイルを指定
%
% 描画内容:
%   - 横軸: gamma (対数スケール)
%   - 縦軸: (1/delta)/hinf_norm
%   - epsilonごとに異なる線（3本の折れ線グラフ）

% ============================================
% 1. データの読み込み
% ============================================
if nargin < 1 || isempty(csv_file)
    csv_file = fullfile(fileparts(mfilename('fullpath')), 'results.csv');
end

fprintf('CSVファイルを読み込み中: %s\n', csv_file);
data = readtable(csv_file);

% infeasibleなデータを除外
data = data(data.infeasible == 0, :);

% NaN値を除外
data = data(~isnan(data.delta) & ~isnan(data.hinf_norm), :);

% (1/delta)/hinf_normを計算
data.ratio = (1 ./ data.delta) ./ data.hinf_norm;

% ============================================
% 2. データの集計（各epsilon、各gammaでの平均値）
% ============================================
% 数値として確実に扱う
data.epsilon = double(data.epsilon);
data.gamma = double(data.gamma);

epsilon_list = unique(data.epsilon);
gamma_list = unique(data.gamma);
gamma_list = sort(gamma_list);  % ソート

% デバッグ: gamma_listの値を確認
fprintf('Gamma値: ');
fprintf('%.0e ', gamma_list);
fprintf('\n');

n_epsilon = length(epsilon_list);
n_gamma = length(gamma_list);

% 各epsilon、各gammaでの平均ratioを計算
mean_ratio = zeros(n_epsilon, n_gamma);

for eps_idx = 1:n_epsilon
    epsilon = epsilon_list(eps_idx);
    eps_data = data(abs(data.epsilon - epsilon) < 1e-10, :);
    
    for g_idx = 1:n_gamma
        gamma = gamma_list(g_idx);
        % 浮動小数点数の比較は許容誤差を使用
        gamma_data = eps_data(abs(eps_data.gamma - gamma) < 1e-10, :);
        
        if ~isempty(gamma_data)
            mean_ratio(eps_idx, g_idx) = mean(gamma_data.ratio);
        else
            mean_ratio(eps_idx, g_idx) = NaN;
        end
    end
end

% ============================================
% 3. 描画
% ============================================
figure('Name', 'Delta/H∞ノルム vs Gamma', 'Position', [100, 100, 1000, 600]);
hold on;

colors = lines(n_epsilon);
markers = {'o', 's', '^'};

for eps_idx = 1:n_epsilon
    epsilon = epsilon_list(eps_idx);
    
    % 有効なデータのみを抽出
    valid_idx = ~isnan(mean_ratio(eps_idx, :));
    gamma_vals = gamma_list(valid_idx);
    ratio_vals = mean_ratio(eps_idx, valid_idx);
    
    % 折れ線グラフを描画
    semilogx(gamma_vals, ratio_vals, ...
        'Color', colors(eps_idx, :), ...
        'LineWidth', 2, ...
        'Marker', markers{mod(eps_idx-1, length(markers))+1}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(eps_idx, :), ...
        'DisplayName', sprintf('ε = %.0e', epsilon));
end

xlabel('Gamma (対数スケール)', 'FontSize', 12);
ylabel('(1/δ) / H∞ノルム', 'FontSize', 12);
title('Gamma値と(1/δ)/H∞ノルムの関係', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
grid on;
grid minor;

% 横軸の範囲を明示的に設定（gammaの最小値と最大値）
% semilogxを使っているので、対数スケールで表示される
if ~isempty(gamma_list) && all(~isnan(gamma_list)) && all(gamma_list > 0)
    xlim([min(gamma_list), max(gamma_list)]);
    % 対数スケールの目盛りを設定
    set(gca, 'XScale', 'log');
end

% 縦軸の範囲を調整
ylim_vals = ylim;
ylim([0, ylim_vals(2) * 1.1]);

hold off;

% ============================================
% 4. 統計情報の表示
% ============================================
fprintf('\n=== 統計情報 ===\n\n');
for eps_idx = 1:n_epsilon
    epsilon = epsilon_list(eps_idx);
    valid_idx = ~isnan(mean_ratio(eps_idx, :));
    fprintf('Epsilon = %.0e:\n', epsilon);
    fprintf('  Gamma範囲: %.0e ～ %.0e\n', min(gamma_list(valid_idx)), max(gamma_list(valid_idx)));
    fprintf('  Ratio範囲: %.4f ～ %.4f\n', ...
        min(mean_ratio(eps_idx, valid_idx)), max(mean_ratio(eps_idx, valid_idx)));
    fprintf('  Ratio平均: %.4f\n', mean(mean_ratio(eps_idx, valid_idx), 'omitnan'));
    fprintf('\n');
end

end

