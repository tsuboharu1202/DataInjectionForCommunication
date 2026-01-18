% plot_results.m - 結果描画テンプレート
%
% このスクリプトは run_experiment.m から呼び出されるか、
% 単独で実行できます（その場合は results/data.mat を読み込みます）
%

%% ========== 結果読み込み（単独実行時） ==========
if ~exist('results', 'var')
    fprintf('結果ファイルを読み込み中...\n');
    load('results/data.mat', 'results', 'cfg', 'rho_rough');
end

%% ========== 描画設定 ==========
% カラー
colors = struct();
colors.baseline = [0.2, 0.4, 0.8];   % 青
colors.proposed = [0.8, 0.3, 0.3];   % 赤
colors.reference = [0.3, 0.7, 0.3];  % 緑

% マーカー
markers = {'o', 's', 'd', '^', 'v'};

% フォントサイズ
fontSize = struct();
fontSize.title = 14;
fontSize.label = 12;
fontSize.tick = 11;
fontSize.legend = 11;

%% ========== グラフ1: サンプル（編集してください） ==========
figure('Name', cfg.experiment_name, 'Position', [100, 100, 800, 600]);

% TODO: ここに描画ロジックを実装
% 
% 例:
% hold on;
% plot(cfg.param_list, results.delta_baseline, 'b-o', 'LineWidth', 2, 'DisplayName', 'baseline');
% plot(cfg.param_list, results.delta_proposed, 'r-s', 'LineWidth', 2, 'DisplayName', 'proposed');
% yline(rho_rough, '--', 'Color', colors.reference, 'LineWidth', 2, 'DisplayName', 'rho\_rough');
% hold off;
% 
% xlabel('パラメータ', 'FontSize', fontSize.label);
% ylabel('delta', 'FontSize', fontSize.label);
% title(cfg.experiment_name, 'FontSize', fontSize.title);
% legend('Location', 'best', 'FontSize', fontSize.legend);
% set(gca, 'XScale', 'log');
% grid on;

% プレースホルダー
text(0.5, 0.5, '描画ロジックを実装してください', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 16);
title(sprintf('%s - 結果', cfg.experiment_name), 'FontSize', fontSize.title);

%% ========== 保存 ==========
if cfg.save_figures
    fig_basename = sprintf('results/figures/%s', cfg.experiment_name);
    
    % PNG
    saveas(gcf, [fig_basename '.png']);
    fprintf('図を保存しました: %s.png\n', fig_basename);
    
    % FIG（MATLABで再編集可能）
    saveas(gcf, [fig_basename '.fig']);
    
    % PDF（論文用）
    if strcmp(cfg.figure_format, 'pdf')
        exportgraphics(gcf, [fig_basename '.pdf'], 'ContentType', 'vector');
        fprintf('図を保存しました: %s.pdf\n', fig_basename);
    end
end

fprintf('描画完了\n');

