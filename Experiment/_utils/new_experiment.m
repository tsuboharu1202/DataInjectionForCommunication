function new_experiment(experiment_name)
% new_experiment - テンプレートから新しい実験フォルダを作成
%
% 使い方:
%   new_experiment('MyExperiment')
%   new_experiment('TradeoffAnalysis_v2')
%
% 引数:
%   experiment_name: 実験名（フォルダ名になります）

if nargin < 1 || isempty(experiment_name)
    experiment_name = sprintf('Experiment_%s', datestr(now, 'yyyymmdd_HHMMSS'));
end

% パス設定
script_dir = fileparts(mfilename('fullpath'));
template_dir = fullfile(script_dir, '_template');
new_dir = fullfile(script_dir, experiment_name);

% テンプレート確認
if ~exist(template_dir, 'dir')
    error('テンプレートフォルダが見つかりません: %s', template_dir);
end

% 既存チェック
if exist(new_dir, 'dir')
    error('同名のフォルダが既に存在します: %s', new_dir);
end

% コピー
fprintf('テンプレートをコピー中...\n');
copyfile(template_dir, new_dir);

% config.m の実験名を更新
config_file = fullfile(new_dir, 'config.m');
config_content = fileread(config_file);
config_content = strrep(config_content, 'Experiment_Template', experiment_name);
fid = fopen(config_file, 'w');
fprintf(fid, '%s', config_content);
fclose(fid);

fprintf('\n========================================\n');
fprintf('新しい実験フォルダを作成しました:\n');
fprintf('  %s\n', new_dir);
fprintf('========================================\n\n');

fprintf('次のステップ:\n');
fprintf('  1. cd %s\n', new_dir);
fprintf('  2. config.m を編集\n');
fprintf('  3. run_experiment.m に実験ロジックを実装\n');
fprintf('  4. run_experiment を実行\n\n');

% フォルダに移動するか確認
fprintf('フォルダに移動しますか？ (y/n): ');
answer = input('', 's');
if strcmpi(answer, 'y')
    cd(new_dir);
    fprintf('移動しました: %s\n', pwd);
end

end

