function extract_results_from_text(text_file, output_file)
% extract_results_from_text - テキスト出力から実験結果を抽出してCSV/テーブルに保存
%
% 使用方法:
%   extract_results_from_text()  % 対話的にファイルを選択
%   extract_results_from_text('output.txt', 'results.csv')  % ファイルを指定
%
% 抽出する情報:
%   - epsilon: 攻撃の大きさ
%   - trial: 試行番号
%   - gamma: 正則化パラメータ
%   - delta: SDPの解のdelta値
%   - hinf_norm: H∞ノルム
%   - success: 攻撃成功フラグ (1/delta < hinf_norm)
%   - infeasible: Infeasibleフラグ (エラーが発生した場合1)

% ============================================
% 1. ファイルの読み込み
% ============================================
if nargin < 1 || isempty(text_file)
    [filename, pathname] = uigetfile({'*.txt;*.log', 'テキストファイル'}, ...
        'テキスト出力ファイルを選択');
    if isequal(filename, 0)
        error('ファイルが選択されませんでした。');
    end
    text_file = fullfile(pathname, filename);
end

fprintf('テキストファイルを読み込み中: %s\n', text_file);
fid = fopen(text_file, 'r');
if fid == -1
    error('ファイルを開けませんでした: %s', text_file);
end

text_content = fread(fid, '*char')';
fclose(fid);

% ============================================
% 2. データの抽出
% ============================================
% 結果を格納する配列
results = [];

% 行ごとに分割
lines = strsplit(text_content, '\n');

% 現在の試行情報を保持
current_epsilon = [];
current_trial = [];

% パターンマッチング用の正規表現
epsilon_pattern = 'Epsilon=([0-9.e+-]+)';
trial_pattern = 'Trial=(\d+)/\d+';
gamma_pattern = 'gamma\s*=\s*([0-9.e+-]+)';
delta_pattern = 'delta=([0-9.e+-]+)';
hinf_pattern = 'hinf_norm=([0-9.e+-]+)';
success_pattern = 'success=(\d+)';
error_pattern = 'エラー|infeasible|Numerical problems';

for i = 1:length(lines)
    line = strtrim(lines{i});
    
    % EpsilonとTrialの情報を抽出
    epsilon_match = regexp(line, epsilon_pattern, 'tokens');
    trial_match = regexp(line, 'Trial=(\d+)/\d+', 'tokens');
    
    if ~isempty(epsilon_match)
        current_epsilon = str2double(epsilon_match{1}{1});
    end
    
    if ~isempty(trial_match)
        current_trial = str2double(trial_match{1}{1});
    end
    
    % gamma行の処理
    if contains(line, 'gamma =') || contains(line, 'gamma=')
        % gamma値を抽出
        gamma_match = regexp(line, gamma_pattern, 'tokens');
        if isempty(gamma_match)
            continue;
        end
        gamma_val = str2double(gamma_match{1}{1});
        
        % エラーチェック
        is_error = ~isempty(regexp(line, error_pattern, 'once'));
        
        if is_error
            % エラーの場合
            results = [results; struct(...
                'epsilon', current_epsilon, ...
                'trial', current_trial, ...
                'gamma', gamma_val, ...
                'delta', NaN, ...
                'hinf_norm', NaN, ...
                'success', 0, ...
                'infeasible', 1)];
        else
            % 正常な結果の場合
            delta_match = regexp(line, delta_pattern, 'tokens');
            hinf_match = regexp(line, hinf_pattern, 'tokens');
            success_match = regexp(line, success_pattern, 'tokens');
            
            if ~isempty(delta_match) && ~isempty(hinf_match)
                delta_val = str2double(delta_match{1}{1});
                hinf_val = str2double(hinf_match{1}{1});
                
                if ~isempty(success_match)
                    success_val = str2double(success_match{1}{1});
                else
                    % successが明示されていない場合は計算
                    success_val = (1/delta_val < hinf_val);
                end
                
                results = [results; struct(...
                    'epsilon', current_epsilon, ...
                    'trial', current_trial, ...
                    'gamma', gamma_val, ...
                    'delta', delta_val, ...
                    'hinf_norm', hinf_val, ...
                    'success', success_val, ...
                    'infeasible', 0)];
            end
        end
    end
end

% ============================================
% 3. テーブルに変換
% ============================================
if isempty(results)
    warning('結果が見つかりませんでした。');
    return;
end

results_table = struct2table(results);

% ============================================
% 4. 結果の表示
% ============================================
fprintf('\n=== 抽出結果 ===\n\n');
fprintf('抽出されたレコード数: %d\n', height(results_table));
fprintf('\n最初の10行:\n');
disp(results_table(1:min(10, height(results_table)), :));

% ============================================
% 5. ファイルに保存
% ============================================
if nargin < 2 || isempty(output_file)
    [pathname, ~, ~] = fileparts(text_file);
    output_file = fullfile(pathname, 'extracted_results.csv');
end

% CSV形式で保存
fprintf('\n結果を保存中: %s\n', output_file);
writetable(results_table, output_file);
fprintf('保存完了！\n');

% MATファイルも保存（MATLABで扱いやすい）
[pathname, filename, ~] = fileparts(output_file);
mat_file = fullfile(pathname, [filename, '.mat']);
save(mat_file, 'results_table', 'results');
fprintf('MATファイルも保存: %s\n', mat_file);

% サマリーを表示
fprintf('\n=== サマリー ===\n');
fprintf('Epsilon値: ');
fprintf('%.0e ', unique(results_table.epsilon));
fprintf('\n');
fprintf('Gamma値: ');
fprintf('%.0e ', unique(results_table.gamma));
fprintf('\n');
fprintf('総試行数: %d\n', length(unique(results_table.trial)));
fprintf('Infeasible数: %d\n', sum(results_table.infeasible == 1));
fprintf('攻撃成功数: %d\n', sum(results_table.success == 1));

end


