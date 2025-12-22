clear; clc; close all;

% === 保存先 ===
thisFileDir = fileparts(mfilename('fullpath'));   % ライブスクリプトなら pwd
datasetDir  = fullfile(thisFileDir, 'DataSet');
if ~exist(datasetDir, 'dir'); mkdir(datasetDir); end

% === 実行メタ ===
dateTag = datestr(now, 'yyyymmdd_HHMM');
csvFile = fullfile(datasetDir, sprintf('attack_log_%s.csv', dateTag));
run_id  = datestr(now,'yyyymmdd_HHMMSS_FFF');
ULIM    = cfg.Const.ATTACKER_UPPERLIMIT;

% === ループ設定 ===
Nlog         = 20;
epsilon_iter = 4;
noise_iter   = 4;

% 基準モデルのパラメータ（ここからノイズを作る）
eps_base   = 0.05;
noise_base = 1e-5;

% 各設計方策用の ε, bound
epsilon_vals = 0.005*(1:epsilon_iter);            % 例: 0.02, 0.04
noise_vals   = 10.^(-(4 + (1:noise_iter)));      % 例: 1e-6, 1e-7

% === ヘッダ作成（初回のみ）===
if ~isfile(csvFile)
    header = ["run_id","index","timestamp","eig1" ,"eig2" ,"eig3", ...
              "rho_base","spectral_base"];       % 基準モデル（ノイズなし）

    header = [header, ...
        sprintf("rho_clean"), ...
        sprintf("spec_clean"), ...
        sprintf("rho_noisy"), ...
        sprintf("spec_noisy")];

    % 各設計方策ごとに「ノイズなし」「ノイズあり」の rho, spectral を横に並べる
    for ei = 1:epsilon_iter
        for bi = 1:noise_iter
            ev = epsilon_vals(ei);
            nv = noise_vals(bi);
            header = [header, ...
                sprintf("rho_clean_e%.3g_b%.0e",   ev, nv), ...
                sprintf("spec_clean_e%.3g_b%.0e",  ev, nv), ...
                sprintf("rho_noisy_e%.3g_b%.0e",   ev, nv), ...
                sprintf("spec_noisy_e%.3g_b%.0e",  ev, nv)];
        end
    end

    header = [header, "ATTACKER_UPPERLIMIT"];

    fid0 = fopen(csvFile,'w');
    fprintf(fid0, "%s\n", strjoin(header, ","));
    fclose(fid0);
end

% === 追記オープン ===
fid = fopen(csvFile,'a');

for index = 1:Nlog
    ts = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

    % --- 基準モデル：これを元にノイズ生成 ---
    sd_ori = make_sd(eps_base, noise_base);   % あなたの make_sd(eps, bound)
    disp('-------------');
    [spectral_base, rho_base] = sd_ori.calc_spectral_rho;

    % 同じ sd_ori からノイズ付きデータを作る（X_adv,Z_adv,U_adv）
    eig_value = eig(sd_ori.A);

    % 1行ぶんのセルを準備
    row_cells = {run_id, num2str(index), ts, ...
                 sprintf('%.16g', eig_value(1)), ...
                 sprintf('%.16g', eig_value(2)), ...
                 sprintf('%.16g', eig_value(3)), ...
                 sprintf('%.16g', rho_base), ...
                 sprintf('%.16g', spectral_base)};

    disp('A');disp(sd_ori.A);
    disp('B');disp(sd_ori.B);
    [X_adv, Z_adv, U_adv] = attack.dgsm(sd_ori);
    visualize.plot_data(sd_ori.X,X_adv);

    sd_copy = sd_ori;
    [sol, K, ~, ~, ~] = solve_quantized_sdp_perturbation(sd_copy);
    rho_temp = sol.rho;
    spectral_temp = max(abs(eig(sd_copy.A + sd_copy.B*K))); 
    row_cells{end+1} = sprintf('%.16g', rho_temp);
    row_cells{end+1} = sprintf('%.16g', spectral_temp);

    sd_target_noisy = sd_copy.withXZU(X_adv, Z_adv, U_adv);
    [sol_noise, K_noise, ~, ~, ~] = solve_quantized_sdp_perturbation(sd_target_noisy);
    rho_noise = sol_noise.rho;
    spectral_noise = max(abs(eig(sd_copy.A + sd_copy.B*K_noise))); 
    row_cells{end+1} = sprintf('%.16g', rho_noise);
    row_cells{end+1} = sprintf('%.16g', spectral_noise);




    % --- 各設計方策 (epsilon, noise) ごとにノイズなし/ありを計算 ---
    for ei = 1:epsilon_iter
        for bi = 1:noise_iter
            eps_val   = epsilon_vals(ei);
            noise_val = noise_vals(bi);

            % この方策の「理論値」（ノイズなし）
            sd_target = change_param(sd_ori, eps_val, noise_val);   % ← パラメータだけ変える関数（想定）
            [spec_clean, rho_clean] = sd_target.calc_spectral_rho;

            % 同じノイズ (X_adv,Z_adv,U_adv) を入れたとき
            sd_target_noisy = sd_target.withXZU(X_adv, Z_adv, U_adv);
            [spec_noisy, rho_noisy] = sd_target_noisy.calc_spectral_rho;

            % 横に追加
            row_cells{end+1} = sprintf('%.16g', rho_clean);
            row_cells{end+1} = sprintf('%.16g', spec_clean);
            row_cells{end+1} = sprintf('%.16g', rho_noisy);
            row_cells{end+1} = sprintf('%.16g', spec_noisy);
        end
    end

    % 最後に上限
    row_cells{end+1} = sprintf('%.16g', ULIM);

    % === CSVへ1行出力 ===
    fprintf(fid, "%s\n", strjoin(row_cells, ","));
end

fclose(fid);

disp("✅ CSV 追記完了: " + csvFile);
