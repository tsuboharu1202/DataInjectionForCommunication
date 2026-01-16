classdef SaveConfig
    % SaveConfig - 結果保存の設定
    %
    % 使用例:
    %   savedir = cfg.SaveConfig.getResultDir('MyExperiment');
    %   filename = cfg.SaveConfig.generateFilename('results', 'mat');
    
    properties (Constant)
        % 結果保存のルートディレクトリ（プロジェクトルートからの相対パス）
        RESULT_ROOT = 'Experiment'
        
        % タイムスタンプ形式
        TIMESTAMP_FORMAT = 'yyyyMMdd_HHmmss'
        
        % 図の保存形式
        FIGURE_FORMATS = {'png', 'fig'}
        
        % デフォルトの図サイズ [width, height]
        FIGURE_SIZE = [800, 600]
        
        % 図のDPI
        FIGURE_DPI = 150
    end
    
    methods (Static)
        function dirpath = getResultDir(experiment_name)
            % 実験結果ディレクトリのパスを取得
            %   存在しない場合は作成する
            dirpath = fullfile(cfg.SaveConfig.RESULT_ROOT, experiment_name, 'results');
            if ~exist(dirpath, 'dir')
                mkdir(dirpath);
            end
        end
        
        function dirpath = getFigureDir(experiment_name)
            % 図の保存ディレクトリのパスを取得
            dirpath = fullfile(cfg.SaveConfig.getResultDir(experiment_name), 'figures');
            if ~exist(dirpath, 'dir')
                mkdir(dirpath);
            end
        end
        
        function filename = generateFilename(prefix, extension)
            % タイムスタンプ付きファイル名を生成
            %   例: 'results_20260116_123456.mat'
            timestamp = datestr(now, cfg.SaveConfig.TIMESTAMP_FORMAT);
            filename = sprintf('%s_%s.%s', prefix, timestamp, extension);
        end
        
        function filepath = getFullPath(experiment_name, filename)
            % 完全なファイルパスを取得
            filepath = fullfile(cfg.SaveConfig.getResultDir(experiment_name), filename);
        end
        
        function saveFigure(fig, experiment_name, filename_without_ext)
            % 図を保存（png と fig の両方）
            figdir = cfg.SaveConfig.getFigureDir(experiment_name);
            for i = 1:length(cfg.SaveConfig.FIGURE_FORMATS)
                fmt = cfg.SaveConfig.FIGURE_FORMATS{i};
                filepath = fullfile(figdir, sprintf('%s.%s', filename_without_ext, fmt));
                if strcmp(fmt, 'fig')
                    savefig(fig, filepath);
                else
                    print(fig, filepath, sprintf('-d%s', fmt), ...
                        sprintf('-r%d', cfg.SaveConfig.FIGURE_DPI));
                end
            end
            fprintf('図を保存しました: %s\n', fullfile(figdir, filename_without_ext));
        end
        
        function saveResults(data, experiment_name, filename)
            % 結果をmatファイルで保存
            filepath = cfg.SaveConfig.getFullPath(experiment_name, filename);
            save(filepath, '-struct', 'data');
            fprintf('結果を保存しました: %s\n', filepath);
        end
    end
end

