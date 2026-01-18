=====================================
実験テンプレート使用方法
=====================================

【新しい実験を始めるとき】

1. テンプレートフォルダをコピー
   >> copyfile('Experiment/_template', 'Experiment/MyExperiment_20260116');

2. コピーしたフォルダに移動
   >> cd Experiment/MyExperiment_20260116

3. config.m を編集して実験パラメータを設定

4. run_experiment.m に実験ロジックを実装

5. plot_results.m に描画ロジックを実装

6. 実験を実行
   >> run_experiment


【ファイル構成】

_template/
├── config.m           # 実験設定（パラメータを編集）
├── run_experiment.m   # 実験実行スクリプト
├── plot_results.m     # 結果描画スクリプト
├── README.txt         # このファイル
└── results/           # 結果保存先
    └── figures/       # 図の保存先


【出力ファイル】

results/
├── data.mat           # 実験結果（results, cfg, rho_rough）
├── log_YYYYMMDD_HHMMSS.txt  # 実験ログ
└── figures/
    ├── ExperimentName.png
    └── ExperimentName.fig


【Tips】

- 乱数シードを固定すると再現性が確保できます
  config.m: cfg.random_seed = 42;

- 既存の実験をベースに新しい実験を作る場合は、
  そのフォルダをコピーして config.m を編集

- 結果を比較するときは、各実験の data.mat を読み込んで
  比較用のスクリプトを別途作成

