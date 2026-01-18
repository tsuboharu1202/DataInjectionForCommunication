=====================================
攻撃データ生成実験
=====================================

【実験の目的】

同一のシステム設定を用いて、次回以降のシミュレーションに使用するための
攻撃データセットを生成する。

【実装概要】

1. データ生成
   - 20セットの (U, X, Z) を生成
   - データ生成方法: datasim.simulate_openloop_stable() を使用
   - 各データセットに対して元のデータでSDPを解き、deltaを記録

2. 攻撃実行
   - 各データセットに対して、複数のepsilon値でIDGSM攻撃を実行
   - 各epsilonに対して、negative と positive の両方向で攻撃
   - 攻撃パラメータ:
     * step_size = epsilon * step_size_coefficient
     * max_iteration = step_max_iter
     * normalize_grad = false（勾配正規化なし）
     * save_history = false（履歴保存なし）

3. 結果保存
   - 元のデータ: U, X, Z, delta, K
   - 攻撃後のデータ: U_adv, X_adv, Z_adv, delta, K
   - 各攻撃の成功/失敗フラグ

【設定ファイル (config.m)】

主要パラメータ:
- cfg.A, cfg.B: システム行列
- cfg.T: サンプル数
- cfg.gamma: 正則化パラメータ（固定値、デフォルト: 1e3）
- cfg.attack_eps_list: 攻撃サイズのリスト（例: [1e-4, 1e-3, 1e-2, 1e-1, 1e0]）
- cfg.step_size_coefficient: ステップサイズ係数（デフォルト: 1e-1）
- cfg.step_max_iter: 最大反復回数（デフォルト: 30）
- cfg.trial: データセット数（デフォルト: 20）
- cfg.save_results: 結果を保存するか（デフォルト: false）

【実行方法】

1. config.m を編集して実験パラメータを設定
   - 特に cfg.save_results = true に設定すると結果が保存されます

2. run_experiment.m を実行
   >> run_experiment

【出力ファイル】

results/
├── attack_data_YYYYMMDD_HHMMSS.mat  # 実験結果（cfg.save_results=trueの場合）
│   └── 構造体:
│       - results.datasets(trial).original: 元のデータ
│       - results.datasets(trial).attacks(attack_idx): 攻撃結果
│       - cfg: 実験設定
│       - rho_rough: 理論的下界
└── log_YYYYMMDD_HHMMSS.txt          # 実験ログ（常に保存）

【データ構造】

results.datasets(trial).original:
  - U, X, Z: 元のデータ
  - delta: 元のdelta値
  - K: 元のコントローラゲイン

results.datasets(trial).attacks(attack_idx):
  - epsilon: 攻撃サイズ
  - step_size: ステップサイズ
  - direction: 'negative' または 'positive'
  - U_adv, X_adv, Z_adv: 攻撃後のデータ
  - delta: 攻撃後のdelta値
  - K: 攻撃後のコントローラゲイン
  - success: 攻撃が成功したか（SDPが解けたか）

【注意事項】

- エラーが発生した攻撃はデータに記録されません（スキップ）
- 攻撃履歴（各反復のdX, dZ, dU）は保存されません（データサイズを抑えるため）
- 結果ファイルは .gitignore によりGit管理から除外されます
- ログファイルは常に保存されます

【統計情報】

実験終了時に以下の統計情報が表示されます:
- 総攻撃試行数
- 成功数と成功率
- Delta変化の統計（平均、最小、最大）
