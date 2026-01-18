# Data-Driven Control with SDP and Data Injection Attack

データ駆動制御における半正定値計画問題（SDP）の解法と、データ注入攻撃のシミュレーションを行うMATLABプロジェクトです。

## 概要

このプロジェクトは、以下の機能を提供します：

- **データ駆動制御のSDP解法**: 論文に基づくSDP問題を解き、制御ゲインを計算
- **Implicit Differentiation**: KKT条件を用いた勾配計算
- **データ注入攻撃のシミュレーション**: データを改ざんして制御性能を劣化させる攻撃の実行
- **H∞ノルム計算**: 閉ループシステムのH∞ノルムを計算

## ディレクトリ構造

```
ForCommunication/
├── +cfg/                       # 設定クラス
│   ├── Const.m                # 定数定義（サンプル数、ステップサイズなど）
│   ├── AttackType.m           # 攻撃タイプの定義
│   ├── System.m               # デフォルトシステム設定
│   ├── SDP.m                  # SDP関連設定
│   └── SaveConfig.m           # 結果保存設定
├── core/                       # コア機能
│   ├── +datasim/              # データ生成・シミュレーション
│   │   ├── make_lti.m         # LTIシステム生成
│   │   ├── simulate_openloop_stable.m  # 安定な開ループシミュレーション
│   │   ├── simulate_openloop.m # 開ループシミュレーション
│   │   └── SystemData.m       # システムデータクラス
│   ├── +input/                # 入力生成
│   │   └── make_inputU.m      # 入力信号生成
│   └── +helper/               # ヘルパー関数
│       └── hinfnorm_AK.m      # H∞ノルム計算
├── methods/                    # SDP解法
│   ├── +baseline/             # ベースライン手法（論文のSDP）
│   │   ├── solve_sdp.m        # SDP問題を解く
│   │   ├── build_lmi_blocks.m # LMIブロック構築
│   │   └── +gradient/         # 勾配計算
│   ├── +proposed/             # 提案手法（正則化付きSDP）
│   │   ├── solve_sdp.m        # 正則化付きSDP
│   │   ├── solve_sdp_with_robust.m  # ロバスト制約付きSDP
│   │   └── +gradient/         # 勾配計算
│   └── +takakiSenpai/         # 量子化SDP
├── attack/                     # データ注入攻撃
│   ├── +algorithms/           # 攻撃アルゴリズム
│   │   ├── execute_attack.m   # 攻撃実行（統一インターフェース）
│   │   ├── dgsm_delta.m       # 直接勾配符号法
│   │   ├── idgsm_delta.m      # 反復勾配符号法
│   │   ├── make_data_adv.m    # 攻撃データ生成
│   │   └── projector.m        # ノルム制約投影
│   └── +common/               # 共通処理
│       └── calc_grad.m        # 勾配計算
├── demos/                      # デモスクリプト
│   ├── demo_sdp.m             # SDP解法のデモ
│   ├── demo_attack.m          # 攻撃シミュレーションのデモ
│   └── demo_implicit_grad.m   # 勾配計算のデモ
├── Experiment/                 # 実験
│   ├── _template/             # 実験テンプレート
│   ├── MethodComparison/      # 手法比較実験
│   └── TradeoffAnalysis/      # トレードオフ分析
├── doc/                        # ドキュメント
├── _archive/                   # アーカイブ（旧コード）
└── startup.m                   # パス設定スクリプト
```

## セットアップ

1. **MATLABの起動**
   - MATLAB R2023a以降を推奨

2. **必要なツールボックス**
   - YALMIP（SDPソルバー）
   - MOSEK（または他のSDPソルバー）
   - Control System Toolbox（H∞ノルム計算用）

3. **パスの設定**
   ```matlab
   % プロジェクトルートで実行
   startup
   ```

## 基本的な使用方法

### 1. SDP解法

```matlab
% システム設定
A = cfg.System.A;
B = cfg.System.B;
[n, m] = cfg.System.getDimensions();
T = 20;  % サンプル数

% データ生成
V = make_inputU(m, T);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);

% SystemData作成（Phi自動生成）
data = datasim.SystemData.create(A, B, X, Z, U);

% ベースラインSDP
[sol, K, ~, ~, ~] = baseline.solve_sdp(data);

% 提案手法SDP（gamma必須）
gamma = 1e3;
[sol, K, delta, ~, ~] = proposed.solve_sdp(data, gamma);
```

### 2. データ注入攻撃

```matlab
% 攻撃オプション（gamma必須）
opts = struct();
opts.gamma = 1e3;
opts.epsilon = 1e-3;
opts.direction = 'positive';  % deltaを大きくする
opts.save_history = true;

% 攻撃実行
[X_adv, Z_adv, U_adv, history] = algorithms.execute_attack(...
    data, cfg.AttackType.IMPLICIT_IDGSM_DELTA, opts);
```

### 3. ノイズ付きシミュレーション

```matlab
% オプションでノイズを指定
sim_opts = struct('noise_std', 0.01);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V, [], sim_opts);
```

## 主要な関数

### SDP解法

| 関数 | 説明 |
|-----|------|
| `baseline.solve_sdp(data)` | 論文のSDP解法 |
| `proposed.solve_sdp(data, gamma)` | 正則化付きSDP（gamma必須） |
| `proposed.solve_sdp_with_robust(data, epsilon, gamma)` | ロバスト制約付きSDP |

### 攻撃

| 関数 | 説明 |
|-----|------|
| `algorithms.execute_attack(data, method, opts)` | 統一攻撃インターフェース |
| `algorithms.dgsm_delta(data, opts)` | 直接勾配符号法 |
| `algorithms.idgsm_delta(data, opts)` | 反復勾配符号法 |

### ヘルパー

| 関数 | 説明 |
|-----|------|
| `helper.hinfnorm_AK(A, B, K)` | H∞ノルム計算 |
| `cfg.System.calcRhoRough(A)` | 理論的下界計算 |

## 設定

### 必須パラメータ（デフォルト値なし）

全ての実験・攻撃で以下のパラメータを**明示的に指定**してください：

| パラメータ | 説明 | 例 |
|-----------|------|-----|
| `T` | サンプル数 | 20 |
| `gamma` | 正則化パラメータ | 1e3 |
| `epsilon` | 攻撃強度 | 1e-3 |
| `alpha` | IDGSMステップサイズ | 1e-4 |
| `max_iteration` | IDGSM最大反復回数 | 30 |
| `direction` | 攻撃方向 | 'positive' or 'negative' |

### cfg.System（システム設定）

デフォルトのシステム行列（A, B）や、Phi行列の生成メソッドを提供。

### cfg.SDP（SDP設定）

| 定数 | 説明 |
|-----|------|
| `TOLERANCE_MARGIN` | 判定マージン（1.01 = 1%） |
| `SOLVER_VERBOSE` | ソルバー出力（0: 非表示） |

**注意**: `gamma` はデフォルト値を提供しません。実験ごとに明示的に指定してください。

## エラー処理

- SDPが解けない場合、`error`を出して処理を中断します（適当な値を返しません）
- 必須パラメータ（gamma等）が指定されていない場合もエラーを出します
- エラーメッセージは日本語で統一されています

## ライセンス

（必要に応じてライセンス情報を追加）

## 作者

（必要に応じて作者情報を追加）
