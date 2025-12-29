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
├── +cfg/                    # 設定クラス
│   ├── Const.m             # 定数定義（サンプル数、ステップサイズなど）
│   └── AttackType.m        # 攻撃タイプの定義
├── basic_src/               # 基本機能
│   ├── +datasim/           # データ生成・シミュレーション
│   │   ├── make_lti.m      # LTIシステム生成
│   │   ├── simulate_openloop_stable.m  # 安定な開ループシミュレーション
│   │   └── SystemData.m    # システムデータクラス
│   └── +visualize/         # 可視化機能
├── com_src/                 # 主要機能
│   ├── +original_thesis/   # 論文のSDP解法
│   │   ├── solve_sdp.m     # SDP問題を解く
│   │   └── build_lmi_blocks.m  # LMIブロック構築
│   ├── +regularization_sdp/  # 正則化付きSDP解法
│   │   └── solve_sdp.m
│   ├── +implicit/           # Implicit Differentiation
│   │   ├── dtDelta_dD.m    # 勾配計算メイン関数
│   │   └── +G_grad/        # KKT条件の各項（G1-G15）
│   ├── +attack/             # データ注入攻撃
│   │   ├── execute_attack.m # 攻撃実行
│   │   └── calc_grad.m      # 勾配計算
│   └── +helper/             # ヘルパー関数
│       └── hinfnorm_AK.m   # H∞ノルム計算
├── scripts/                 # デモスクリプト
│   ├── demo_sdp.m          # SDP解法のデモ
│   ├── demo_implicit_grad.m # 勾配計算のデモ
│   └── demo_attack.m       # 攻撃シミュレーションのデモ
├── doc/                     # ドキュメント
│   └── kkt_implicit_differentiation.tex  # KKT条件とImplicit Differentiationの理論
├── Experiment/              # 実験用スクリプトとデータ
└── startup.m                # パス設定スクリプト
```

## セットアップ

1. **MATLABの起動**
   - MATLAB R2025a以降を推奨

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

### 1. SDP解法のデモ

```matlab
% scripts/demo_sdp.m を実行
demo_sdp
```

このスクリプトは：
- LTIシステムを生成
- データを生成
- 元の論文のSDP解法と正則化付きSDP解法を実行
- H∞ノルムを計算して比較

### 2. Implicit Differentiationによる勾配計算

```matlab
% scripts/demo_implicit_grad.m を実行
demo_implicit_grad
```

このスクリプトは：
- SDPを解く
- KKT条件を用いて勾配を計算
- 行列のランクを確認

### 3. データ注入攻撃のシミュレーション

```matlab
% scripts/demo_attack.m を実行
demo_attack
```

このスクリプトは：
- 元のシステムでSDPを解く
- データ注入攻撃を実行
- 攻撃後のシステム性能を評価

## 主要な関数

### SDP解法

- `original_thesis.solve_sdp(data, opts)`: 論文のSDP解法
- `regularization_sdp.solve_sdp(data, opts)`: 正則化付きSDP解法

### 勾配計算

- `implicit.dtDelta_dD(...)`: Implicit Differentiationによる勾配計算

### 攻撃実行

- `attack.execute_attack(data, attack_type, ...)`: データ注入攻撃の実行

### ヘルパー関数

- `helper.hinfnorm_AK(A, B, K)`: 閉ループシステム `G_{A,K}(z) = K(zI - A - BK)^{-1}B` のH∞ノルムを計算

## 設定

`+cfg/Const.m`で以下の定数を設定できます：

- `SAMPLE_COUNT`: サンプル数（デフォルト: 20）
- `FD_STEP`: 数値微分のステップサイズ
- `ATTACKER_UPPERLIMIT`: 攻撃側の制約上限
- `EPSILON`: 許容誤差
- `NOISE_BOUND`: ノイズの上限

## エラー処理

- SDPが解けない場合、`error`を出して処理を中断します（適当な値を返しません）
- Dual変数の取得に失敗した場合も同様に`error`を出します

## 理論的背景

詳細な理論は `doc/kkt_implicit_differentiation.tex` を参照してください。このドキュメントには以下が含まれます：

- SDP問題の定式化
- KKT条件
- Implicit Differentiationによる勾配計算
- 各LMIブロック（F1, F2, F3）の定義と微分

## 注意事項

- プロジェクトはMATLABのパッケージ（`+`で始まるディレクトリ）を使用しています
- `startup.m`を実行してパスを設定してください
- SDPソルバー（MOSEKなど）が正しく設定されていることを確認してください

## ライセンス

（必要に応じてライセンス情報を追加）

## 作者

（必要に応じて作者情報を追加）

