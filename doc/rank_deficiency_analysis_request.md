# Rank Deficiency問題の原因探索依頼

## 1. やりたいこと（目的）

半正定値計画問題（SDP）の最適解に対するImplicit Differentiationを実装していますが、**KKT条件の線形化システム `Matrix_H * dx/dD = -Matrix_y` において、`Matrix_H`がrank deficient（ランク不足）になっています**。

この問題の原因を特定し、解決策を提案していただきたいです。

---

## 2. 問題の詳細

### 2.1 問題の概要

- **システム**: `Matrix_H` は `[行数, 列数] = [630, 324]` のサイズ
- **期待されるrank**: 324（列数と等しい）
- **実際のrank**: 321（または323など、列数より1〜3小さい）
- **自由度**: 列数 - rank = 1〜3

つまり、線形システム `Matrix_H * x = -Matrix_y` の解が一意に定まらない状態です。

### 2.2 観察された現象

実際の実行結果（例）：
```
Matrix_H: size=[630,324], rank=321
解は一意ではありません (自由度: 3 = 列数 324 - rank 321)

=== Null Space分析 ===
SVD: 最小特異値 = 8.715856e-10, 最大特異値 = 1.521942e+06, 条件数 = 1.746176e+15
下位10個の特異値:
  1.153033e-03
  8.320529e-04
  3.535804e-04
  2.505978e-04
  1.435300e-04
  4.397260e-05
  1.975690e-05
  9.661836e-09
  1.375024e-09
  8.715856e-10

Null space次元: 3 (tol=1.000000e-06)

各変数ブロックでのnull space成分のノルム:
Null vector 1:
  dL:        6.398317e-05
  dY:        1.000000e+00  ← 主要な成分がdYに集中
  dAlpha:    4.180560e-08
  dBeta:     3.368175e-11
  dtDelta:   8.451825e-12
  dLambda1:  1.935975e-07
  dLambda3:  5.116533e-08
  dLambda_alpha: 4.488990e-16
  dLambda_beta:  3.182878e-09
  dLambda_tDelta: 1.153662e-08
  dLambda_Y: 1.036083e-09

Null vector 2:
  dL:        7.040087e-05
  dY:        1.000000e+00  ← 主要な成分がdYに集中
  ...

Null vector 3:
  dL:        3.050619e-05
  dY:        1.000000e+00  ← 主要な成分がdYに集中
  ...
```

**重要な観察**:
- Null spaceの主要な成分が`dY`（Yの微分）に集中している
- 自由度は1〜3程度（システムサイズによって変動）
- 条件数が非常に大きい（1e15程度）

### 2.3 関連する情報

- **正則化版（regularization_sdp）では問題が解決している**: 同じ問題設定で正則化項を追加したバージョンでは、rankが満たされている
- **SDPの目的関数に正則化項を追加すると解決**: `obj = -tDelta + 1e-8*norm(Y, 'fro')` のような正則化項を追加すると、rankが満たされる

---

## 3. 現状のコードと実装

### 3.1 SDP問題の定式化

**ファイル**: `com_src/calc_sdp/+original_thesis/solve_sdp.m`

```matlab
% 決定変数
Y      = sdpvar(n,n,'symmetric');  % n×n 対称行列
L      = sdpvar(m,n,'full');       % m×n 行列
alpha  = sdpvar(1,1);              % スカラー
beta   = sdpvar(1,1);              % スカラー
tDelta = sdpvar(1,1);              % スカラー（δ²）

% LMIブロック（F1, F2, F3）
[F1, F2, F3] = original_thesis.build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,data);

% 制約
const_1 = [F1 - alpha*F2 >= 0];              % Lambda1に対応
const_2 = [F3 >= tolerance*eye(n+m)];        % Lambda3に対応
const_alpha = [alpha >= 0];                  % Lambda_alphaに対応
const_beta = [beta >= tolerance];            % Lambda_betaに対応
const_tDelta_lower = [tDelta >= tolerance];  % Lambda_tDeltaに対応
const_Y = [Y >= tolerance*eye(n)];          % Lambda_Yに対応
const_tDelta_upper = [(1 - tolerance)^2 >= tDelta];  % Lambda_tDelta_upperに対応

% 目的関数
obj = -tDelta + 1e-6*norm(Y, 'fro') + 1e-6*norm(L, 'fro');
```

**注意**: 目的関数に小さな正則化項（`1e-6*norm(Y, 'fro')`）が含まれていますが、これでもrank不足が発生しています。

### 3.2 LMIブロックの定義

**ファイル**: `com_src/calc_sdp/+original_thesis/build_lmi_blocks.m`

```matlab
function [F1, F2, F3] = build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,sd)
    % F1: (3n+m) × (3n+m) 行列
    row1 = [ Y - tDelta*(B*B') - beta*eye(n),  Znn,  B*L,    Znm ];
    row2 = [ Znn,                          Znn,  Y,    Znm ];
    row3 = [ L'*B',                          Y',  Y,    L' ];
    row4 = [ Zmn,                          Zmn,  L, eye(m) ];
    F1 = [ row1; row2; row3; row4 ];
    
    % F2: データ依存項（alphaを含まない）
    F2 = G*Phi*G';
    
    % F3: (n+m) × (n+m) 行列
    F3 = [Y , L'; L, eye(m) ];
end
```

### 3.3 KKT条件の線形化システム

**ファイル**: `com_src/attck_side/+implicit/dtDelta_dD.m`

KKT条件をデータ`D = [Z', X', U']'`で微分した結果、以下の線形システムが得られます：

```
Matrix_H * [dL; dY; dAlpha; dBeta; dtDelta; dLambda1; dLambda3; dLambda_alpha; dLambda_beta; dLambda_tDelta; dLambda_Y] = -Matrix_y
```

ここで、`Matrix_H`は以下の15個の`G_merge`関数の出力を縦に連結したものです：

```matlab
Matrix_H = [
    G5_sol.G5_row_without_Data;  % ∂L/∂tDelta = 0
    G1_sol.G1_row_without_Data;  % ∂L/∂L = 0
    G2_sol.G2_row_without_Data;  % ∂L/∂Y = 0
    G3_sol.G3_row_without_Data;  % ∂L/∂alpha = 0
    G4_sol.G4_row_without_Data;  % ∂L/∂beta = 0
    G6_sol.G6_row_without_Data;  % 相補性: (F1-α*F2)*Lambda1' = 0
    G7_sol.G7_row_without_Data;  % 相補性: F3*Lambda3' = 0
    G8_sol.G8_row_without_Data;  % 相補性: alpha*Lambda_alpha = 0
    G9_sol.G9_row_without_Data;  % 相補性: (beta-ε)*Lambda_beta = 0
    G10_sol.G10_row_without_Data; % 相補性: (tDelta-ε)*Lambda_tDelta = 0
    G11_sol.G11_row_without_Data; % Lambda1の対称性
    G12_sol.G12_row_without_Data; % Lambda3の対称性
    G13_sol.G13_row_without_Data; % Yの対称性: Y - Y' = 0
    G14_sol.G14_row_without_Data; % 相補性: Y*Lambda_Y' = 0
    G15_sol.G15_row_without_Data; % Lambda_Yの対称性: Lambda_Y - Lambda_Y' = 0
];
```

### 3.4 各G_merge関数の役割

#### G1: LagrangianをLで微分
```matlab
% G1_merge.m
sol.dLambda1 = (-implicit.helper.dF1_dL(n,m, B))';
sol.dLambda3 = (-implicit.helper.dF3_dL(n,m))';
```

#### G2: LagrangianをYで微分
```matlab
% G2_merge.m
sol.dLambda1 = (-implicit.helper.dF1_dY(n,m))';
sol.dLambda3 = (-implicit.helper.dF3_dY(n,m))';
sol.Lambda_Y = speye(n*n);  % Yの制約から
```

#### G5: LagrangianをtDeltaで微分
```matlab
% G5_merge.m
E = [B; zeros(2*n+m,m)];
M = E*E';
sol.dLambda1 = M(:)';  % F1にtDeltaが含まれるため
sol.Lambda_tDelta = -1;  % tDeltaの制約から
```

#### G6: 相補性条件 (F1-α*F2)*Lambda1' = 0
```matlab
% G6_merge.m
% vec((F1-alpha*F2)*Lambda1') = (Lambda1 ⊗ I) vec(F1-alpha*F2)
sub_matrix = kron(Lambda1, eye(3*n+m));
sol.dL = sub_matrix * implicit.helper.dF1_dL(n,m,B);
sol.dY = sub_matrix * implicit.helper.dF1_dY(n,m);
sol.dAlpha = -Lambda1F2(:);  % -vec(Lambda1*F2)
sol.dBeta = -Lambda1M11(:);  % -vec(Lambda1*M11)
sol.dtDelta = -Lambda1M11_withB(:);  % -vec(Lambda1*M11_withB)
sol.dLambda1 = kron(eye, F1_minus_alphaF2) * commutation(...);
```

#### G13: Yの対称性 Y - Y' = 0
```matlab
% G13_merge.m
sol.dY = speye(n*n) - implicit.helper.commutation(n,n);
```

#### G14: 相補性条件 Y*Lambda_Y' = 0
```matlab
% G14_merge.m
sol.dY = kron(Lambda_Y, eye(n));
sol.Lambda_Y = kron(speye(n), Y) * implicit.helper.commutation(n,n);
```

#### G15: Lambda_Yの対称性 Lambda_Y - Lambda_Y' = 0
```matlab
% G15_merge.m
sol.Lambda_Y = speye(n*n) - implicit.helper.commutation(n,n);
```

### 3.5 変数の次元

- `dL`: `n*m`次元
- `dY`: `n*n`次元（対称性により`n*(n+1)/2`に削減される可能性があるが、実装では`n*n`）
- `dAlpha`: 1次元
- `dBeta`: 1次元
- `dtDelta`: 1次元
- `dLambda1`: `(3*n+m)^2`次元
- `dLambda3`: `(n+m)^2`次元
- `dLambda_alpha`: 1次元
- `dLambda_beta`: 1次元
- `dLambda_tDelta`: 1次元
- `dLambda_Y`: `n*n`次元

**合計列数**: `n*m + n*n + 1 + 1 + 1 + (3*n+m)^2 + (n+m)^2 + 1 + 1 + 1 + n*n = 324` (n=4, m=3の場合)

### 3.6 重要な実装の詳細

1. **相補性条件の形式**: `(制約) * Lambda^T = 0` の形式を使用
   - 例: `(F1 - alpha*F2) * Lambda1' = 0`
   - 例: `Y * Lambda_Y' = 0`

2. **対称性の扱い**: 
   - `Y`は`sdpvar(n,n,'symmetric')`として定義されているが、微分では`n*n`次元のベクトルとして扱う
   - `G13`で`Y - Y' = 0`の制約を追加
   - `G15`で`Lambda_Y - Lambda_Y' = 0`の制約を追加

3. **Commutation行列**: `K_{p,q}`は`vec(A^T) = K_{p,q} * vec(A)`を満たす行列

---

## 4. 推測される原因の候補

### 4.1 Yに関する制約の重複

Null space分析で`dY`に主要な成分が集中していることから、以下の可能性が考えられます：

1. **G2とG13の相互作用**: 
   - `G2`: `∂L/∂Y = 0`（Yに関する最適性条件）
   - `G13`: `Y - Y' = 0`（Yの対称性）
   - これらが線形従属になっている可能性

2. **G13とG14の相互作用**:
   - `G13`: `Y - Y' = 0`
   - `G14`: `Y * Lambda_Y' = 0`（相補性条件）
   - `Y`が対称であることを考慮すると、`G14`の制約が冗長になっている可能性

3. **G2とG14の相互作用**:
   - `G2`: `∂L/∂Y = 0`には`Lambda_Y`の項が含まれる
   - `G14`: `Y * Lambda_Y' = 0`
   - これらが線形従属になっている可能性

### 4.2 SDPの最適解の非一意性

- 目的関数に`1e-6*norm(Y, 'fro')`という小さな正則化項があるが、それでも`Y`の最適解が一意に定まらない可能性
- より大きな正則化項（`1e-8*norm(Y, 'fro')`など）を追加するとrankが満たされることから、SDPの最適解自体が非一意である可能性が高い

### 4.3 対称性制約の扱い

- `Y`は`sdpvar(n,n,'symmetric')`として定義されているが、微分では`n*n`次元のベクトルとして扱っている
- `G13`で`Y - Y' = 0`の制約を追加しているが、これが他の制約と線形従属になっている可能性

---

## 5. 質問事項

以下の点について、ご意見をお聞かせください：

1. **Yに関する制約の分析**:
   - `G2`（`∂L/∂Y = 0`）、`G13`（`Y - Y' = 0`）、`G14`（`Y * Lambda_Y' = 0`）の間で線形従属が発生している可能性はありますか？
   - 特に、`Y`が対称であることを考慮すると、これらの制約のうちどれかが冗長になっていませんか？

2. **SDPの最適解の一意性**:
   - 目的関数に小さな正則化項があるにもかかわらず、`Y`の最適解が一意に定まらない原因は何でしょうか？
   - より大きな正則化項を追加するとrankが満たされることから、SDPの最適解自体が非一意である可能性が高いですが、これは理論的に正しい解釈でしょうか？

3. **対称性制約の実装**:
   - `Y`を`sdpvar(n,n,'symmetric')`として定義し、さらに`G13`で`Y - Y' = 0`の制約を追加するのは冗長ではありませんか？
   - 微分計算において、対称行列を`n*n`次元のベクトルとして扱う場合、対称性をどのように考慮すべきでしょうか？

4. **Lambda_Yに関する制約**:
   - `G14`（`Y * Lambda_Y' = 0`）と`G15`（`Lambda_Y - Lambda_Y' = 0`）の間で線形従属が発生している可能性はありますか？
   - `Y`が対称で、`Lambda_Y`も対称である場合、これらの制約はどのように相互作用しますか？

5. **解決策の提案**:
   - このrank不足問題を解決するには、どのような方法が考えられますか？
   - 正則化項を追加する以外に、KKT条件の線形化システム自体を修正する方法はありますか？
   - 理論的には、どの制約が冗長であるべきで、それをどのように除去または修正すべきでしょうか？

6. **数値的な問題**:
   - 条件数が非常に大きい（1e15程度）ことから、数値的な不安定性も懸念されます。これはrank不足の原因ではなく、結果として現れている現象でしょうか？

---

## 6. 参考情報

### 6.1 関連ファイル

- SDP解法: `com_src/calc_sdp/+original_thesis/solve_sdp.m`
- LMIブロック構築: `com_src/calc_sdp/+original_thesis/build_lmi_blocks.m`
- Implicit Differentiation: `com_src/attck_side/+implicit/dtDelta_dD.m`
- G_merge関数群: `com_src/attck_side/+implicit/+G_grad/`
- 理論的説明: `doc/kkt_implicit_differentiation.tex`

### 6.2 テストスクリプト

- `scripts/demo_implicit_grad.m`: Implicit Differentiationのテストスクリプト

### 6.3 正則化版との比較

正則化版（`regularization_sdp`）では、異なるSDP定式化を使用しており、rank不足が発生していません。比較することで、問題の原因を特定できる可能性があります。

---

## 7. 現在の実装の詳細（各G_merge関数の計算内容）

以下、各G_merge関数がどのような計算を行っているかを詳しく説明します。

### 7.1 G1: LagrangianをLで微分（∂L/∂L = 0）

**ファイル**: `G1_merge.m`

```matlab
sol.dLambda1 = (-implicit.helper.dF1_dL(n,m, B))';
sol.dLambda3 = (-implicit.helper.dF3_dL(n,m))';
```

**計算内容**:
- LagrangianをLで微分すると、`∂L/∂L = -dF1_dL' * vec(Lambda1) - dF3_dL' * vec(Lambda3) = 0`
- ここで、`dF1_dL`は`vec(F1)`を`vec(L)`で微分した行列（サイズ: `(3n+m)^2 × nm`）
- `dF3_dL`は`vec(F3)`を`vec(L)`で微分した行列（サイズ: `(n+m)^2 × nm`）

**dF1_dLの計算** (`dF1_dL.m`):
```matlab
% F1の各ブロックがLに依存する部分を抽出
E1 = [speye(n); zeros(2*n+m,n)];  % F1の(1,1)ブロック用
E2 = [zeros(2*n,n); speye(n); zeros(m,n)];  % F1の(2,3)ブロック用
E3 = [zeros(3*n,m); speye(m)];  % F1の(4,3)ブロック用

% Kronecker積とcommutation行列を使用して微分を計算
term1 = kron(E1, E2*B);  % B*Lの項
term2 = kron(E2*B, E1) * commutation(m,n);  % L'*B'の項（対称性のため）
term3 = kron(E3, E2);  % Lの項
term4 = kron(E2, E3) * commutation(n,m);  % L'の項（対称性のため）

dF1_dL = term1 + term2 + term3 + term4;
```

### 7.2 G2: LagrangianをYで微分（∂L/∂Y = 0）

**ファイル**: `G2_merge.m`

```matlab
sol.dLambda1 = (-implicit.helper.dF1_dY(n,m))';
sol.dLambda3 = (-implicit.helper.dF3_dY(n,m))';
sol.Lambda_Y = speye(n*n);  % Y >= tolerance*eye(n)の制約から
```

**計算内容**:
- LagrangianをYで微分すると、`∂L/∂Y = -dF1_dY' * vec(Lambda1) - dF3_dY' * vec(Lambda3) + vec(Lambda_Y) = 0`
- `dF1_dY`は`vec(F1)`を`vec(Y)`で微分した行列（サイズ: `(3n+m)^2 × n^2`）
- `dF3_dY`は`vec(F3)`を`vec(Y)`で微分した行列（サイズ: `(n+m)^2 × n^2`）
- `Lambda_Y`は`Y >= tolerance*eye(n)`の制約に対するラグランジュ乗数

**dF1_dYの計算** (`dF1_dY.m`):
```matlab
% F1の各ブロックがYに依存する部分を抽出
E1 = [speye(n); zeros(2*n+m,n)];  % F1の(1,1)ブロック: Y - tDelta*(B*B') - beta*I
E2 = [zeros(n,n); speye(n); zeros(n+m,n)];  % F1の(2,3)ブロック: Y
E3 = [zeros(n,n); zeros(n,n); speye(n); zeros(m,n)];  % F1の(3,2)ブロック: Y'

% Kronecker積とcommutation行列を使用
dF1_dY = kron(E1, E1) + kron(E2, E2) + kron(E3, E2) + kron(E2, E3) * commutation(n,n);
```

### 7.3 G3: Lagrangianをalphaで微分（∂L/∂alpha = 0）

**ファイル**: `G3_merge.m`

```matlab
sol.dLambda1 = (F2(:))';  % F1 - alpha*F2の項から
sol.Lambda_Alpha = -1;  % alpha >= 0の制約から
sol.Data = Lambda1(:)' * implicit.helper.dF2_dD(n,m,B,T,G,Phi);  % データ依存項
```

**計算内容**:
- Lagrangianをalphaで微分すると、`∂L/∂alpha = vec(F2)' * vec(Lambda1) - Lambda_alpha = 0`
- `F2 = G*Phi*G'`はデータ依存項なので、`dF2/dD`の項が`Data`に含まれる

### 7.4 G4: Lagrangianをbetaで微分（∂L/∂beta = 0）

**ファイル**: `G4_merge.m`

```matlab
E = [speye(n); zeros(2*n+m,n)];
M = E*E';  % M = [I_n; 0; 0; 0] * [I_n; 0; 0; 0]'
sol.dLambda1 = M(:)';  % F1のbeta*I_n項から
sol.Lambda_Beta = -1;  % beta >= toleranceの制約から
```

**計算内容**:
- Lagrangianをbetaで微分すると、`∂L/∂beta = -vec(M)' * vec(Lambda1) - Lambda_beta = 0`
- ここで、`M`は`F1`の`beta*I_n`項に対応する行列

### 7.5 G5: LagrangianをtDeltaで微分（∂L/∂tDelta = 0）

**ファイル**: `G5_merge.m`

```matlab
E = [B; zeros(2*n+m,m)];
M = E*E';  % M = [B; 0; 0; 0] * [B; 0; 0; 0]'
sol.dLambda1 = M(:)';  % F1のtDelta*(B*B')項から
sol.Lambda_tDelta = -1;  % tDelta >= toleranceの制約から
```

**計算内容**:
- LagrangianをtDeltaで微分すると、`∂L/∂tDelta = -vec(M)' * vec(Lambda1) - Lambda_tDelta - 1 = 0`
- 最後の`-1`は目的関数`-tDelta`から来る

### 7.6 G6: 相補性条件 (F1-α*F2)*Lambda1' = 0 をデータで微分

**ファイル**: `G6_merge.m`

**計算内容**:
相補性条件`(F1 - alpha*F2) * Lambda1' = 0`をデータ`D`で微分すると：

```matlab
% vec((F1-alpha*F2)*Lambda1') = (Lambda1 ⊗ I) * vec(F1-alpha*F2)
sub_matrix = kron(Lambda1, speye(3*n+m));

% d/dL [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) * dvec(F1)/dL
sol.dL = sub_matrix * implicit.helper.dF1_dL(n,m,B);

% d/dY [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) * dvec(F1)/dY
sol.dY = sub_matrix * implicit.helper.dF1_dY(n,m);

% d/dalpha [vec((F1-alpha*F2)*Lambda1')] = -(Lambda1 ⊗ I) * vec(F2) = -vec(Lambda1*F2)
Lambda1F2 = Lambda1 * F2;
sol.dAlpha = -Lambda1F2(:);

% d/dbeta [vec((F1-alpha*F2)*Lambda1')] = -(Lambda1 ⊗ I) * vec(M11) = -vec(Lambda1*M11)
E11 = [speye(n); zeros(2*n+m,n)];
M11 = E11*E11';
Lambda1M11 = Lambda1 * M11;
sol.dBeta = -Lambda1M11(:);

% d/dtDelta [vec((F1-alpha*F2)*Lambda1')] = -(Lambda1 ⊗ I) * vec(M11_withB) = -vec(Lambda1*M11_withB)
E11_withB = [B; zeros(2*n+m,m)];
M11_withB = E11_withB*E11_withB';
Lambda1M11_withB = Lambda1 * M11_withB;
sol.dtDelta = -Lambda1M11_withB(:);

% d/dLambda1 [vec((F1-alpha*F2)*Lambda1')] = (I ⊗ (F1-alpha*F2)) * K
% ここで、Kはcommutation行列（vec(Lambda1') = K * vec(Lambda1)）
sol.dLambda1 = kron(speye(3*n+m), F1_minus_alphaF2) * commutation(3*n+m, 3*n+m);

% データ依存項: d/dD [vec((F1-alpha*F2)*Lambda1')] = -alpha * (Lambda1 ⊗ I) * dvec(F2)/dD
sol.Data = -alpha * sub_matrix * implicit.helper.dF2_dD(n,m,B,T,G,Phi);
```

### 7.7 G7: 相補性条件 F3*Lambda3' = 0 をデータで微分

**ファイル**: `G7_merge.m`

**計算内容**:
相補性条件`F3 * Lambda3' = 0`をデータ`D`で微分すると：

```matlab
% vec(F3*Lambda3') = (Lambda3 ⊗ I) * vec(F3)
sub_matrix = kron(Lambda3, speye(n+m));

% d/dL [vec(F3*Lambda3')] = (Lambda3 ⊗ I) * dvec(F3)/dL
sol.dL = sub_matrix * implicit.helper.dF3_dL(n,m);

% d/dY [vec(F3*Lambda3')] = (Lambda3 ⊗ I) * dvec(F3)/dY
sol.dY = sub_matrix * implicit.helper.dF3_dY(n,m);

% d/dLambda3 [vec(F3*Lambda3')] = (I ⊗ F3) * K
sol.dLambda3 = kron(speye(n+m), F3) * commutation(n+m, n+m);
```

### 7.8 G8: 相補性条件 alpha*Lambda_alpha = 0 をデータで微分

**ファイル**: `G8_merge.m`

```matlab
sol.dAlpha = Lambda_alpha;  % d/dalpha [alpha*Lambda_alpha] = Lambda_alpha
sol.Lambda_Alpha = alpha;  % d/dLambda_alpha [alpha*Lambda_alpha] = alpha
```

**計算内容**:
- 相補性条件`alpha * Lambda_alpha = 0`をデータで微分すると、`Lambda_alpha * dalpha + alpha * dLambda_alpha = 0`

### 7.9 G9: 相補性条件 (beta-ε)*Lambda_beta = 0 をデータで微分

**ファイル**: `G9_merge.m`

```matlab
sol.dBeta = Lambda_beta;  % d/dbeta [(beta-ε)*Lambda_beta] = Lambda_beta
sol.Lambda_Beta = beta;  % d/dLambda_beta [(beta-ε)*Lambda_beta] = beta-ε ≈ beta
```

### 7.10 G10: 相補性条件 (tDelta-ε)*Lambda_tDelta = 0 をデータで微分

**ファイル**: `G10_merge.m`

```matlab
sol.dtDelta = Lambda_tDelta;  % d/dtDelta [(tDelta-ε)*Lambda_tDelta] = Lambda_tDelta
sol.Lambda_tDelta = tDelta;  % d/dLambda_tDelta [(tDelta-ε)*Lambda_tDelta] = tDelta-ε ≈ tDelta
```

### 7.11 G11: Lambda1の対称性 Lambda1 - Lambda1' = 0

**ファイル**: `G11_merge.m`

```matlab
sol.dLambda1 = speye((3*n+m)^2) - commutation(3*n+m, 3*n+m);
```

**計算内容**:
- `vec(Lambda1 - Lambda1') = vec(Lambda1) - vec(Lambda1') = (I - K) * vec(Lambda1) = 0`
- ここで、`K`はcommutation行列（`vec(Lambda1') = K * vec(Lambda1)`）

### 7.12 G12: Lambda3の対称性 Lambda3 - Lambda3' = 0

**ファイル**: `G12_merge.m`

```matlab
sol.dLambda3 = speye((n+m)^2) - commutation(n+m, n+m);
```

**計算内容**: G11と同様

### 7.13 G13: Yの対称性 Y - Y' = 0

**ファイル**: `G13_merge.m`

```matlab
sol.dY = speye(n*n) - commutation(n, n);
```

**計算内容**:
- `vec(Y - Y') = vec(Y) - vec(Y') = (I - K) * vec(Y) = 0`
- ただし、`Y`は`sdpvar(n,n,'symmetric')`として定義されているため、この制約は理論的には冗長

### 7.14 G14: 相補性条件 Y*Lambda_Y' = 0 をデータで微分

**ファイル**: `G14_merge.m`

**計算内容**:
相補性条件`Y * Lambda_Y' = 0`をデータ`D`で微分すると：

```matlab
% vec(Y*Lambda_Y') = (Lambda_Y ⊗ I) * vec(Y)
sol.dY = kron(Lambda_Y, speye(n));  % d/dY [vec(Y*Lambda_Y')] = Lambda_Y ⊗ I

% d/dLambda_Y [vec(Y*Lambda_Y')] = (I ⊗ Y) * K
sol.Lambda_Y = kron(speye(n), Y) * commutation(n, n);
```

### 7.15 G15: Lambda_Yの対称性 Lambda_Y - Lambda_Y' = 0

**ファイル**: `G15_merge.m`

```matlab
sol.Lambda_Y = speye(n*n) - commutation(n, n);
```

**計算内容**: G11, G12, G13と同様

### 7.16 Commutation行列の役割

**ファイル**: `commutation.m`

Commutation行列`K_{p,q}`は、`vec(A^T) = K_{p,q} * vec(A)`を満たす行列です。

```matlab
% Aがp×q行列の場合、K_{p,q}は(pq)×(pq)行列
% vec(A)のインデックス: i + (j-1)*p (i行j列)
% vec(A^T)のインデックス: j + (i-1)*q (j行i列)
```

この行列は、相補性条件の微分において`vec(Lambda^T)`を`vec(Lambda)`で表現するために使用されます。

### 7.17 Matrix_Hの構築

最終的に、`Matrix_H`は以下のように構築されます：

```matlab
Matrix_H = [
    G5_sol.G5_row_without_Data;  % 1行
    G1_sol.G1_row_without_Data;  % n*m行
    G2_sol.G2_row_without_Data;  % n*n行
    G3_sol.G3_row_without_Data;  % 1行
    G4_sol.G4_row_without_Data;  % 1行
    G6_sol.G6_row_without_Data;  % (3*n+m)^2行
    G7_sol.G7_row_without_Data;  % (n+m)^2行
    G8_sol.G8_row_without_Data;  % 1行
    G9_sol.G9_row_without_Data;  % 1行
    G10_sol.G10_row_without_Data; % 1行
    G11_sol.G11_row_without_Data; % (3*n+m)^2行
    G12_sol.G12_row_without_Data; % (n+m)^2行
    G13_sol.G13_row_without_Data; % n*n行
    G14_sol.G14_row_without_Data; % n*n行
    G15_sol.G15_row_without_Data; % n*n行
];
```

**合計行数**: `1 + n*m + n*n + 1 + 1 + (3*n+m)^2 + (n+m)^2 + 1 + 1 + 1 + (3*n+m)^2 + (n+m)^2 + n*n + n*n + n*n = 630` (n=4, m=3の場合)

**列数**: `n*m + n*n + 1 + 1 + 1 + (3*n+m)^2 + (n+m)^2 + 1 + 1 + 1 + n*n = 324` (n=4, m=3の場合)

---

## 8. お願い

この問題の原因を特定し、理論的根拠に基づいた解決策を提案していただけますでしょうか。特に、以下の点について詳しくご説明いただけると助かります：

1. **理論的な分析**: なぜrank不足が発生するのか、どの制約が冗長なのか
2. **実装上の問題**: 現在の実装に問題があるのか、それとも理論的な問題なのか
3. **解決策**: 正則化項を追加する以外の方法があるのか、あるいは正則化項の追加が唯一の解決策なのか

ご回答をお待ちしております。

