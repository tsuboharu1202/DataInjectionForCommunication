# 2つのSDP実装の相違点

## 1. 提供されたコード（Xingchen Liの論文実装）

## 2. original_thesis/solve_sdp.m

---

## 主な相違点

### 1. 変数名の違い

| Xingchen Liのコード | original_thesis/solve_sdp.m | 説明 |
|-------------------|---------------------------|------|
| `X` | `L` | 制御ゲイン行列（m×n） |
| `a` | `alpha` | スカラー変数 |
| `b` | `beta` | スカラー変数 |
| `d` | `tDelta` | δ²（delta squared） |

**注意**: `Y`は両方で同じ（n×n対称行列）

### 2. データの扱い方

#### Xingchen Liのコード:
```matlab
% データ生成関数 gen_data を使用
[Xm, Xp, U, W] = gen_data(A, B, T, w_max);
% Xm: X_minus (過去の状態)
% Xp: X_plus (次の状態)
```

#### original_thesis/solve_sdp.m:
```matlab
% SystemDataオブジェクトから取得
Z = data.Z;     % n×T (Xpに相当)
Xm = data.X;    % n×T (Xmに相当)
Um = data.U;    % m×T (Uに相当)
```

**相違点**: 
- Xingchen Liのコード: `Xp`と`Xm`を明示的に区別
- original_thesis: `Z`（Xpに相当）と`X`（Xmに相当）を使用

### 3. G行列（データ依存項）の定義

#### Xingchen Liのコード:
```matlab
N = [eye(n), Xp-B*U;
     zeros(n,n), -Xm;
     zeros(n+m,n+T)] * ...
     [phi11 phi12; phi12' phi22] * ...
     [eye(n), Xp-B*U;
      zeros(n,n), -Xm;
      zeros(n+m,n+T)]';
```
**サイズ**: `(n+n+n+m) × (n+T)` = `(3n+m) × (n+T)`

#### original_thesis/solve_sdp.m:
```matlab
G = [ eye(n)  ,  Z-B*Um ;      % 1: n
      zeros(n,n) , -Xm ;        % 2: n
      zeros(n,n) , zeros(n,T) ; % 3: n
      zeros(m,n+T)];            % 4: m
```
**サイズ**: `(n+n+n+m) × (n+T)` = `(3n+m) × (n+T)`

**相違点**: 
- 基本的に同じ構造だが、Xingchen Liのコードでは`N = G*Phi*G'`を直接計算
- original_thesisでは`G`と`Phi`を分離し、`F2 = G*Phi*G'`として後で計算

### 4. M行列（F1に相当）の定義

#### Xingchen Liのコード:
```matlab
M = [Y-d*(B*B')-b*eye(n), zeros(n), B*X, zeros(n,m);
     zeros(n), zeros(n), Y, zeros(n,m);
     X'*B', Y', Y, X';
     zeros(m,n), zeros(m,n), X, eye(m)];
```

#### original_thesis/build_lmi_blocks.m:
```matlab
row1 = [ Y - tDelta*(B*B') - beta*eye(n),  Znn,  B*L,    Znm ];
row2 = [ Znn,                          Znn,  Y,    Znm ];
row3 = [ L'*B',                          Y',  Y,    L' ];
row4 = [ Zmn,                          Zmn,  L, eye(m) ];
F1 = [ row1; row2; row3; row4 ];
```

**相違点**: 
- **実質的に同じ**。変数名のみ異なる（`X` vs `L`, `d` vs `tDelta`, `b` vs `beta`）

### 5. 制約条件

#### Xingchen Liのコード:
```matlab
constraints = [Y>=0, a>=0, b>=0, [Y X'; X eye(m)]>=0, M-a*N>=0];
```

制約の内訳:
1. `Y >= 0` (Yは半正定値)
2. `a >= 0` (alpha >= 0)
3. `b >= 0` (beta >= 0)
4. `[Y X'; X eye(m)] >= 0` (F3に相当)
5. `M - a*N >= 0` (F1 - alpha*F2 >= 0)

#### original_thesis/solve_sdp.m:
```matlab
const_1 = [F1 - alpha*F2 >= 0];
const_2 = [F3 >= tolerance*eye(n+m)];
const_alpha = [alpha >= 0];
const_beta = [beta >= tolerance];
const_tDelta_lower = [tDelta >= tolerance];
const_Y = [Y >= tolerance*eye(n)];
const_tDelta_upper = [(1 - tolerance)^2 >= tDelta];
```

**主な相違点**:

1. **betaの制約**:
   - Xingchen Li: `b >= 0` (非負のみ)
   - original_thesis: `beta >= tolerance` (正の下界、tolerance=1e-8)

2. **Yの制約**:
   - Xingchen Li: `Y >= 0` (半正定値のみ)
   - original_thesis: `Y >= tolerance*eye(n)` (正定値、数値的安定性のため)

3. **tDelta (d)の制約**:
   - Xingchen Li: 明示的な制約なし（暗黙的に`d >= 0`）
   - original_thesis: 
     - `tDelta >= tolerance` (下界)
     - `tDelta <= (1 - tolerance)^2` (上界)

4. **F3の制約**:
   - Xingchen Li: `[Y X'; X eye(m)] >= 0` (半正定値のみ)
   - original_thesis: `F3 >= tolerance*eye(n+m)` (正定値、数値的安定性のため)

### 6. 目的関数

#### Xingchen Liのコード:
```matlab
solution = optimize(constraints, -d, options);
```
**目的関数**: `-d` (delta²を最大化)

#### original_thesis/solve_sdp.m:
```matlab
obj = -tDelta;
diagnostics = optimize(constr, obj, params);
```
**目的関数**: `-tDelta` (delta²を最大化)

**相違点**: **同じ**（変数名のみ異なる）

### 7. 出力と後処理

#### Xingchen Liのコード:
```matlab
K = value(X)/value(Y);
fprintf('The value of max(abs(eig(A+B*K))) is: %f\n', max(abs(eig(A+B*K))));
fprintf('The value of d is: %f\n', value(d));
```

#### original_thesis/solve_sdp.m:
```matlab
K = sol.L/sol.Y;
delta = sqrt(value(tDelta));  % tDelta = δ²なので、δ = sqrt(tDelta)
sol.delta = delta;
sol.rho = (1-delta)/(1+delta);  % Theorem 3: ρ = (1-δ)/(1+δ)
```

**相違点**:
- Xingchen Li: `d`（delta²）を直接出力
- original_thesis: `delta = sqrt(tDelta)`としてdeltaを計算し、さらに`rho = (1-delta)/(1+delta)`を計算

### 8. Dual変数（ラグランジュ乗数）の取得

#### Xingchen Liのコード:
```matlab
% Dual変数の取得なし
```

#### original_thesis/solve_sdp.m:
```matlab
sol.Lambda1 = dual(const_1);
sol.Lambda3 = dual(const_2);
sol.Lambda_alpha = dual(const_alpha);
sol.Lambda_beta = dual(const_beta);
sol.Lambda_tDelta = dual(const_tDelta_lower);
sol.Lambda_Y = dual(const_Y);
sol.Lambda_tDelta_upper = dual(const_tDelta_upper);
```

**相違点**: 
- Xingchen LiのコードはDual変数を取得しない
- original_thesisはImplicit DifferentiationのためにDual変数を取得

### 9. エラーハンドリング

#### Xingchen Liのコード:
```matlab
if solution.problem == 0
    fprintf('The solution is feasible\n');
else
    fprintf('The solution is infeasible\n');
end
```

#### original_thesis/solve_sdp.m:
```matlab
if diagnostics.problem ~= 0
    error('solve_sdp:OptimizationFailed', ...
        '最適化が解けませんでした。status: %d, info: %s', ...
        diagnostics.problem, diagnostics.info);
end
```

**相違点**:
- Xingchen Li: 警告を出すのみ
- original_thesis: エラーをthrowして処理を中断

### 10. 正則化項

#### Xingchen Liのコード:
```matlab
% 目的関数に正則化項なし
obj = -d;
```

#### original_thesis/solve_sdp.m:
```matlab
% コメントアウトされているが、以前は正則化項があった
% obj = -tDelta + 1e-6*norm(Y, 'fro') + 1e-6*norm(L, 'fro');
obj = -tDelta;
```

**相違点**: 
- 現在の実装では両方とも正則化項なし
- original_thesisには過去に正則化項があった痕跡（コメントアウト）

---

## まとめ

### 本質的に同じ部分
1. **SDP問題の構造**: 両方とも同じSDP問題を解いている
2. **M行列（F1）**: 実質的に同じ
3. **N行列（F2）**: `G*Phi*G'`として同じ
4. **目的関数**: delta²を最大化

### 主な相違点
1. **数値的安定性**: original_thesisは`tolerance`を使用して正定値制約を追加
2. **制約の厳しさ**: original_thesisはより厳しい制約（下界・上界）
3. **Dual変数の取得**: original_thesisはImplicit DifferentiationのためにDual変数を取得
4. **出力形式**: original_thesisはより詳細な出力（delta, rhoなど）
5. **エラーハンドリング**: original_thesisはより厳格

### 実装上の違い
- **変数名**: `X` vs `L`, `a` vs `alpha`, `b` vs `beta`, `d` vs `tDelta`
- **データ構造**: Xingchen Liは関数引数、original_thesisはSystemDataオブジェクト
- **コード構造**: original_thesisはよりモジュール化（build_lmi_blocks関数を使用）

---

## 結論

両方の実装は**同じSDP問題を解いている**が、original_thesisの実装は以下の点で改善されている：

1. **数値的安定性**: toleranceを使用した正定値制約
2. **拡張性**: Dual変数の取得によりImplicit Differentiationに対応
3. **保守性**: モジュール化されたコード構造
4. **エラーハンドリング**: より厳格なエラー処理

ただし、Xingchen Liのコードは**論文の実装そのまま**であり、理論的な正しさは保証されている。



