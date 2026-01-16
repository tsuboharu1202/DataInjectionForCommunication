# Original Thesis版SDPのImplicit Differentiation実装の検証依頼

## 1. やりたいこと（目的）

論文「Data-driven Quantized Control of Partially Unknown Linear Systems with Noises」に基づくSDP問題に対して、**KKT条件を用いたImplicit Differentiation**により、**delta（量子化パラメータ）のデータに対する勾配**を計算する実装を行っています。

具体的には、データ`D = [Z', X', U']'`（Z: 次状態、X: 状態、U: 入力）に対する`delta`の勾配`d(delta)/dD`を計算したいです。

この実装の**計算やコーディングミスがないか**を確認していただきたいです。

---

## 2. SDP問題の定式化

### 2.1 決定変数

- `Y ∈ ℝ^(n×n)`: 対称正定値行列
- `L ∈ ℝ^(m×n)`: 制御ゲイン行列
- `alpha ∈ ℝ`: スカラー（alpha >= 0）
- `beta ∈ ℝ`: スカラー（beta > 0）
- `tDelta ∈ ℝ`: スカラー（tDelta > 0）、実際には`delta²`（delta > 0）

### 2.2 制約条件

1. **LMI制約**: `F1 - alpha*F2 >= 0`
2. **LMI制約**: `F3 >= tolerance*eye(n+m)`
3. **非負制約**: `alpha >= 0`
4. **正定値制約**: `beta >= tolerance`
5. **正定値制約**: `tDelta >= tolerance`
6. **正定値制約**: `Y >= tolerance*eye(n)`
7. **上界制約**: `tDelta <= (1 - tolerance)^2`

### 2.3 目的関数

```
minimize: -tDelta
```

**注意**: 目的関数に正則化項は含まれていない（コメントアウトされている）

### 2.4 LMIブロックの定義

#### F1: (3n+m) × (3n+m) 行列

```matlab
row1 = [ Y - tDelta*(B*B') - beta*eye(n),  Znn,  B*L,    Znm ];
row2 = [ Znn,                          Znn,  Y,    Znm ];
row3 = [ L'*B',                          Y',  Y,    L' ];
row4 = [ Zmn,                          Zmn,  L, eye(m) ];
F1 = [ row1; row2; row3; row4 ];
```

#### F2: データ依存項（alphaを含まない）

```matlab
F2 = G*Phi*G';
```

ここで：
```matlab
G = [ eye(n)  ,  Z-B*Um ;
     zeros(n,n) , -Xm ;
     zeros(n,n) , zeros(n,T) ;
     zeros(m,n+T)];
Phi = [Phi11 Phi12; Phi12' Phi22];
```

#### F3: (n+m) × (n+m) 行列

```matlab
F3 = [Y , L';
      L, eye(m) ];
```

---

## 3. 実装コード

### 3.1 SDP解法: `original_thesis/solve_sdp.m`

```matlab
function [sol, K, Y, L, diagnostics] = solve_sdp(data, opts)

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'verbose'),   opts.verbose   = 0;     end
if ~isfield(opts,'solver'),    opts.solver    = 'mosek'; end

% -------------------------
% Unpack data
% -------------------------
B = data.B;
Z     = data.Z;     % n×T
Xm     = data.X;    % n×T
Um     = data.U;    % m×T  (m=1 でも可)

Gamma_Matrix = [Um;Xm];

Phi11  = data.Phi11;
Phi12  = data.Phi12;
Phi22  = data.Phi22;

[n, T] = size(Z);
m = size(Um,1);

% -------------------------
% Decision variables
% -------------------------
Y      = sdpvar(n,n,'symmetric');
L      = sdpvar(m,n,'full');
alpha  = sdpvar(1,1);
beta   = sdpvar(1,1);
tDelta = sdpvar(1,1);     % = δ²

% -------------------------
G = [ eye(n)  ,  Z-B*Um ;
     zeros(n,n) , -Xm ;
     zeros(n,n) , zeros(n,T) ;
     zeros(m,n+T)];

Phi  = [Phi11 Phi12; Phi12' Phi22];

% -------------------------
% LMI blocks (25)(26)(27)
% -------------------------
[F1, F2, F3] = original_thesis.build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,data);

% -------------------------
% Constraints
% -------------------------

tolerance = 1e-9;

% 制約を個別に定義（dual取得のため）
const_1 = [F1 - alpha*F2 >= 0];  % Lambda1に対応 (F2 = G*Phi*G'はalphaを含まない)
const_2 = [F3 >= tolerance*eye(n+m)];  % Lambda3に対応
const_alpha = [alpha >= 0];  % Lambda_alphaに対応
const_beta = [beta >= tolerance];  % Lambda_betaに対応
const_tDelta_lower = [tDelta >= tolerance];  % Lambda_tDelta_lowerに対応
const_Y = [Y >= tolerance*eye(n)];  % Yの制約（必要に応じて）
const_tDelta_upper = [(1 - tolerance)^2 >= tDelta];  % tDeltaの上界制約

constr  = [];
constr  = [constr, const_1];
constr  = [constr, const_2];
constr  = [constr, const_alpha];
constr  = [constr, const_beta];
constr  = [constr, const_tDelta_lower];
constr  = [constr, const_Y];
constr  = [constr, const_tDelta_upper];

% Objective: maximize δ  <=>  minimize tDelta
obj = -tDelta;

params = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
diagnostics = optimize(constr, obj, params);

% -------------------------
% 最適化結果のチェック
% -------------------------
sol.status = diagnostics.problem;

if diagnostics.problem ~= 0
    error('solve_sdp:OptimizationFailed', ...
        '最適化が解けませんでした。status: %d, info: %s', ...
        diagnostics.problem, diagnostics.info);
end

% -------------------------
% Output pack
% -------------------------
sol.Y         = value(Y);
sol.L         = value(L);
sol.alpha     = value(alpha);
sol.beta      = value(beta);
% tDelta = δ^{2}なので、δ = (tDelta)^(1/2)
delta = sqrt(value(tDelta));
sol.delta     = delta;
sol.rho = (1-delta)/(1+delta);  % Theorem 3: ρ = (1-δ)/(1+δ)
sol.objective = value(obj);

K = sol.L/sol.Y;
sol.K = K;

% -------------------------
% Dual values (Lagrange multipliers)
% -------------------------
try
    sol.Lambda1 = dual(const_1);  % F1 - α*F2 >= 0 のdual（サイズ: (3*n+m) × (3*n+m)）
    sol.Lambda3 = dual(const_2);  % F3 >= tolerance*eye(n+m) のdual（サイズ: (n+m) × (n+m)）
    sol.Lambda_alpha = dual(const_alpha);  % alpha >= 0 のdual（スカラー）
    sol.Lambda_beta = dual(const_beta);  % beta >= tolerance のdual（スカラー）
    sol.Lambda_tDelta = dual(const_tDelta_lower);  % tDelta >= tolerance のdual（スカラー）
    sol.Lambda_Y = dual(const_Y);  % Y >= tolerance*eye(n) のdual（サイズ: n × n）
    sol.Lambda_tDelta_upper = dual(const_tDelta_upper);  % (1 - tolerance)^2 >= tDelta のdual（スカラー）
catch ME
    error('solve_sdp:DualFailed', ...
        'Dual変数の取得に失敗しました: %s', ME.message);
end
end
```

### 3.2 LMIブロック構築: `original_thesis/build_lmi_blocks.m`

```matlab
function [F1, F2, F3] = build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,sd)
% sizes
n  = size(Y,1);
m  = size(L,1);

B = sd.B;

% shorthand zeros
Znn  = zeros(n,n);   Znm  = zeros(n,m);
Zmn  = zeros(m,n);

% ---------- (12) Left big symmetric block ----------
row1 = [ Y - tDelta*(B*B') - beta*eye(n),  Znn,  B*L,    Znm ];
row2 = [ Znn,                          Znn,  Y,    Znm ];
row3 = [ L'*B',                          Y',  Y,    L' ];
row4 = [ Zmn,                          Zmn,  L, eye(m) ];

F1 = [ row1;
    row2;
    row3;
    row4 ];

% data-dependent term: G * Phi * G'
F2 = G*Phi*G';

F3 = [Y , L';
    L, eye(m) ];
end
```

---

### 3.3 Implicit Differentiation: `implicit/dtDelta_dD.m`

```matlab
function dtDelta_dD = dtDelta_dD(n,m,T,B,G,Phi,Lambda1,F1,F2,Lambda3,F3,alpha,Lambda_alpha,beta,Lambda_beta,tDelta,Lambda_tDelta,Lambda_Y,Y)
% Dの形状は D = [Z',X', U']'

G1_sol = implicit.G_grad.G1_merge(n,m,B,T);
G2_sol = implicit.G_grad.G2_merge(n,m,B,T);
G3_sol = implicit.G_grad.G3_merge(n,m,T,Lambda1,F2,B,G,Phi);
G4_sol = implicit.G_grad.G4_merge(n,m,T);
G5_sol = implicit.G_grad.G5_merge(n,m,T,B);
G6_sol = implicit.G_grad.G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi);
G7_sol = implicit.G_grad.G7_merge(n,m,T,Lambda3,F3);
G8_sol = implicit.G_grad.G8_merge(n,m,T,alpha,Lambda_alpha);
G9_sol = implicit.G_grad.G9_merge(n,m,T,beta,Lambda_beta);
G10_sol = implicit.G_grad.G10_merge(n,m,T,tDelta,Lambda_tDelta);
G11_sol = implicit.G_grad.G11_merge(n,m,T);
G12_sol = implicit.G_grad.G12_merge(n,m,T);
G13_sol = implicit.G_grad.G13_merge(n,m,T);
G14_sol = implicit.G_grad.G14_merge(n,m,T,Lambda_Y,Y);
G15_sol = implicit.G_grad.G15_merge(n,m,T,Lambda_Y);

Matrix_H = [
    G5_sol.G5_row_without_Data;
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    G8_sol.G8_row_without_Data;
    G9_sol.G9_row_without_Data;
    G10_sol.G10_row_without_Data;
    G11_sol.G11_row_without_Data;
    G12_sol.G12_row_without_Data;
    G13_sol.G13_row_without_Data;
    G14_sol.G14_row_without_Data;
    G15_sol.G15_row_without_Data];

Mtarix_y = [
    G5_sol.Data;
    G1_sol.Data;
    G2_sol.Data;
    G3_sol.Data;
    G4_sol.Data;
    G6_sol.Data;
    G7_sol.Data;
    G8_sol.Data;
    G9_sol.Data;
    G10_sol.Data;
    G11_sol.Data;
    G12_sol.Data;
    G13_sol.Data;
    G14_sol.Data;
    G15_sol.Data];

% sdpvarが含まれている可能性があるため、value()で数値型に変換
Matrix_H_full = double(full(Matrix_H));
Matrix_y_full = double(full(Mtarix_y));

rank_H = double(rank(Matrix_H_full));
n_rows = double(size(Matrix_H_full, 1));
n_cols = double(size(Matrix_H_full, 2));
fprintf('Matrix_H: size=[%d,%d], rank=%d\n', n_rows, n_cols, rank_H);
if rank_H == n_cols
    fprintf('  解は一意です (rank = 列数)\n');
else
    fprintf('  解は一意ではありません (自由度: %d = 列数 %d - rank %d)\n', n_cols - rank_H, n_cols, rank_H);
    
    % Null space分析
    fprintf('\n=== Null Space分析 ===\n');
    s = svd(Matrix_H_full);
    fprintf('SVD: 最小特異値 = %e, 最大特異値 = %e, 条件数 = %e\n', min(s), max(s), max(s)/min(s));
    fprintf('下位10個の特異値:\n');
    fprintf('  %e\n', s(end-9:end));
    
    % Null spaceを計算
    tol_null = 1e-6;
    nullvec = null(Matrix_H_full, tol_null);
    if ~isempty(nullvec)
        fprintf('Null space次元: %d (tol=%e)\n', size(nullvec, 2), tol_null);
        
        % 各変数ブロックのサイズ
        n_L = n*m;
        n_Y = n*n;
        n_alpha = 1;
        n_beta = 1;
        n_tDelta = 1;
        n_Lambda1 = (3*n+m)*(3*n+m);
        n_Lambda3 = (n+m)*(n+m);
        n_Lambda_alpha = 1;
        n_Lambda_beta = 1;
        n_Lambda_tDelta = 1;
        n_Lambda_Y = n*n;
        
        % 各ブロックのインデックス
        idx_L = 1:n_L;
        idx_Y = n_L + (1:n_Y);
        idx_alpha = n_L + n_Y + 1;
        idx_beta = n_L + n_Y + n_alpha + 1;
        idx_tDelta = n_L + n_Y + n_alpha + n_beta + 1;
        idx_Lambda1 = n_L + n_Y + n_alpha + n_beta + n_tDelta + (1:n_Lambda1);
        idx_Lambda3 = idx_Lambda1(end) + (1:n_Lambda3);
        idx_Lambda_alpha = idx_Lambda3(end) + 1;
        idx_Lambda_beta = idx_Lambda_alpha + 1;
        idx_Lambda_tDelta = idx_Lambda_beta + 1;
        idx_Lambda_Y = idx_Lambda_tDelta + (1:n_Lambda_Y);
        
        % 各null vectorの各ブロックでのノルム
        fprintf('\n各変数ブロックでのnull space成分のノルム:\n');
        for i = 1:size(nullvec, 2)
            fprintf('  Null vector %d:\n', i);
            fprintf('    dL:        %e\n', norm(nullvec(idx_L, i)));
            fprintf('    dY:        %e\n', norm(nullvec(idx_Y, i)));
            fprintf('    dAlpha:    %e\n', norm(nullvec(idx_alpha, i)));
            fprintf('    dBeta:     %e\n', norm(nullvec(idx_beta, i)));
            fprintf('    dtDelta:   %e\n', norm(nullvec(idx_tDelta, i)));
            fprintf('    dLambda1:  %e\n', norm(nullvec(idx_Lambda1, i)));
            fprintf('    dLambda3:  %e\n', norm(nullvec(idx_Lambda3, i)));
            fprintf('    dLambda_alpha: %e\n', norm(nullvec(idx_Lambda_alpha, i)));
            fprintf('    dLambda_beta:  %e\n', norm(nullvec(idx_Lambda_beta, i)));
            fprintf('    dLambda_tDelta: %e\n', norm(nullvec(idx_Lambda_tDelta, i)));
            fprintf('    dLambda_Y: %e\n', norm(nullvec(idx_Lambda_Y, i)));
        end
    else
        fprintf('Null spaceは空です (tol=%e)\n', tol_null);
    end
end
matrix_x = Matrix_H_full \ (-Matrix_y_full);

% 1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
end
```

**変数の順序**: `[dL, dY, dAlpha, dBeta, dtDelta, dLambda1, dLambda3, dLambda_alpha, dLambda_beta, dLambda_tDelta, dLambda_Y]`
- `dL`: `n*m`次元
- `dY`: `n*n`次元
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

---

### 3.4 G1: LagrangianをLで微分（∂L/∂L = 0）

```matlab
function sol = G1_merge(n,m, B,T)
% G1 は Lagrangian をLで変微分したもの
nm = n*m;
sol = struct();
sol.dL = sparse(nm,nm);
sol.dY = sparse(nm,n*n);
sol.dAlpha = sparse(nm,1);
sol.dBeta = sparse(nm,1);
sol.dtDelta = sparse(nm,1);

sol.dLambda1 = (-implicit.helper.dF1_dL(n,m, B))';
sol.dLambda3 = (-implicit.helper.dF3_dL(n,m))';
sol.Lambda_Alpha = sparse(nm,1);
sol.Lambda_Beta = sparse(nm,1);
sol.Lambda_tDelta = sparse(nm,1);
sol.Lambda_Y = sparse(nm,n*n);
sol.Data = sparse(nm,T*(2*n+m));

sol.G1_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- LagrangianをLで微分: `∂L/∂L = -dF1_dL' * vec(Lambda1) - dF3_dL' * vec(Lambda3) = 0`

---

### 3.5 G2: LagrangianをYで微分（∂L/∂Y = 0）

```matlab
function sol = G2_merge(n,m, B,T)
% G2 は Lagrangian をYで変微分したもの
nn = n*n;
nm = n*m;
sol = struct();
sol.dL = sparse(nn,nm);
sol.dY = sparse(nn,nn);
sol.dAlpha = sparse(nn,1);
sol.dBeta = sparse(nn,1);
sol.dtDelta = sparse(nn,1);

sol.dLambda1 = (-implicit.helper.dF1_dY(n,m))';
sol.dLambda3 = (-implicit.helper.dF3_dY(n,m))';
sol.Lambda_Alpha = sparse(nn,1);
sol.Lambda_Beta = sparse(nn,1);
sol.Lambda_tDelta = sparse(nn,1);
sol.Lambda_Y = speye(nn);
sol.Data = sparse(nn,T*(2*n+m));

sol.G2_row_without_Data  = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- LagrangianをYで微分: `∂L/∂Y = -dF1_dY' * vec(Lambda1) - dF3_dY' * vec(Lambda3) + vec(Lambda_Y) = 0`

---

### 3.6 G3: Lagrangianをalphaで微分（∂L/∂alpha = 0）

```matlab
function sol = G3_merge(n,m,T,Lambda1,F2,B,G,Phi)
% G3 は Lagrangian をAlphaで変微分したもの
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = sparse(1,1);

sol.dLambda1 = (F2(:))';
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = -1;
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = sparse(1,1);
sol.Lambda_Y = sparse(1,n*n);
sol.Data = Lambda1(:)'*implicit.helper.dF2_dD(n,m,B,T,G,Phi);

sol.G3_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- Lagrangianをalphaで微分: `∂L/∂alpha = vec(F2)' * vec(Lambda1) - Lambda_alpha = 0`
- データ依存項: `dF2/dD`の項が`Data`に含まれる

---

### 3.7 G4: Lagrangianをbetaで微分（∂L/∂beta = 0）

```matlab
function sol = G4_merge(n,m,T)
% G4 は Lagrangian をBetaで変微分したもの
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = sparse(1,1);

E = [speye(n);
    zeros(2*n+m,n)];
M = E*E';
sol.dLambda1 = M(:)';
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = -1;
sol.Lambda_tDelta = sparse(1,1);
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G4_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- Lagrangianをbetaで微分: `∂L/∂beta = -vec(M)' * vec(Lambda1) - Lambda_beta = 0`
- ここで、`M`は`F1`の`beta*I_n`項に対応する行列

---

### 3.8 G5: LagrangianをtDeltaで微分（∂L/∂tDelta = 0）

```matlab
function sol = G5_merge(n,m,T,B)
% G5 は Lagrangian をtDeltaで変微分したもの
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = sparse(1,1);

E = [B;
    zeros(2*n+m,m)];
M = E*E';
sol.dLambda1 = M(:)';
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = -1;
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G5_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- LagrangianをtDeltaで微分: `∂L/∂tDelta = -vec(M)' * vec(Lambda1) - Lambda_tDelta - 1 = 0`
- 最後の`-1`は目的関数`-tDelta`から来る

---

### 3.9 G6: 相補性条件 (F1-α*F2)*Lambda1' = 0 をデータで微分

```matlab
function sol = G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi)
% G6 は 相補性条件における (F1-αF2)*Lambda1' = 0
% 統一形式: (制約) * Λ^T = 0
n1_2 = (3*n+m)*(3*n+m);
n2_2 = (n+m)*(n+m);
eye = speye(3*n+m);
F1_minus_alphaF2 = F1 - alpha*F2;

sol = struct();

% d/dL [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) dvec(F1)/dL
sub_matrix = kron(Lambda1,eye);
sol.dL = sub_matrix*implicit.helper.dF1_dL(n,m,B);
sol.dY = sub_matrix*implicit.helper.dF1_dY(n,m);

% d/dalpha [vec((F1-alpha*F2)*Lambda1')] = -(Lambda1 ⊗ I) vec(F2) = -vec(F2*Lambda1')
Lambda1F2 = F2*Lambda1';
sol.dAlpha = -Lambda1F2(:);

E11 = [speye(n);
    zeros(2*n+m,n)];
M11 = E11*E11';
% d/dbeta [vec((F1-alpha*F2)*Lambda1')] = -vec(M11*Lambda1')
Lambda1M11 = M11*Lambda1';
sol.dBeta = -Lambda1M11(:);

E11_withB = [B;
    zeros(2*n+m,m)];
M11_withB = E11_withB*E11_withB';
% d/dtDelta [vec((F1-alpha*F2)*Lambda1')] = -vec(M11_withB*Lambda1')
Lambda1M11_withB = M11_withB*Lambda1';
sol.dtDelta = -Lambda1M11_withB(:);

% d/dLambda1 [vec((F1-alpha*F2)*Lambda1')] = (I ⊗ (F1-alpha*F2)) K
sol.dLambda1 = kron(eye, F1_minus_alphaF2)*implicit.helper.commutation(3*n+m,3*n+m);
sol.dLambda3 = sparse(n1_2,n2_2);
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n*n);
sol.Data = -alpha*sub_matrix*implicit.helper.dF2_dD(n,m,B,T,G,Phi);

sol.G6_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**計算内容**:
- 相補性条件`(F1 - alpha*F2) * Lambda1' = 0`をデータで微分
- **注意**: `Lambda1F2 = F2*Lambda1'`（`F2*Lambda1'`の形式）
- `d/dalpha [vec((F1-alpha*F2)*Lambda1')] = -vec(F2*Lambda1')`
- `d/dbeta [vec((F1-alpha*F2)*Lambda1')] = -vec(M11*Lambda1')`
- `d/dtDelta [vec((F1-alpha*F2)*Lambda1')] = -vec(M11_withB*Lambda1')`

---

### 3.10 G7: 相補性条件 F3*Lambda3' = 0 をデータで微分

```matlab
function sol = G7_merge(n,m,T,Lambda3,F3)
% G7 は 相補性条件における F3*Lambda3' = 0
% 統一形式: (制約) * Λ^T = 0
n2_2 = (n+m)*(n+m);
eye = speye(n+m);

sol = struct();

% vec(F3*Lambda3') = (Lambda3 ⊗ I) vec(F3)
sub_matrix = kron(Lambda3,eye);
sol.dL = sub_matrix*implicit.helper.dF3_dL(n,m);
sol.dY = sub_matrix*implicit.helper.dF3_dY(n,m);
sol.dAlpha = sparse(n2_2,1);
sol.dBeta = sparse(n2_2,1);
sol.dtDelta = sparse(n2_2,1);

sol.dLambda1 = sparse(n2_2,(3*n+m)*(3*n+m));
% d/dLambda3 [vec(F3*Lambda3')] = (I ⊗ F3) K
sol.dLambda3 = kron(eye,F3)*implicit.helper.commutation(n+m,n+m);
sol.Lambda_Alpha = sparse(n2_2,1);
sol.Lambda_Beta = sparse(n2_2,1);
sol.Lambda_tDelta = sparse(n2_2,1);
sol.Lambda_Y = sparse(n2_2,n*n);
sol.Data = sparse(n2_2,T*(2*n+m));

sol.G7_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.11 G8: 相補性条件 alpha*Lambda_alpha = 0 をデータで微分

```matlab
function sol = G8_merge(n,m,T,alpha,Lambda_alpha)
% G8 は 相補性条件における αΛ_alpha = 0
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = Lambda_alpha;
sol.dBeta = sparse(1,1);
sol.dtDelta = sparse(1,1);

sol.dLambda1 = sparse(1,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = alpha;
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = sparse(1,1);
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G8_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.12 G9: 相補性条件 (beta-ε)*Lambda_beta = 0 をデータで微分

```matlab
function sol = G9_merge(n,m,T,beta,Lambda_beta)
% G9 は 相補性条件における βΛ_beta = 0
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = Lambda_beta;
sol.dtDelta = sparse(1,1);

sol.dLambda1 = sparse(1,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = beta;
sol.Lambda_tDelta = sparse(1,1);
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G9_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.13 G10: 相補性条件 (tDelta-ε)*Lambda_tDelta = 0 をデータで微分

```matlab
function sol = G10_merge(n,m,T,tDelta,Lambda_tDelta)
% G10 は 相補性条件における tΔΛ_tΔ = 0
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = Lambda_tDelta;

sol.dLambda1 = sparse(1,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = tDelta;
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G10_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.14 G11: Lambda1の対称性 Lambda1 - Lambda1' = 0

```matlab
function sol = G11_merge(n,m,T)
% G11 は 双対Λの対称性 Λ1 - Λ1' = 0
n1_2 = (3*n+m)*(3*n+m);
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = sparse(n1_2,n*n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = speye(n1_2) - implicit.helper.commutation(3*n+m,3*n+m);
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n*n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G11_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.15 G12: Lambda3の対称性 Lambda3 - Lambda3' = 0

```matlab
function sol = G12_merge(n,m,T)
% G12 は 双対Λの対称性 Λ3 - Λ3' = 0
n1_2 = (n+m)*(n+m);
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = sparse(n1_2,n*n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = speye(n1_2)-implicit.helper.commutation(n+m,n+m);
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n*n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G12_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.16 G13: Yの対称性 Y - Y' = 0

```matlab
function sol = G13_merge(n,m,T)
% G13 は 双対Λの対称性 Y - Y' = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = speye(n1_2)-implicit.helper.commutation(n,n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n1_2);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G13_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

**注意**: `Y`は`sdpvar(n,n,'symmetric')`として定義されているため、この制約は理論的には冗長

---

### 3.17 G14: 相補性条件 Y*Lambda_Y' = 0 をデータで微分

```matlab
function sol = G14_merge(n,m,T,Lambda_Y,Y)
% G14 は 相補性条件における Y*Lambda_Y' = 0
% 統一形式: (制約) * Λ^T = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);

% d/dY [vec(Y*Lambda_Y')] = Lambda_Y ⊗ I
sol.dY = kron(Lambda_Y,eye(n));
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
% d/dLambda_Y [vec(Y*Lambda_Y')] = (I ⊗ Y) K
sol.Lambda_Y = kron(speye(n),Y)*implicit.helper.commutation(n,n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G14_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.18 G15: Lambda_Yの対称性 Lambda_Y - Lambda_Y' = 0

```matlab
function sol = G15_merge(n,m,T,Lambda_Y)
% G15 は 双対Λの対称性 Lambda_Y - Lambda_Y' = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = sparse(n1_2,n*n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = speye(n1_2)-implicit.helper.commutation(n,n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G15_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end
```

---

### 3.19 Helper関数: dF1_dL.m

```matlab
function dF1_dL = dF1_dL(n,m, B)
E1 = [speye(n);
    zeros(2*n+m,n)];

E3= [zeros(2*n,n);
    speye(n);
    zeros(m,n)];

E4 = [zeros(3*n,m);
    speye(m)];

term1 = kron(E3,E1*B);
term2 = kron(E1*B,E3)*implicit.helper.commutation(m,n);
term3 = kron(E3,E4);
term4 = kron(E4,E3)*implicit.helper.commutation(m,n);

dF1_dL = term1 + term2 + term3 + term4;
end
```

**計算内容**:
- `vec(F1)`を`vec(L)`で微分した行列
- F1の各ブロックがLに依存する部分（`B*L`, `L'*B'`, `L`, `L'`）を抽出

---

### 3.20 Helper関数: dF1_dY.m

```matlab
function dF1_dY = dF1_dY(n,m)
E1 = [speye(n);
    zeros(2*n+m,n)];

E2 = [zeros(n,n);
    speye(n);
    zeros(n+m,n)];

E3 = [zeros(n,n);
    zeros(n,n);
    speye(n);
    zeros(m,n)];
dF1_dY = kron(E1,E1) + kron(E3,E3) + kron(E3,E2) + kron(E2,E3)*implicit.helper.commutation(n,n);
end
```

**計算内容**:
- `vec(F1)`を`vec(Y)`で微分した行列
- F1の各ブロックがYに依存する部分を抽出

---

### 3.21 Helper関数: dF2_dD.m

```matlab
function dF2_dD = dF2_dD(n,m,B,T,G,Phi)
Ez = [speye(n), zeros(n,n), zeros(n,m)];
Ex = [zeros(n,n), speye(n), zeros(n,m)];
Eu = [zeros(m,n), zeros(m,n), speye(m)];

n1 = 3*n+m;
n2 = n+T;

Eleft = [sparse(n,T);
    speye(T)];

Eright = [Ez - B*Eu;
    -Ex;
    sparse(n+m,2*n+m)];

GPhi = G*Phi;

dG_dD = kron(Eleft,Eright);

dF2_dD = kron(GPhi,speye(n1))*dG_dD + kron(speye(n1), GPhi)*implicit.helper.commutation(n1,n2)*dG_dD;
end
```

**計算内容**:
- `F2 = G*Phi*G'`をデータ`D`で微分
- `dF2/dD = d(G*Phi*G')/dD = (G*Phi ⊗ I) * dG/dD + (I ⊗ G*Phi) * K * dG/dD`

---

### 3.22 Helper関数: dF3_dL.m

```matlab
function dF3_dL = dF3_dL(n,m)

E1 = [speye(n);
    zeros(m,n)];

E2 = [zeros(n,m);
    speye(m)];

dF3_dL = kron(E1, E2) + kron(E2, E1) * implicit.helper.commutation(m,n);
end
```

---

### 3.23 Helper関数: dF3_dY.m

```matlab
function dF3_dY = dF3_dY(n,m)
E1 = [speye(n);
    zeros(m,n)];

dF3_dY = kron(E1,E1);
end
```

---

## 4. KKT条件の構造

### 4.1 変数の順序

`Matrix_H`の列は以下の順序で並んでいます：

1. `dL`: `n*m`次元
2. `dY`: `n*n`次元
3. `dAlpha`: 1次元
4. `dBeta`: 1次元
5. `dtDelta`: 1次元
6. `dLambda1`: `(3*n+m)^2`次元
7. `dLambda3`: `(n+m)^2`次元
8. `dLambda_alpha`: 1次元
9. `dLambda_beta`: 1次元
10. `dLambda_tDelta`: 1次元
11. `dLambda_Y`: `n*n`次元

**合計列数**: `n*m + n*n + 1 + 1 + 1 + (3*n+m)^2 + (n+m)^2 + 1 + 1 + 1 + n*n = 324` (n=4, m=3の場合)

### 4.2 各行の対応

| 行 | KKT条件 | 行数 |
|---|---------|------|
| G5 | `∂L/∂tDelta = 0` | 1 |
| G1 | `∂L/∂L = 0` | `n*m` |
| G2 | `∂L/∂Y = 0` | `n*n` |
| G3 | `∂L/∂alpha = 0` | 1 |
| G4 | `∂L/∂beta = 0` | 1 |
| G6 | `(F1-α*F2)*Lambda1' = 0` (相補性) | `(3*n+m)^2` |
| G7 | `F3*Lambda3' = 0` (相補性) | `(n+m)^2` |
| G8 | `alpha*Lambda_alpha = 0` (相補性) | 1 |
| G9 | `(beta-ε)*Lambda_beta = 0` (相補性) | 1 |
| G10 | `(tDelta-ε)*Lambda_tDelta = 0` (相補性) | 1 |
| G11 | `Lambda1 - Lambda1' = 0` (対称性) | `(3*n+m)^2` |
| G12 | `Lambda3 - Lambda3' = 0` (対称性) | `(n+m)^2` |
| G13 | `Y - Y' = 0` (対称性) | `n*n` |
| G14 | `Y*Lambda_Y' = 0` (相補性) | `n*n` |
| G15 | `Lambda_Y - Lambda_Y' = 0` (対称性) | `n*n` |

**合計行数**: `1 + n*m + n*n + 1 + 1 + (3*n+m)^2 + (n+m)^2 + 1 + 1 + 1 + (3*n+m)^2 + (n+m)^2 + n*n + n*n + n*n = 630` (n=4, m=3の場合)

---

## 5. 確認したいこと

以下の点について、計算やコーディングミスがないか確認していただきたいです：

### 5.1 数学的な正しさ

1. **G6 (相補性条件 (F1-α*F2)*Lambda1' = 0)**:
   - `Lambda1F2 = F2*Lambda1'`は正しいか？（`F2*Lambda1'`の形式）
   - `Lambda1M11 = M11*Lambda1'`は正しいか？
   - `Lambda1M11_withB = M11_withB*Lambda1'`は正しいか？
   - これらは`(F1-α*F2)*Lambda1' = 0`の形式と整合しているか？

2. **G1, G2 (最適性条件)**:
   - `dF1_dL`と`dF1_dY`の計算は正しいか？
   - `dF3_dL`と`dF3_dY`の計算は正しいか？

3. **G3 (∂L/∂alpha = 0)**:
   - `sol.dLambda1 = (F2(:))'`は正しいか？
   - データ依存項`dF2_dD`の計算は正しいか？

4. **G5 (∂L/∂tDelta = 0)**:
   - `sol.dLambda1 = M(:)'`は正しいか？（Mは`B*B'`の項に対応）

5. **G13, G14, G15 (Yに関する制約)**:
   - `Y`が`sdpvar(n,n,'symmetric')`として定義されている場合、`G13`の`Y - Y' = 0`は冗長ではないか？
   - `G14`と`G15`の相互作用は正しいか？

### 5.2 データ依存項の計算

1. **dF2_dD.m**:
   - `F2 = G*Phi*G'`の微分は正しいか？
   - `dG/dD`の計算は正しいか？
   - Commutation行列の使用は正しいか？

### 5.3 次元の整合性

1. **各G_merge関数の出力サイズ**:
   - すべての行の列数が一致しているか？
   - `Data`項のサイズが`T*(2*n+m)`になっているか？

### 5.4 符号の確認

1. **Lagrangianの符号**:
   - `L = -tDelta - tr((F1-alpha*F2)*Lambda1') - tr((F3-ε*I)*Lambda3') - ...`
   - 各項の符号は正しいか？

2. **微分の符号**:
   - `∂L/∂L = -dF1_dL' * vec(Lambda1) - dF3_dL' * vec(Lambda3)`の符号は正しいか？

### 5.5 Commutation行列の使用

1. **正しい使用**:
   - `vec(A^T) = K * vec(A)`の関係を正しく使用しているか？
   - Commutation行列のサイズは正しいか？

### 5.6 対称性の扱い

1. **Yの対称性**:
   - `G13`で`Y - Y' = 0`を追加するのは冗長ではないか？（`Y`は既に対称として定義されている）

2. **Lambda1, Lambda3, Lambda_Yの対称性**:
   - 各Lambdaの対称性制約は正しく実装されているか？

---

## 6. 参考情報

### 6.1 理論的背景

詳細な理論的説明は`doc/kkt_implicit_differentiation.tex`を参照してください。

### 6.2 重要な関係式

1. **Kronecker積とvec演算**:
   - `vec(ABC) = (C^T ⊗ A) * vec(B)`
   - `vec(AB^T) = (B ⊗ A) * vec(I)`（適切なcommutation行列を使用）

2. **Commutation行列**:
   - `vec(A^T) = K_{p,q} * vec(A)`（Aがp×q行列の場合）

3. **相補性条件の微分**:
   - `d/dx [vec((Constraint)*Lambda^T)] = (Lambda ⊗ I) * dvec(Constraint)/dx + (I ⊗ Constraint) * K * dvec(Lambda)/dx`

---

## 7. お願い

この実装の**計算やコーディングミスがないか**を確認していただきたいです。特に、以下の点について詳しくご確認いただけると助かります：

1. **数学的な正しさ**: 各G_merge関数の計算式が理論的に正しいか
2. **次元の整合性**: すべての行列のサイズが一致しているか
3. **符号の確認**: Lagrangianの符号と微分の符号が正しいか
4. **Commutation行列の使用**: 正しく使用されているか
5. **データ依存項の計算**: `dF2_dD`などの計算が正しいか
6. **G6の計算**: `F2*Lambda1'`の形式が正しいか（`Lambda1*F2`ではないか？）

ご回答をお待ちしております。





