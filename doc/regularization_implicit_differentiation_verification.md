# 正則化付きSDPのImplicit Differentiation実装の検証依頼

## 1. やりたいこと（目的）

正則化付き量子化制御器設計のSDP問題に対して、**KKT条件を用いたImplicit Differentiation**により、**delta（量子化パラメータ）のデータに対する勾配**を計算する実装を行っています。

具体的には、データ`D = [Z', X', U']'`（Z: 次状態、X: 状態、U: 入力）に対する`delta`の勾配`d(delta)/dD`を計算したいです。

この実装の**計算やコーディングミスがないか**を確認していただきたいです。

---

## 2. SDP問題の定式化

### 2.1 決定変数

- `L ∈ ℝ^(T×n)`: 制御パラメータ行列
- `delta ∈ ℝ`: 量子化パラメータ（スカラー、delta > 0）

### 2.2 制約条件

1. **LMI制約**: `F1 >= 0`
2. **正定値制約**: `(X*L)_sym >= tolerance*eye(n)` （`(X*L)_sym = (X*L + (X*L)')/2`）
3. **対称性制約**: `X*L == (X*L)'` （等号制約）

### 2.3 目的関数

```
minimize: -delta + gamma * ||Pi*L||_F^2
```

ここで：
- `Pi = eye(T) - pinv(Gamma)*Gamma`
- `Gamma = [U; X]`
- `gamma > 0`: 正則化パラメータ（通常1e3〜1e5）

### 2.4 LMIブロック F1 の定義

```matlab
Xm_L_Sym = (Xm*L + (Xm*L)')/2;
const_mat = [Xm_L_Sym, (Z*L)', zeros(n,m), (Um*L)';
             Z*L, Xm_L_Sym, delta*B, zeros(n,m);
             zeros(m,n), delta*B', eye(m), zeros(m,m);
             Um*L, zeros(m,n), zeros(m,m), eye(m)];
```

**サイズ**: `(2n+2m) × (2n+2m)`

---

## 3. 実装コード

### 3.1 SDP解法: `regularization_sdp/solve_sdp.m`

```matlab
function [sol, K, delta_val, L_val, diagnostics] = solve_sdp(data, gamma, opts)

if nargin < 2 || isempty(gamma), gamma = 1e5; end
if nargin < 3 || isempty(opts), opts = struct(); end
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

[n, T] = size(Z);
m = size(Um,1);

% -------------------------
% Decision variables
% -------------------------
L      = sdpvar(T,n,'full');
delta = sdpvar(1,1);     % = δ

% -------------------------
% LMI blocks
% -------------------------
Xm_L_Sym = (Xm*L + (Xm*L)')/2;
const_mat = [Xm_L_Sym, (Z*L)', zeros(n,m),(Um*L)';
    Z*L, Xm_L_Sym, delta*B, zeros(n,m);
    zeros(m,n), delta*B', eye(m), zeros(m,m);
    Um*L, zeros(m,n), zeros(m,m), eye(m)];

% -------------------------
% Constraints
% -------------------------

tolerance = 0;

% 制約を個別に定義（dual取得のため）
const_1 = [const_mat >= 0];
const_2 = [Xm_L_Sym >= (tolerance)*eye(n)];
const_3 = [Xm*L == (Xm*L)'];

constr  = [];
constr  = [constr, const_1];
constr  = [constr, const_2];
constr  = [constr, const_3];
% constr  = [constr, delta >= tolerance];
% Objective: maximize δ  <=>  minimize -delta
pi_matrix = eye(T) - pinv(Gamma_Matrix)*Gamma_Matrix;
obj = -delta + gamma*norm(pi_matrix*L, 'fro')^2;

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
L_val = value(L);
delta_val = value(delta);
sol.L = L_val;
rho_val = (1-delta_val)/(1+delta_val);
sol.rho = rho_val;
sol.objective = value(obj);
sol.delta = delta_val;
sol.F1 = value(const_mat);

K = Um*L_val*(Xm*L_val)^(-1);
sol.K = K;

% -------------------------
% Dual values (Lagrange multipliers)
% -------------------------
try
    sol.Lambda1 = dual(const_1);  % F1 >= 0 のdual（サイズ: (2n+2m) × (2n+2m)）
    sol.Lambda2 = dual(const_2);  % Xm_L_Sym >= tolerance*eye(n) のdual（サイズ: n × n）
    sol.Lambda3 = dual(const_3);  % Xm*L == (Xm*L)' のdual（サイズ: n × n、等号制約なので符号なし）
catch ME
    error('solve_sdp:DualFailed', ...
        'Dual変数の取得に失敗しました: %s', ME.message);
end
end
```

**注意**: 
- `Lambda2`は`Lambda_P`としても参照される（相補性条件で使用）
- `Lambda3`は等号制約の双対変数（符号の制約なし）

---

### 3.2 Implicit Differentiation: `implicit_regularization/dtDelta_dD.m`

```matlab
function dtDelta_dD = dtDelta_dD(n,m,T,B,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,F1,Lambda3)
% Dの形状は D = [Z',X', U']'

G1_sol = implicit_regularization.G_grad.G1_merge(n,m,B,T);
G2_sol = implicit_regularization.G_grad.G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,Lambda3);
G3_sol = implicit_regularization.G_grad.G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B);
G4_sol = implicit_regularization.G_grad.G4_merge(n,m,T);
G5_sol = implicit_regularization.G_grad.G5_merge(n,m,T,X,Lambda_P,L);
G6_sol = implicit_regularization.G_grad.G6_merge(n,m,T);
G7_sol = implicit_regularization.G_grad.G7_merge(n,m,T,X,L);
G8_sol = implicit_regularization.G_grad.G8_merge(n,m,T);

Matrix_H = [
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G5_sol.G5_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    G8_sol.G8_row_without_Data;
    ];

Mtarix_y = [
    G1_sol.Data;
    G2_sol.Data;
    G3_sol.Data;
    G4_sol.Data;
    G5_sol.Data;
    G6_sol.Data;
    G7_sol.Data;
    G8_sol.Data;
    ];

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
end

matrix_x = Matrix_H_full \ (-Matrix_y_full);

% 1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
end
```

**変数の順序**: `[dtDelta, dL, dLambda, dLambda_P, dLambda3]`
- `dtDelta`: 1次元
- `dL`: `T*n`次元（vec(L)）
- `dLambda`: `(2n+2m)^2`次元（vec(Lambda1)）
- `dLambda_P`: `n*n`次元（vec(Lambda_P)）
- `dLambda3`: `n*n`次元（vec(Lambda3)）

**合計列数**: `1 + T*n + (2n+2m)^2 + n*n + n*n`

---

### 3.3 G1: Lagrangianをdeltaで微分（∂L/∂delta = 0）

```matlab
function sol = G1_merge(n,m, B,T)
% G1 は Lagrangian をdeltaで変微分したもの
sol = struct();
sol.dL = sparse(1,n*T);
dF_dDelta = [sparse(n,n), sparse(n,n), sparse(n,m), sparse(n,m);
    sparse(n,n), sparse(n,n), B, sparse(n,m);
    sparse(m,n), B', sparse(m,m), sparse(m,m);
    sparse(m,n), sparse(m,n), sparse(m,m), sparse(m,m)];

sol.dLambda = -dF_dDelta(:)';
sol.dtDelta = 0;
sol.dLambda_P = sparse(1,n*n);
sol.dLambda3 = sparse(1, n*n);
sol.Data = sparse(1,(2*n+m)*T);
sol.G1_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- Lagrangianをdeltaで微分: `∂L/∂delta = -1 - vec(dF_dDelta)' * vec(Lambda) = 0`
- ここで、`dF_dDelta`は`F1`のdelta依存部分（`delta*B`の項）

---

### 3.4 G2: LagrangianをLで微分（∂L/∂L = 0）

```matlab
function sol = G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,Lambda3)
% G2 は Lagrangian をLで変微分したもの
sol = struct();
sol.dL = 2*gamma*kron(speye(n),Pi);
dF_dL = implicit_regularization.helper.dF_dL(n,m,T,X,Z,U);
dF_dLambda = -dF_dL';

sol.dLambda = dF_dLambda;
sol.dtDelta = sparse(T*n,1);

sol.dLambda_P = -kron(speye(n),X');
% 後で
dPi_dD = implicit_regularization.helper.dPi_dD(n,m,T,Gamma);
term1 = 2*gamma*kron(L',speye(T))*dPi_dD;
term2 = -implicit_regularization.helper.dVec_Lambda_TdFdL_dD(n,m,T,Lambda);

C2nmT = implicit.helper.commutation(2*n+m,T);
Ex = [sparse(n,n), speye(n),sparse(n,m)];
term3 = - (kron(Lambda_P'*Ex,speye(T)))*C2nmT;

Knn = implicit.helper.commutation(n,n);
sol.dLambda3 = kron(speye(n), X') * (speye(n*n) - Knn);

sol.Data = term1 + term2 + term3;
% --- new: Data term due to X dependence of X'*(Lambda3-Lambda3') ---
M   = Lambda3 - Lambda3';
KTx = implicit.helper.commutation(T,n);        % vec(dX') = KTx * vec(dX)
Ex  = [sparse(n,n), speye(n), sparse(n,m)];    % picks X from D=[Z;X;U]
extra_Data = kron(M', speye(T)) * KTx * kron(speye(T), Ex);

sol.Data = sol.Data + extra_Data;

sol.G2_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- LagrangianをLで微分: `∂L/∂L = 2*gamma*Pi*L - dF_dL' * vec(Lambda) - X' * vec(Lambda_P) + X'*(Lambda3 - Lambda3') = 0`
- 等号制約`X*L = (X*L)'`の双対変数`Lambda3`による項: `X'*(Lambda3 - Lambda3')`
- `vec(X'*(Lambda3 - Lambda3')) = kron(I_n, X') * (I - K_nn) * vec(Lambda3)`
- データ依存項: `dPi/dD`、`dF/dD`、および`X`依存項（`Lambda3`の項から）

---

### 3.5 G3: 相補性条件 F1*Lambda1' = 0 をデータで微分

```matlab
function sol = G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B)
% G3 は 相補性条件 をdeltaで変微分したもの
n1 = 2*n+2*m;
sol = struct();
Left = kron(Lambda',speye(n1));
dF_dL = implicit_regularization.helper.dF_dL(n,m,T,X,Z,U);
sol.dL = Left*dF_dL;

sol.dLambda = kron(speye(n1),F1)*implicit.helper.commutation(n1,n1);

% dF/dDelta（スカラーなので vec(∂F/∂δ) を作る）
dF_dDelta = implicit_regularization.helper.dF_dDelta(n,m,B); % (n1^2 x 1)
sol.dtDelta = Left * dF_dDelta;  % ★ここが本質

sol.dLambda_P = sparse(n1*n1,n*n);
sol.dLambda3 = sparse(n1*n1, n*n);
sol.Data = Left*implicit_regularization.helper.dF_dD(n,m,T,L);
sol.G3_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- 相補性条件`F1 * Lambda1' = 0`をデータで微分
- `vec(F1*Lambda1') = (Lambda1 ⊗ I) * vec(F1)`
- `d/ddelta [vec(F1*Lambda1')] = (Lambda1 ⊗ I) * dF_dDelta`

---

### 3.6 G4: Lambda1の対称性 Lambda1 - Lambda1' = 0

```matlab
function sol = G4_merge(n,m,T)
% G4 は Lambda - Lambda' = 0 をLで変微分したもの
n1 = 2*n+2*m;
n1_2 = n1*n1;
sol = struct();
sol.dL = sparse(n1_2,n*T);
sol.dLambda = speye(n1_2) - implicit.helper.commutation(n1,n1);
sol.dtDelta = sparse(n1_2,1);
sol.dLambda_P = sparse(n1_2,n*n);
sol.dLambda3 = sparse(n1*n1, n*n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G4_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- `vec(Lambda1 - Lambda1') = (I - K) * vec(Lambda1) = 0`
- ここで、`K`はcommutation行列

---

### 3.7 G5: 相補性条件 (X*L)_sym * Lambda_P' = 0 をデータで微分

```matlab
function sol = G5_merge(n,m,T,X,Lambda_P,L)
% G5 は Lambda_P の相補性条件
n1 = 2*n+2*m;
sol = struct();
sol.dL = 0.5*kron(Lambda_P,X) + 0.5*kron(Lambda_P*X,speye(n))*implicit.helper.commutation(T,n);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = kron(speye(n),X*L)*implicit.helper.commutation(n,n);

Ex = [sparse(n,n), speye(n),sparse(n,m)];
sol.Data = 0.5*kron(Lambda_P*L',Ex) + 0.5*kron(Lambda_P*Ex,L')*implicit.helper.commutation(2*n+m,T);
sol.dLambda3 = sparse(n*n, n*n);
sol.G5_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- 相補性条件`(X*L)_sym * Lambda_P' = 0`をデータで微分
- `(X*L)_sym = (X*L + (X*L)')/2`の対称化を考慮

---

### 3.8 G6: Lambda_Pの対称性 Lambda_P - Lambda_P' = 0

```matlab
function sol = G6_merge(n,m,T)
% G6 は Lambda_P = Lambda_P'
n1 = 2*n+2*m;
sol = struct();
sol.dL = sparse(n*n,n*T);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = speye(n*n)-implicit.helper.commutation(n,n);
sol.dLambda3 = sparse(n*n, n*n);
sol.Data = sparse(n*n,T*(2*n+m));
sol.G6_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

---

### 3.9 G7: 対称性制約 X*L = (X*L)' をデータで微分

```matlab
function sol = G7_merge(n,m,T,X,L)
% G7 は XL = (XL)'
n1 = 2*n+2*m;
sol = struct();
sol.dL = kron(speye(n),X) - kron(X,speye(n))*implicit.helper.commutation(T,n);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = sparse(n*n,n*n);
C2nmT = implicit.helper.commutation(2*n+m,T);
Ex = [sparse(n,n), speye(n),sparse(n,m)];
sol.Data = kron(L',Ex)- kron(Ex,L')*C2nmT;
sol.dLambda3 = sparse(n*n, n*n);
sol.G7_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- 等号制約`X*L = (X*L)'`をデータで微分
- `vec(X*L - (X*L)') = (I ⊗ X - X ⊗ I * K) * vec(L) + データ依存項 = 0`
- この制約の双対変数は`Lambda3`（G2で使用）

---

### 3.9 G8: Lambda3のゲージ固定 Lambda3 + Lambda3' = 0

```matlab
function sol = G8_merge(n, m, T)
sol = struct();

Knn = implicit.helper.commutation(n,n);
r = n*n; % rows

sol.dtDelta   = sparse(r,1);
sol.dL        = sparse(r,T*n);

n1 = 2*n+2*m;
sol.dLambda   = sparse(r,n1*n1);
sol.dLambda_P = sparse(r,n*n);

sol.dLambda3  = speye(r) + Knn;

sol.Data      = sparse(r, T*(2*n+m));

sol.G8_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
```

**計算内容**:
- `Lambda3`のゲージ固定条件: `Lambda3 + Lambda3' = 0`（歪対称）
- `vec(Lambda3 + Lambda3') = (I + K_nn) * vec(Lambda3) = 0`
- この条件により、`Lambda3`の対称成分の自由度を除去し、rank deficiencyを解消

---

### 3.10 G8: Lambda3のゲージ固定 Lambda3 + Lambda3' = 0

**注意**: G8は`Lambda3`の対称成分を固定するためのゲージ条件です。等号制約`X*L = (X*L)'`の双対変数`Lambda3`は対称成分が自由（ゲージ自由度）であるため、rank deficiencyが発生します。G8によりこの自由度を除去します。

---

### 3.11 Helper関数: dF_dL.m

```matlab
function dF_dL = dF_dL(n,m,T,X,Z,U)
E1 = [speye(n);
    sparse(n+2*m,n)];
E2 = [sparse(n,n);
    speye(n);
    sparse(2*m,n)];
X1 = [X;
    sparse(2*m+n,T)];
X2 = [sparse(n,T);
    X;
    sparse(2*m,T)];
Z2 = [sparse(n,T);
    Z;
    sparse(2*m,T)];
U4 = [sparse(2*n+m,T);
    U];
dF_dL = kron(E1,X1)+ kron(E2,X2) + kron(E2,Z2) +...
    kron(E1,Z2) +kron(Z2,E1)*implicit.helper.commutation(T,n) +...
    kron(E1,U4) +kron(U4,E1)*implicit.helper.commutation(T,n);
end
```

**計算内容**:
- `vec(F1)`を`vec(L)`で微分した行列
- F1の各ブロックがLに依存する部分（`X*L`, `Z*L`, `U*L`など）を抽出

---

### 3.12 Helper関数: dF_dDelta.m

```matlab
function v = dF_dDelta(n,m,B)
n1 = 2*n+2*m;

dF = sparse(n1,n1);

r2 = (n+1):(2*n);
c2 = (n+1):(2*n);
r3 = (2*n+1):(2*n+m);
c3 = (2*n+1):(2*n+m);

% δB block: (row2, col3)
dF(r2, c3) = B;

% δB^T block: (row3, col2)
dF(r3, c2) = B';

v = dF(:); % vec(∂F/∂δ)
end
```

**計算内容**:
- `F1`のdelta依存部分（`delta*B`の項）を抽出
- F1の(2,3)ブロックと(3,2)ブロックのみがdeltaに依存

---

### 3.13 Helper関数: dF_dD.m

```matlab
function dF_dD = dF_dD(n,m,T,L)
n1 = 2*n+m;

L1 = [L';
    sparse(n+2*m,T)];
L2 = [sparse(n,T);
    L';
    sparse(2*m,T)];

Ez = [speye(n,n),sparse(n,n+m)];
Ex = [sparse(n,n),speye(n),sparse(n,m)];
Eu = [sparse(m,n),sparse(m,n),speye(m)];

Ex1 = [Ex;
    sparse(n+m*2,n1)];
Ex2 = [sparse(n,n1);
    Ex;
    sparse(m+m,n1)];
Ez2 = [sparse(n,n1);
    Ez;
    sparse(m+m,n1)];
Eu4 = [sparse(2*n+m,n1);
    Eu];

term1 = kron(L1,Ex1);
term2 = kron(L2,Ex2);
term3 = kron(L1,Ez2);
term4 = kron(Ez2,L1)*implicit.helper.commutation(n1,T);
term5 = kron(L1,Eu4);
term6 = kron(Eu4,L1)*implicit.helper.commutation(n1,T);

dF_dD = term1 + term2 + term3 + term4 + term5 + term6;
end
```

**計算内容**:
- `vec(F1)`をデータ`D = [Z', X', U']'`で微分した行列
- F1の各ブロックがデータに依存する部分（`X*L`, `Z*L`, `U*L`など）を抽出

---

### 3.14 Helper関数: dPi_dXU.m

```matlab
function [dPi_dX,dPi_dU] = dPi_dXU(n,m,T,Gamma)
tempX = [sparse(m,n);speye(n)];
tempU = [speye(m);sparse(n,m)];
dGamma_dX = kron(speye(T),tempX);
dGamma_dU = kron(speye(T),tempU);

pinvGamma = pinv(Gamma);
dpinvGamma_dGamma = implicit_regularization.helper.dpinvGamma_dGamma(n,m,T,Gamma);

termX1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dX;
termX2 = -kron(speye(T),pinvGamma)*dGamma_dX;
dPi_dX = termX1 + termX2;

termU1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dU;
termU2 = -kron(speye(T),pinvGamma)*dGamma_dU;
dPi_dU = termU1 + termU2;
end
```

**計算内容**:
- `Pi = eye(T) - pinv(Gamma)*Gamma`を`X`と`U`で微分
- `Gamma = [U; X]`なので、`dGamma/dX`と`dGamma/dU`を計算
- `dPi/dGamma = -d(pinv(Gamma))/dGamma * Gamma - pinv(Gamma) * dGamma/dGamma`

---

### 3.15 Helper関数: dVec_Lambda_TdFdL_dD.m

```matlab
function dVec_Lambda_TdFdL_dD = dVec_Lambda_TdFdL_dD(n,m,T,Lambda)

n1 = 2*n+m;
Lambda11 = Lambda(1:n,1:n);
Lambda21 = Lambda(n+1:n+n,1:n);
Lambda22 = Lambda(n+1:n+n,n+1:n+n);
Lambda41 = Lambda(2*n+m+1:2*n+2*m,1:n);

Ez = [speye(n,n),sparse(n,n+m)];
Ex = [sparse(n,n),speye(n),sparse(n,m)];
Eu = [sparse(m,n),sparse(m,n),speye(m)];

LambdaX = (Lambda11 + Lambda22)'*Ex;
LambdaZ = (Lambda21')*Ez;
LambdaU = (Lambda41')*Eu;

I_T = speye(T);

Cn1T = implicit.helper.commutation(n1,T);

term1 = kron(LambdaX,I_T)*Cn1T;
term2 = 2*kron(LambdaZ,I_T)*Cn1T;
term3 = 2*kron(LambdaU,I_T)*Cn1T;

dVec_Lambda_TdFdL_dD = term1 + term2 + term3;
end
```

**計算内容**:
- `vec(Lambda' * dF_dL)`をデータ`D`で微分
- `Lambda`をブロックに分解し、各ブロックがデータに依存する部分を抽出

---

## 4. KKT条件の構造

### 4.1 変数の順序

`Matrix_H`の列は以下の順序で並んでいます：

1. `dtDelta`: 1次元
2. `dL`: `T*n`次元（`vec(L)`）
3. `dLambda`: `(2n+2m)^2`次元（`vec(Lambda1)`）
4. `dLambda_P`: `n*n`次元（`vec(Lambda_P)`）
5. `dLambda3`: `n*n`次元（`vec(Lambda3)`）

**合計列数**: `1 + T*n + (2n+2m)^2 + n*n + n*n`

### 4.2 各行の対応

| 行 | KKT条件 | 行数 |
|---|---------|------|
| G1 | `∂L/∂delta = 0` | 1 |
| G2 | `∂L/∂L = 0` | `T*n` |
| G3 | `F1 * Lambda1' = 0` (相補性) | `(2n+2m)^2` |
| G4 | `Lambda1 - Lambda1' = 0` (対称性) | `(2n+2m)^2` |
| G5 | `(X*L)_sym * Lambda_P' = 0` (相補性) | `n*n` |
| G6 | `Lambda_P - Lambda_P' = 0` (対称性) | `n*n` |
| G7 | `X*L = (X*L)'` (等号制約) | `n*n` |
| G8 | `Lambda3 + Lambda3' = 0` (ゲージ固定) | `n*n` |

**合計行数**: `1 + T*n + (2n+2m)^2 + (2n+2m)^2 + n*n + n*n + n*n + n*n`

---

## 5. 確認したいこと

以下の点について、計算やコーディングミスがないか確認していただきたいです：

### 5.1 数学的な正しさ

1. **G1 (∂L/∂delta = 0)**:
   - `dF_dDelta`の定義は正しいか？（F1のdelta依存部分の抽出）
   - `sol.dLambda = -dF_dDelta(:)'`は正しいか？

2. **G2 (∂L/∂L = 0)**:
   - 正則化項`gamma * ||Pi*L||_F^2`の微分`2*gamma*Pi*L`は正しいか？
   - `dF_dL`の計算は正しいか？
   - `dPi/dD`の計算（`dPi_dXU.m`）は正しいか？
   - `term3`の計算（`Lambda_P`に関する項）は正しいか？
   - `Lambda3`に関する項`X'*(Lambda3 - Lambda3')`の計算は正しいか？
   - `sol.dLambda3 = kron(speye(n), X') * (speye(n*n) - Knn)`は正しいか？
   - `extra_Data`項（`X`依存項）の計算は正しいか？

3. **G3 (相補性条件 F1*Lambda1' = 0)**:
   - `Left = kron(Lambda',speye(n1))`は正しいか？（`vec(F1*Lambda1') = (Lambda1 ⊗ I) * vec(F1)`）
   - `sol.dtDelta = Left * dF_dDelta`は正しいか？（deltaに関する項）
   - `sol.dLambda = kron(speye(n1),F1)*commutation(n1,n1)`は正しいか？（`vec(F1*Lambda1')`の`Lambda1`に関する微分）

4. **G5 (相補性条件 (X*L)_sym * Lambda_P' = 0)**:
   - `(X*L)_sym = (X*L + (X*L)')/2`の微分は正しいか？
   - `sol.dL`の計算（0.5の係数）は正しいか？

5. **G7 (等号制約 X*L = (X*L)')**:
   - `vec(X*L - (X*L)')`の計算は正しいか？
   - `sol.dL = kron(speye(n),X) - kron(X,speye(n))*commutation(T,n)`は正しいか？

6. **G8 (ゲージ固定 Lambda3 + Lambda3' = 0)**:
   - `vec(Lambda3 + Lambda3') = (I + K_nn) * vec(Lambda3)`の計算は正しいか？
   - `sol.dLambda3 = speye(n*n) + Knn`は正しいか？
   - この条件により`Lambda3`の対称成分の自由度が除去されるか？

### 5.2 データ依存項の計算

1. **dF_dD.m**:
   - F1の各ブロック（`X*L`, `Z*L`, `U*L`）のデータ依存部分の抽出は正しいか？
   - Commutation行列の使用は正しいか？

2. **dPi_dXU.m**:
   - `Pi = eye(T) - pinv(Gamma)*Gamma`の微分は正しいか？
   - `dpinvGamma_dGamma`の計算は正しいか？

3. **dVec_Lambda_TdFdL_dD.m**:
   - `vec(Lambda' * dF_dL)`のデータ依存部分の抽出は正しいか？
   - Lambdaのブロック分解は正しいか？

### 5.3 次元の整合性

1. **各G_merge関数の出力サイズ**:
   - `G1_row_without_Data`: `1 × (1 + T*n + (2n+2m)^2 + n*n + n*n)`
   - `G2_row_without_Data`: `T*n × (1 + T*n + (2n+2m)^2 + n*n + n*n)`
   - `G3_row_without_Data`: `(2n+2m)^2 × (1 + T*n + (2n+2m)^2 + n*n + n*n)`
   - `G8_row_without_Data`: `n*n × (1 + T*n + (2n+2m)^2 + n*n + n*n)`
   - など、すべての行の列数が一致しているか？

2. **Data項のサイズ**:
   - すべての`Data`項が`T*(2*n+m)`列（データ`D`のサイズ）になっているか？

### 5.4 符号の確認

1. **Lagrangianの符号**:
   - `L = -delta + gamma*||Pi*L||_F^2 - tr(F1*Lambda1') - tr((X*L)_sym*Lambda_P') + tr(Lambda3'*(X*L - (X*L)'))`
   - 各項の符号は正しいか？
   - 等号制約`X*L = (X*L)'`の双対変数`Lambda3`の符号は正しいか？

2. **微分の符号**:
   - `∂L/∂L = 0`の計算で、`-dF_dL' * vec(Lambda)`の符号は正しいか？
   - 相補性条件の微分で、符号は正しいか？

### 5.5 Commutation行列の使用

1. **正しい使用**:
   - `vec(A^T) = K * vec(A)`の関係を正しく使用しているか？
   - Commutation行列のサイズは正しいか？

### 5.6 対称性の扱い

1. **Lambda1の対称性**:
   - `G4`で`Lambda1 - Lambda1' = 0`を正しく実装しているか？

2. **Lambda_Pの対称性**:
   - `G6`で`Lambda_P - Lambda_P' = 0`を正しく実装しているか？

3. **X*Lの対称性**:
   - `G7`で`X*L = (X*L)'`を正しく実装しているか？

4. **Lambda3のゲージ固定**:
   - `G8`で`Lambda3 + Lambda3' = 0`（歪対称）を正しく実装しているか？
   - この条件により、`Lambda3`の対称成分の自由度が除去され、rank deficiencyが解消されるか？

---

## 6. 参考情報

### 6.1 理論的背景

詳細な理論的説明は`doc/regularization_implicit_differentiation.tex`を参照してください。

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
5. **データ依存項の計算**: `dF_dD`, `dPi_dXU`などの計算が正しいか

ご回答をお待ちしております。

