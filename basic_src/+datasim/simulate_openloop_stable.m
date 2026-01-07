function [X,Z,U] = simulate_openloop_stable(A,B,V,x0)
% X: n×T,  Z: n×T  (Zは1ステップ先)
% 論文の「安定なサンプリング」方針:
% 事前に安定化フィードバック u = -K0 x + v を入れて
%   x_{t+1} = (A - B K0) x_t + B v_t
% として、探索信号 v_t (=V) はPEだが状態は暴れないようにする。
%
% 注意: MATLABの dlqr は u = -Kx を返す（標準）ので、その約束に合わせる。

if nargin < 5 || isempty(w_max), w_max = 1e-2; end
[n, T] = deal(size(A,1), size(V,2));
m = size(B,2);

% 事前安定化ゲイン K0 を設計
K0 = zeros(m,n);
if max(abs(eig(A))) < 0.99
    Ac = A;
else
    % 1) まずは dlqr/dare を試す（モデル既知でデータ生成する段階なのでOK）
    Q0 = eye(n);
    R0 = eye(m);
    K0 = dlqr(A,B,Q0,R0);      % u = -K0 x
    Ac = A - B*K0;
end

% 初期値
if nargin<4 || isempty(x0)
    x0 = randn(n,1);
end

% シミュレーション実行
X = zeros(n,T); Z = zeros(n,T);
x = x0;
for t = 1:T
    X(:,t) = x;
    x = Ac*x + B*V(:,t) + w_max*randn(n,1);
    % x = Ac*x + B*V(:,t);
    Z(:,t) = x;
end

% Input actually applied to the original system: u_t = -K0 x_t + v_t
U = -K0*X + V;
end
