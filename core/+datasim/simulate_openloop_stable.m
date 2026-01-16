function [X,Z,U] = simulate_openloop_stable(A,B,V,x0,opts)
% simulate_openloop_stable: 安定化済みシステムでのオープンループシミュレーション
%
% 出力:
%   X: n×T 状態軌道
%   Z: n×T 1ステップ先の状態（Z(:,t) = x(t+1)）
%   U: m×T 実際に適用された入力
%
% 引数:
%   A, B: システム行列
%   V: m×T 探索信号（入力の摂動成分）
%   x0: (optional) 初期状態（デフォルト: randn(n,1)）
%   opts: (optional) オプション構造体
%     - noise_std: プロセスノイズの標準偏差（デフォルト: 0 = ノイズなし）
%
% 論文の「安定なサンプリング」方針:
%   事前に安定化フィードバック u = -K0 x + v を入れて
%   x_{t+1} = (A - B K0) x_t + B v_t + w_t
%   として、探索信号 v_t (=V) はPEだが状態は暴れないようにする。
%
% 注意: MATLABの dlqr は u = -Kx を返す（標準）ので、その約束に合わせる。

[n, T] = deal(size(A,1), size(V,2));
m = size(B,2);

% オプションのデフォルト値
if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'noise_std')
    opts.noise_std = 0;  % デフォルト: ノイズなし
end

% 事前安定化ゲイン K0 を設計
K0 = zeros(m,n);
if max(abs(eig(A))) < 0.99
    Ac = A;
else
    % dlqr を使って安定化ゲインを設計
    Q0 = eye(n);
    R0 = eye(m);
    K0 = dlqr(A,B,Q0,R0);  % u = -K0 x
    Ac = A - B*K0;
end

% 初期値
if nargin < 4 || isempty(x0)
    x0 = randn(n,1);
end

% シミュレーション実行
X = zeros(n,T); Z = zeros(n,T);
x = x0;
for t = 1:T
    X(:,t) = x;
    % ノイズを追加（opts.noise_std > 0 の場合のみ）
    if opts.noise_std > 0
        noise = opts.noise_std * randn(n,1);
    else
        noise = zeros(n,1);
    end
    x = Ac*x + B*V(:,t) + noise;
    Z(:,t) = x;
end

% 実際に適用された入力: u_t = -K0 x_t + v_t
U = -K0*X + V;
end
