function [gradX, gradZ, gradU] = calc_grad_implicit(sd, minimize_delta)
% calc_grad_implicit: Calculate gradient using implicit differentiation
%
% Inputs:
%   sd: SystemData object
%   minimize_delta: (optional) If true, return -dtDelta/dD (to minimize Delta).
%                   If false (default), return dtDelta/dD (to maximize Delta).
%
% Outputs:
%   gradX, gradZ, gradU: Gradients with respect to X, Z, U
%                        If minimize_delta=true, these are negated to minimize Delta.

if nargin < 2
    minimize_delta = false;
end

% SDPを解いて解とラグランジュ乗数を取得
[sol, K, Y_val, L_val, diagnostics] = original_thesis.solve_sdp(sd);

% システムサイズ
n = size(sd.A, 1);
m = size(sd.B, 2);
T = cfg.Const.SAMPLE_COUNT;

% G, Phiを構築
G = [eye(n), sd.Z - sd.B*sd.U;
    zeros(n,n), -sd.X;
    zeros(n,n), zeros(n,T);
    zeros(m,n+T)];

Phi = [sd.Phi11, sd.Phi12; sd.Phi12', sd.Phi22];

% F1, F2, F3を構築
[F1, F2, F3] = original_thesis.build_lmi_blocks(Y_val, L_val, sol.alpha, sol.beta, sol.tDelta, G, Phi, sd);

% ラグランジュ乗数を取得
Lambda1 = sol.Lambda1;
Lambda3 = sol.Lambda3;
Lambda_alpha = sol.Lambda_alpha;
Lambda_beta = sol.Lambda_beta;
Lambda_tDelta = sol.Lambda_tDelta;
Lambda_Y = sol.Lambda_Y;

% Implicit differentiationで勾配を計算
dtDelta_dD = implicit.dtDelta_dD(n, m, T, sd.B, G, Phi, ...
    Lambda1, F1, F2, Lambda3, F3, ...
    sol.alpha, Lambda_alpha, sol.beta, Lambda_beta, sol.tDelta, Lambda_tDelta, Lambda_Y, Y_val);

% dtDelta_dDは (2*n+m) x T の形状
% D = [Z', X', U']' なので、dtDelta_dDの行は [d/dZ, d/dX, d/dU] に対応
gradZ = dtDelta_dD(1:n, :);
gradX = dtDelta_dD(n+1:2*n, :);
gradU = dtDelta_dD(2*n+1:2*n+m, :);

% Deltaを小さくする方向の勾配が必要な場合、符号を反転
if minimize_delta
    gradZ = -gradZ;
    gradX = -gradX;
    gradU = -gradU;
end

end

