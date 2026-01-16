function [gradX_delta, gradZ_delta, gradU_delta] = calc_grad(sd, opts)
% calc_grad: Delta勾配を計算（Implicit Differentiation使用）
%
% 引数:
%   sd: SystemData object
%   opts: オプション構造体
%     - gamma: 正則化パラメータ（必須）
%
% 戻り値:
%   gradX_delta, gradZ_delta, gradU_delta: 各データに対するdeltaの勾配

if nargin < 2
    opts = struct();
end

% gamma は必須
if ~isfield(opts, 'gamma')
    error('calc_grad:MissingGamma', ...
        'opts.gamma は必須です。正則化パラメータを指定してください。');
end
gamma = opts.gamma;

% システムサイズ
[n, T] = size(sd.Z);
m = size(sd.U, 1);
B = sd.B;

% SDPを解く
[sol, ~, ~, ~, ~] = proposed.solve_sdp(sd, gamma);

% Implicit differentiationで勾配計算
Gamma = [sd.U; sd.X];
Pi = eye(T) - pinv(Gamma)*Gamma;

dtDelta_dD = implicit_regularization.dtDelta_dD(n, m, T, B, sd.X, sd.Z, sd.U, ...
    Pi, gamma, Gamma, sol.L, sol.Lambda1, sol.Lambda2, sol.F1, sol.Lambda3);

% dtDelta_dDの形状: (2*n+m) × T
% D = [Z', X', U']' なので、dtDelta_dD = [d(delta)/dZ; d(delta)/dX; d(delta)/dU]
gradZ_delta = dtDelta_dD(1:n, :);
gradX_delta = dtDelta_dD(n+1:2*n, :);
gradU_delta = dtDelta_dD(2*n+1:end, :);

end
