% function [gradX_delta, gradZ_delta, gradU_delta] = calc_grad(sd, opts)
% % calc_grad: Delta勾配を計算（Implicit Differentiation使用、regularizationのみ）
% %   sd: SystemData object
% %   opts: (optional) 追加オプション（gammaなど）

% if nargin < 2
%     opts = struct();
% end

% % システムサイズ
% [n, T] = size(sd.Z);
% m = size(sd.U, 1);
% B = sd.B;

% % SDPを解く（regularizationのみ）
% gamma = 1e3;  % デフォルト値（optsから取得可能にする場合はここで処理）
% if isfield(opts, 'gamma')
%     gamma = opts.gamma;
% end

% [sol, ~, ~, ~, ~] = regularization_sdp.solve_sdp(sd, gamma);

% % Implicit differentiationで勾配計算
% Gamma = [sd.U; sd.X];
% Pi = eye(T) - pinv(Gamma)*Gamma;

% dtDelta_dD = implicit_regularization.dtDelta_dD(n, m, T, B, sd.X, sd.Z, sd.U, ...
%     Pi, gamma, Gamma, sol.L, sol.Lambda1, sol.Lambda2, sol.F1,sol.Lambda3);

% % dtDelta_dDの形状: (2*n+m) × T
% % D = [Z', X', U']' なので、dtDelta_dD = [d(delta)/dZ; d(delta)/dX; d(delta)/dU]
% gradZ_delta = dtDelta_dD(1:n, :);
% gradX_delta = dtDelta_dD(n+1:2*n, :);
% gradU_delta = dtDelta_dD(2*n+1:end, :);

% end
