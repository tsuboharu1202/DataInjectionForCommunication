function hinf_norm = hinfnorm_AK(A, B, K)
%HINFNORM_AK H∞ノルムを計算する
%   hinf_norm = hinfnorm_AK(A, B, K) は、離散時間システム
%   G_{A,K}(z) = K(zI - A - BK)^{-1}B
%   のH∞ノルムを計算します。
%
%   入力:
%       A: n×n システム行列
%       B: n×m 入力行列
%       K: m×n フィードバックゲイン行列
%
%   出力:
%       hinf_norm: H∞ノルム（スカラー）

% サイズ確認
[n, ~] = size(A);
[~, m] = size(B);
[m_k, ~] = size(K);

assert(size(A, 2) == n, 'A must be square');
assert(size(B, 1) == n, 'B rows must match A');
assert(size(K, 1) == m, 'K rows must match input dimension');
assert(size(K, 2) == n, 'K cols must match state dimension');

% 閉ループシステム行列
A_cl = A + B * K;

% 離散時間状態空間モデルを作成
% G(z) = K(zI - A_cl)^{-1}B
% 状態空間表現: (A_cl, B, K, 0)
sys = ss(A_cl, B, K, zeros(m_k, m), -1);  % -1は離散時間を意味
hinf_norm = hinfnorm(sys);
end
