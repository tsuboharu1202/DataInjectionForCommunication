function U = make_inputU(m, T)
% make_inputU: ランダムウォークで入力信号を生成
%
% 引数:
%   m: 入力次元（必須）
%   T: サンプル数（必須）
%
% 出力:
%   U: m × T の入力信号行列

if nargin < 2
    error('make_inputU:MissingArgs', ...
        'm（入力次元）と T（サンプル数）は必須です。');
end

U = zeros(m, T);
U(:, 1) = randn(m, 1);
for t = 2:T
    st = 0.3;
    U(:, t) = U(:, t-1) + st * randn(m, 1);
end
end
