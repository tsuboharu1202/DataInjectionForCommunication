function [dX, dZ, dU, checkX, checkZ, checkU] = projector(dX, dZ, dU, checkX, checkZ, checkU, epsilon)
% projector: ノルム制約を適用する投影関数
%
% 引数（全て必須）:
%   dX, dZ, dU: ノイズ
%   checkX, checkZ, checkU: 処理済みフラグ
%   epsilon: ノルム制約の上限
%
% 出力:
%   dX, dZ, dU: 投影後のノイズ
%   checkX, checkZ, checkU: 更新された処理済みフラグ

if nargin < 7
    error('projector:MissingArgs', ...
        '全ての引数（dX, dZ, dU, checkX, checkZ, checkU, epsilon）は必須です。');
end

% X
if ~isempty(dX)
    mX = ~checkX & (abs(dX) > epsilon);    % まだ未処理 かつ |d|>ε
    dX(mX)     = epsilon .* sign(dX(mX));  % ε·sign に張り付け
    checkX(mX) = true;                     % 処理済みにする
end
% Z
if ~isempty(dZ)
    mZ = ~checkZ & (abs(dZ) > epsilon);
    dZ(mZ)     = epsilon .* sign(dZ(mZ));
    checkZ(mZ) = true;
end
% U
if ~isempty(dU)
    mU = ~checkU & (abs(dU) > epsilon);
    dU(mU)     = epsilon .* sign(dU(mU));
    checkU(mU) = true;
end
end
