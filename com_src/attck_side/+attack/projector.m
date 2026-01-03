function [dX, dZ, dU, checkX, checkZ, checkU] = projector(dX, dZ, dU, checkX, checkZ, checkU, epsilon)
% projector: ノルム制約を適用する投影関数
%
% 入力:
%   dX, dZ, dU: ノイズ
%   checkX, checkZ, checkU: 処理済みフラグ
%   epsilon: (オプション) ノルム制約の上限。指定しない場合はcfg.Const.ATTACKER_UPPERLIMITを使用
%
% 出力:
%   dX, dZ, dU: 投影後のノイズ
%   checkX, checkZ, checkU: 更新された処理済みフラグ

if nargin < 7 || isempty(epsilon)
    epsilon = cfg.Const.ATTACKER_UPPERLIMIT;
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

