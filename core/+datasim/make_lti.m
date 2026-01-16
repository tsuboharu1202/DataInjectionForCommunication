function [A,B] = make_lti(n,m)
% デモ用の安定系を作る。実系があるなら読み込みに差し替えてOK
% ランダム生成 + 安定化

% ランダムなシステムを生成
% 不安定なシステムを作る
A = (rand(n,n)-0.5*ones(n,n))*2 + eye(n);
B = (rand(n,m)-0.5*ones(n,m))*2;

% システムを安定化（すべての固有値の実部を負にする）
% A_eig = eig(A);
% max_real = max(real(A_eig));
% if max_real >= 0
%     % 不安定な場合、シフトして安定化
%     A = A - (max_real + 0.1) * eye(size(A));
% end

% % 安定性を確認
% A_eig_final = eig(A);
% assert(all(real(A_eig_final) < 0), 'System A is not stable after stabilization');

end
