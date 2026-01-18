% function F1 = F1_of_L(n, m, T, B, X, Z, U, L, delta)
% %F1_OF_L F1行列を数値計算する（最適化なし）
% %   F1 = F1_of_L(n, m, T, B, X, Z, U, L, delta)
% %
% %   入力:
% %       n: 状態次元
% %       m: 入力次元
% %       T: 時間ステップ数
% %       B: n×m 入力行列
% %       X: n×T 状態行列
% %       Z: n×T 次状態行列
% %       U: m×T 入力行列
% %       L: T×n 制御パラメータ行列
% %       delta: スカラー（量子化パラメータ）
% %
% %   出力:
% %       F1: (2n+2m)×(2n+2m) 行列
% %
% %   注意:
% %       - solve_sdp内のconst_matと完全に同じ定義
% %       - すべてdouble数値で動作（sdpvar不要）
% %
% %   使用例:
% %       F1 = implicit_regularization.helper.F1_of_L(n, m, T, B, X, Z, U, L_val, delta_val);

% % solve_sdp内の定義と完全に同じ
% Xm_L_Sym = (X*L + (X*L)')/2;
% F1 = [Xm_L_Sym, (Z*L)', zeros(n,m), (U*L)';
%     Z*L,      Xm_L_Sym, delta*B,  zeros(n,m);
%     zeros(m,n), delta*B', eye(m), zeros(m,m);
%     U*L, zeros(m,n), zeros(m,m), eye(m)];

% end

