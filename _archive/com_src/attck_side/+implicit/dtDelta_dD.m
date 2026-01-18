% function dtDelta_dD = dtDelta_dD(n,m,T,B,G,Phi,Lambda1,F1,F2,Lambda3,F3,alpha,Lambda_alpha,beta,Lambda_beta,tDelta,Lambda_tDelta,Lambda_Y,Y)
% % Dの形状は D = [Z',X', U']'

% G1_sol = implicit.G_grad.G1_merge(n,m,B,T);
% G2_sol = implicit.G_grad.G2_merge(n,m,B,T);
% G3_sol = implicit.G_grad.G3_merge(n,m,T,Lambda1,F2,B,G,Phi);
% G4_sol = implicit.G_grad.G4_merge(n,m,T);
% G5_sol = implicit.G_grad.G5_merge(n,m,T,B);
% G6_sol = implicit.G_grad.G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi);
% G7_sol = implicit.G_grad.G7_merge(n,m,T,Lambda3,F3);
% G8_sol = implicit.G_grad.G8_merge(n,m,T,alpha,Lambda_alpha);
% G9_sol = implicit.G_grad.G9_merge(n,m,T,beta,Lambda_beta);
% G10_sol = implicit.G_grad.G10_merge(n,m,T,tDelta,Lambda_tDelta);
% G11_sol = implicit.G_grad.G11_merge(n,m,T);
% G12_sol = implicit.G_grad.G12_merge(n,m,T);
% G13_sol = implicit.G_grad.G13_merge(n,m,T);
% G14_sol = implicit.G_grad.G14_merge(n,m,T,Lambda_Y,Y);
% G15_sol = implicit.G_grad.G15_merge(n,m,T,Lambda_Y);

% Matrix_H = [
%     G5_sol.G5_row_without_Data;
%     G1_sol.G1_row_without_Data;
%     G2_sol.G2_row_without_Data;
%     G3_sol.G3_row_without_Data;
%     G4_sol.G4_row_without_Data;
%     G6_sol.G6_row_without_Data;
%     G7_sol.G7_row_without_Data;
%     G8_sol.G8_row_without_Data;
%     G9_sol.G9_row_without_Data;
%     G10_sol.G10_row_without_Data;
%     G11_sol.G11_row_without_Data;
%     G12_sol.G12_row_without_Data;
%     G13_sol.G13_row_without_Data;
%     G14_sol.G14_row_without_Data;
%     G15_sol.G15_row_without_Data];

% Mtarix_y = [
%     G5_sol.Data;
%     G1_sol.Data;
%     G2_sol.Data;
%     G3_sol.Data;
%     G4_sol.Data;
%     G6_sol.Data;
%     G7_sol.Data;
%     G8_sol.Data;
%     G9_sol.Data;
%     G10_sol.Data;
%     G11_sol.Data;
%     G12_sol.Data;
%     G13_sol.Data;
%     G14_sol.Data;
%     G15_sol.Data];

% % sdpvarが含まれている可能性があるため、value()で数値型に変換
% Matrix_H_full = double(full(Matrix_H));
% Matrix_y_full = double(full(Mtarix_y));

% rank_H = double(rank(Matrix_H_full));
% n_rows = double(size(Matrix_H_full, 1));
% n_cols = double(size(Matrix_H_full, 2));
% fprintf('Matrix_H: size=[%d,%d], rank=%d\n', n_rows, n_cols, rank_H);
% if rank_H == n_cols
%     fprintf('  解は一意です (rank = 列数)\n');
% else
%     fprintf('  解は一意ではありません (自由度: %d = 列数 %d - rank %d)\n', n_cols - rank_H, n_cols, rank_H);

%     % Null space分析
%     fprintf('\n=== Null Space分析 ===\n');
%     s = svd(Matrix_H_full);
%     fprintf('SVD: 最小特異値 = %e, 最大特異値 = %e, 条件数 = %e\n', min(s), max(s), max(s)/min(s));
%     fprintf('下位10個の特異値:\n');
%     fprintf('  %e\n', s(end-9:end));

%     % Null spaceを計算
%     tol_null = 1e-6;
%     nullvec = null(Matrix_H_full, tol_null);
%     if ~isempty(nullvec)
%         fprintf('Null space次元: %d (tol=%e)\n', size(nullvec, 2), tol_null);

%         % 各変数ブロックのサイズ
%         n_L = n*m;
%         n_Y = n*n;
%         n_alpha = 1;
%         n_beta = 1;
%         n_tDelta = 1;
%         n_Lambda1 = (3*n+m)*(3*n+m);
%         n_Lambda3 = (n+m)*(n+m);
%         n_Lambda_alpha = 1;
%         n_Lambda_beta = 1;
%         n_Lambda_tDelta = 1;
%         n_Lambda_Y = n*n;

%         % 各ブロックのインデックス
%         idx_L = 1:n_L;
%         idx_Y = n_L + (1:n_Y);
%         idx_alpha = n_L + n_Y + 1;
%         idx_beta = n_L + n_Y + n_alpha + 1;
%         idx_tDelta = n_L + n_Y + n_alpha + n_beta + 1;
%         idx_Lambda1 = n_L + n_Y + n_alpha + n_beta + n_tDelta + (1:n_Lambda1);
%         idx_Lambda3 = idx_Lambda1(end) + (1:n_Lambda3);
%         idx_Lambda_alpha = idx_Lambda3(end) + 1;
%         idx_Lambda_beta = idx_Lambda_alpha + 1;
%         idx_Lambda_tDelta = idx_Lambda_beta + 1;
%         idx_Lambda_Y = idx_Lambda_tDelta + (1:n_Lambda_Y);

%         % 各null vectorの各ブロックでのノルム
%         fprintf('\n各変数ブロックでのnull space成分のノルム:\n');
%         for i = 1:size(nullvec, 2)
%             fprintf('  Null vector %d:\n', i);
%             fprintf('    dL:        %e\n', norm(nullvec(idx_L, i)));
%             fprintf('    dY:        %e\n', norm(nullvec(idx_Y, i)));
%             fprintf('    dAlpha:    %e\n', norm(nullvec(idx_alpha, i)));
%             fprintf('    dBeta:     %e\n', norm(nullvec(idx_beta, i)));
%             fprintf('    dtDelta:   %e\n', norm(nullvec(idx_tDelta, i)));
%             fprintf('    dLambda1:  %e\n', norm(nullvec(idx_Lambda1, i)));
%             fprintf('    dLambda3:  %e\n', norm(nullvec(idx_Lambda3, i)));
%             fprintf('    dLambda_alpha: %e\n', norm(nullvec(idx_Lambda_alpha, i)));
%             fprintf('    dLambda_beta:  %e\n', norm(nullvec(idx_Lambda_beta, i)));
%             fprintf('    dLambda_tDelta: %e\n', norm(nullvec(idx_Lambda_tDelta, i)));
%             fprintf('    dLambda_Y: %e\n', norm(nullvec(idx_Lambda_Y, i)));
%         end
%     else
%         fprintf('Null spaceは空です (tol=%e)\n', tol_null);
%     end
% end
% matrix_x = Matrix_H_full \ (-Matrix_y_full);

% % 回の1行目の要素がdtDelta_dDの値
% dtDelta_dD = matrix_x(1,:);
% dtDelta_dD = dtDelta_dD';

% % 形状を整える
% dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
% end
