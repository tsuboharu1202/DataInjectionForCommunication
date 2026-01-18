function dtDelta_dD = dtDelta_dD(n,m,T,B,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,F1,Lambda3)
% Dの形状は D = [Z',X', U']'

G1_sol = proposed.gradient.G_grad.G1_merge(n,m,B,T);
G2_sol = proposed.gradient.G_grad.G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,Lambda3);
G3_sol = proposed.gradient.G_grad.G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B);
G4_sol = proposed.gradient.G_grad.G4_merge(n,m,T);
G5_sol = proposed.gradient.G_grad.G5_merge(n,m,T,X,Lambda_P,L);
G6_sol = proposed.gradient.G_grad.G6_merge(n,m,T);
G7_sol = proposed.gradient.G_grad.G7_merge(n,m,T,X,L);
G8_sol = proposed.gradient.G_grad.G8_merge(n,m,T);

Matrix_H = [
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G5_sol.G5_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    G8_sol.G8_row_without_Data;
    ];

Mtarix_y = [
    G1_sol.Data;
    G2_sol.Data;
    G3_sol.Data;
    G4_sol.Data;
    G5_sol.Data;
    G6_sol.Data;
    G7_sol.Data;
    G8_sol.Data;
    ];
% Matrix_H_full = double(full(Matrix_H));
% [U,S,V] = svd(Matrix_H_full,'econ');
% sv = diag(S);
% disp('sv(end-20:end):');
% disp(sv(end-20:end))   % 小さい特異値を見る
% vnull = V(:,end);      % 一番小さい方向
% disp('vnull:');
% disp(vnull);




% vnull を各ブロックに分解して、どこが自由度かを見る



% sdpvarが含まれている可能性があるため、value()で数値型に変換
Matrix_H_full = double(full(Matrix_H));
Matrix_y_full = double(full(Mtarix_y));

rank_H = double(rank(Matrix_H_full));
n_rows = double(size(Matrix_H_full, 1));
n_cols = double(size(Matrix_H_full, 2));
% if rank_H ~= n_cols
%     fprintf('  解は一意ではありません (自由度: %d = 列数 %d - rank %d)\n', n_cols - rank_H, n_cols, rank_H);

% end
% matrix_x = Matrix_H_full \ (-Matrix_y_full);
matrix_x = lsqminnorm(Matrix_H_full, -Matrix_y_full);

% 回の1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
end
