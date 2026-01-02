function dtDelta_dD = dtDelta_dD(n,m,T,B,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,F1)
% Dの形状は D = [Z',X', U']'

G1_sol = implicit_regularization.G_grad.G1_merge(n,m,B,T);
G2_sol = implicit_regularization.G_grad.G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P);
G3_sol = implicit_regularization.G_grad.G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B);
G4_sol = implicit_regularization.G_grad.G4_merge(n,m,T);
G5_sol = implicit_regularization.G_grad.G5_merge(n,m,T,X,Lambda_P,L);
G6_sol = implicit_regularization.G_grad.G6_merge(n,m,T);
G7_sol = implicit_regularization.G_grad.G7_merge(n,m,T,X,L);
fprintf('G1_sol.G1_row_without_Data: %dx%d\n', size(G1_sol.G1_row_without_Data,1), size(G1_sol.G1_row_without_Data,2));
fprintf('G2_sol.G2_row_without_Data: %dx%d\n', size(G2_sol.G2_row_without_Data,1), size(G2_sol.G2_row_without_Data,2));
fprintf('G3_sol.G3_row_without_Data: %dx%d\n', size(G3_sol.G3_row_without_Data,1), size(G3_sol.G3_row_without_Data,2));
fprintf('G4_sol.G4_row_without_Data: %dx%d\n', size(G4_sol.G4_row_without_Data,1), size(G4_sol.G4_row_without_Data,2));

Matrix_H = [
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G5_sol.G5_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    ];

Mtarix_y = [
    G1_sol.Data;
    G2_sol.Data;
    G3_sol.Data;
    G4_sol.Data;
    G5_sol.Data;
    G6_sol.Data;
    G7_sol.Data;
    ];

% sdpvarが含まれている可能性があるため、value()で数値型に変換
Matrix_H_full = double(full(Matrix_H));
Matrix_y_full = double(full(Mtarix_y));

rank_H = double(rank(Matrix_H_full));
n_rows = double(size(Matrix_H_full, 1));
n_cols = double(size(Matrix_H_full, 2));
fprintf('Matrix_H: size=[%d,%d], rank=%d\n', n_rows, n_cols, rank_H);
if rank_H == n_cols
    fprintf('  解は一意です (rank = 列数)\n');
else
    fprintf('  解は一意ではありません (自由度: %d = 列数 %d - rank %d)\n', n_cols - rank_H, n_cols, rank_H);
    
end
matrix_x = Matrix_H_full \ (-Matrix_y_full);

% 回の1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
end
