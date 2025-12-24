function dtDelta_dD = dtDelta_dD(n,m,T,B,G,Phi,Lambda1,F1,F2,Lambda3,F3,alpha,Lambda_alpha,beta,Lambda_beta,tDelta,Lambda_tDelta,Lambda_Y,Y)
% Dの形状は D = [Z',X', U']'

G1_sol = implicit.G_grad.G1_merge(n,m,B,T);
G2_sol = implicit.G_grad.G2_merge(n,m,B,T);
G3_sol = implicit.G_grad.G3_merge(n,m,T,Lambda1,F2,B,G,Phi);
G4_sol = implicit.G_grad.G4_merge(n,m,T);
G5_sol = implicit.G_grad.G5_merge(n,m,T,B);
G6_sol = implicit.G_grad.G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi);
G7_sol = implicit.G_grad.G7_merge(n,m,T,Lambda3,F3);
G8_sol = implicit.G_grad.G8_merge(n,m,T,alpha,Lambda_alpha);
G9_sol = implicit.G_grad.G9_merge(n,m,T,beta,Lambda_beta);
G10_sol = implicit.G_grad.G10_merge(n,m,T,tDelta,Lambda_tDelta);
G11_sol = implicit.G_grad.G11_merge(n,m,T);
G12_sol = implicit.G_grad.G12_merge(n,m,T);
G13_sol = implicit.G_grad.G13_merge(n,m,T);
G14_sol = implicit.G_grad.G14_merge(n,m,T,Lambda_Y,Y);
G15_sol = implicit.G_grad.G15_merge(n,m,T,Lambda_Y);

Matrix_H = [
    G5_sol.G5_row_without_Data;
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    G8_sol.G8_row_without_Data;
    G9_sol.G9_row_without_Data;
    G10_sol.G10_row_without_Data;
    G11_sol.G11_row_without_Data;
    G12_sol.G12_row_without_Data;
    G13_sol.G13_row_without_Data;
    G14_sol.G14_row_without_Data;
    G15_sol.G15_row_without_Data];

Mtarix_y = [
    G5_sol.Data;
    G1_sol.Data;
    G2_sol.Data;
    G3_sol.Data;
    G4_sol.Data;
    G6_sol.Data;
    G7_sol.Data;
    G8_sol.Data;
    G9_sol.Data;
    G10_sol.Data;
    G11_sol.Data;
    G12_sol.Data;
    G13_sol.Data;
    G14_sol.Data;
    G15_sol.Data];

Matrix_H_full = full(Matrix_H);
matrix_x = Matrix_H_full \ (-full(Mtarix_y));

% 回の1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
end
