function dtDelta_dD = dtDelta_dD(n,m,T,B,G,Phi,Lambda1,F2,Lambda3,F3,alpha,Lambda_alpha,beta,Lambda_beta,tDelta,Lambda_tDelta)
% Dの形状は D = [Z',X', U']'
G1_sol = implicit.G_grad.G1_merge(n,m,B,T);
G2_sol = implicit.G_grad.G2_merge(n,m,B,T);
G3_sol = implicit.G_grad.G3_merge(n,m,T,Lambda1,F2,B,G,Phi);
G4_sol = implicit.G_grad.G4_merge(n,m,T);
G5_sol = implicit.G_grad.G5_merge(n,m,T,B);
G6_sol = implicit.G_grad.G6_merge(n,m,T,Lambda1,alpha,F2,B,G,Phi);
G7_sol = implicit.G_grad.G7_merge(n,m,T,Lambda3,F3);
G8_sol = implicit.G_grad.G8_merge(n,m,T,alpha,Lambda_alpha);
G9_sol = implicit.G_grad.G9_merge(n,m,T,beta,Lambda_beta);
G10_sol = implicit.G_grad.G10_merge(n,m,T,tDelta,Lambda_tDelta);
G11_sol = implicit.G_grad.G11_merge(n,m,T);
G12_sol = implicit.G_grad.G12_merge(n,m,T);
G13_sol = implicit.G_grad.G13_merge(n,m,T);

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
    G13_sol.G13_row_without_Data];

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
    G13_sol.Data];

% デバッグ用：Matrix_Hのrank確認
Matrix_H_full = full(Matrix_H);
[r_H, c_H] = size(Matrix_H_full);
rank_H = rank(Matrix_H_full);
fprintf('Matrix_H size: %dx%d, rank: %d\n', r_H, c_H, rank_H);

% 特異性チェック
if rank_H < min(r_H, c_H)
    warning('Matrix_H is rank deficient! rank=%d < min(%d,%d)', rank_H, r_H, c_H);
    fprintf('  Matrix_H condition number: %e\n', cond(Matrix_H_full));
    
    % 解の自由度を確認
    if rank_H < c_H
        fprintf('  解の自由度: %d (列数 %d - rank %d)\n', c_H - rank_H, c_H, rank_H);
        fprintf('  → 解が一意に決まらない可能性があります\n');
    elseif rank_H < r_H
        fprintf('  行の自由度: %d (行数 %d - rank %d)\n', r_H - rank_H, r_H, rank_H);
        fprintf('  → 重複した方程式がありますが、解は一意に決まる可能性があります\n');
    end
end

matrix_x = Matrix_H \ (-Mtarix_y);
% 解の一意性チェック（rank不足の場合）
if rank_H < c_H
    % 解が一意でない場合、最小ノルム解が返される
    % 残差をチェックして、解が妥当か確認
    residual = Matrix_H_full * matrix_x - (-full(Mtarix_y));
    fprintf('  解の残差ノルム: %e (小さいほど良い)\n', norm(residual(:)));
    if norm(residual(:)) > 1e-6
        warning('  残差が大きいです。解が正確でない可能性があります。');
    end
end

% 回の1行目の要素がdtDelta_dDの値
dtDelta_dD = matrix_x(1,:);
dtDelta_dD = dtDelta_dD';

% 形状を整える
dtDelta_dD = reshape(dtDelta_dD, [2*n+m,T]);
