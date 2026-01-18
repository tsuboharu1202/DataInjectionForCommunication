% function sol = G7_merge(n,m,T,Lambda3,F3)
% % G7 は 相補性条件における F3*Lambda3' = 0
% % 統一形式: (制約) * Λ^T = 0
% n2_2 = (n+m)*(n+m);
% eye = speye(n+m);

% sol = struct();

% % d/dL [vec(F3*Lambda3')] = d/dL [vec(F3)*Lambda3']
% % vec(F3*Lambda3') = (Lambda3 ⊗ I) vec(F3)
% % したがって、d/dL [vec(F3*Lambda3')] = (Lambda3 ⊗ I) dvec(F3)/dL
% sub_matrix = kron(Lambda3,eye);
% sol.dL = sub_matrix*implicit.helper.dF3_dL(n,m);
% sol.dY = sub_matrix*implicit.helper.dF3_dY(n,m);
% sol.dAlpha = sparse(n2_2,1);
% sol.dBeta = sparse(n2_2,1);
% sol.dtDelta = sparse(n2_2,1);

% sol.dLambda1 = sparse(n2_2,(3*n+m)*(3*n+m));
% % d/dLambda3 [vec(F3*Lambda3')]
% % vec(F3*Lambda3') = (I ⊗ F3) vec(Lambda3')
% % vec(Lambda3') = K vec(Lambda3)
% % したがって、d/dLambda3 [vec(F3*Lambda3')] = (I ⊗ F3) K
% % 実際のコードでは K (I ⊗ F3) の形式を使用
% sol.dLambda3 = kron(eye,F3)*implicit.helper.commutation(n+m,n+m);
% sol.Lambda_Alpha = sparse(n2_2,1);
% sol.Lambda_Beta = sparse(n2_2,1);
% sol.Lambda_tDelta = sparse(n2_2,1);
% sol.Lambda_Y = sparse(n2_2,n*n);
% sol.Data = sparse(n2_2,T*(2*n+m));

% sol.G7_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
% end