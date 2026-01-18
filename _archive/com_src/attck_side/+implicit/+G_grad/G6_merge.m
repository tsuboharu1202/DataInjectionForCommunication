% function sol = G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi)
% % G6 は 相補性条件における (F1-αF2)*Lambda1' = 0
% % 統一形式: (制約) * Λ^T = 0
% n1_2 = (3*n+m)*(3*n+m);
% n2_2 = (n+m)*(n+m);
% eye = speye(3*n+m);
% F1_minus_alphaF2 = F1 - alpha*F2;

% sol = struct();

% % d/dL [vec((F1-alpha*F2)*Lambda1')] = d/dL [vec(F1-alpha*F2)*Lambda1']
% % vec((F1-alpha*F2)*Lambda1') = (Lambda1 ⊗ I) vec(F1-alpha*F2)
% % したがって、d/dL [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) dvec(F1)/dL
% sub_matrix = kron(Lambda1,eye);
% sol.dL = sub_matrix*implicit.helper.dF1_dL(n,m,B);
% sol.dY = sub_matrix*implicit.helper.dF1_dY(n,m);

% % d/dalpha [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) d/dalpha [vec(F1-alpha*F2)]
% % d/dalpha [F1-alpha*F2] = -F2
% % したがって、d/dalpha [vec((F1-alpha*F2)*Lambda1')] = -(Lambda1 ⊗ I) vec(F2) = -vec(Lambda1*F2)
% Lambda1F2 = F2*Lambda1';
% sol.dAlpha = -Lambda1F2(:);

% E11 = [speye(n);
%     zeros(2*n+m,n)];
% M11 = E11*E11';
% % d/dbeta [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) d/dbeta [vec(F1)]
% % dF1/dbeta = -M11
% % したがって、d/dbeta [vec((F1-alpha*F2)*Lambda1')] = -vec(Lambda1*M11)
% Lambda1M11 = M11*Lambda1';
% sol.dBeta = -Lambda1M11(:);

% E11_withB = [B;
%     zeros(2*n+m,m)];
% M11_withB = E11_withB*E11_withB';
% % d/dtDelta [vec((F1-alpha*F2)*Lambda1')] = (Lambda1 ⊗ I) d/dtDelta [vec(F1)]
% % dF1/dtDelta = -M11_withB
% % したがって、d/dtDelta [vec((F1-alpha*F2)*Lambda1')] = -vec(Lambda1*M11_withB)
% Lambda1M11_withB = M11_withB*Lambda1';
% sol.dtDelta = -Lambda1M11_withB(:);

% % d/dLambda1 [vec((F1-alpha*F2)*Lambda1')]
% % vec((F1-alpha*F2)*Lambda1') = (I ⊗ (F1-alpha*F2)) vec(Lambda1')
% % vec(Lambda1') = K vec(Lambda1)
% % したがって、d/dLambda1 [vec((F1-alpha*F2)*Lambda1')] = (I ⊗ (F1-alpha*F2)) K
% % 実際のコードでは K (I ⊗ (F1-alpha*F2)) の形式を使用
% sol.dLambda1 = kron(eye, F1_minus_alphaF2)*implicit.helper.commutation(3*n+m,3*n+m);
% sol.dLambda3 = sparse(n1_2,n2_2);
% sol.Lambda_Alpha = sparse(n1_2,1);
% sol.Lambda_Beta = sparse(n1_2,1);
% sol.Lambda_tDelta = sparse(n1_2,1);
% sol.Lambda_Y = sparse(n1_2,n*n);
% sol.Data = -alpha*sub_matrix*implicit.helper.dF2_dD(n,m,B,T,G,Phi);

% sol.G6_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
% end