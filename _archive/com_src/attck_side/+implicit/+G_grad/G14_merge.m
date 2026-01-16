function sol = G14_merge(n,m,T,Lambda_Y,Y)
% G14 は 相補性条件における Y*Lambda_Y' = 0
% 統一形式: (制約) * Λ^T = 0
% 制約: Y >= 0 (数値的安定性のため tolerance*eye を使うが、本質的には >= 0)
% 相補性: Y * Lambda_Y' = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);

% d/dY [vec(Y*Lambda_Y')] = d/dY [vec(Y)*Lambda_Y']
% vec(Y*Lambda_Y') = (Lambda_Y ⊗ I) vec(Y)
% したがって、d/dY [vec(Y*Lambda_Y')] = Lambda_Y ⊗ I
sol.dY = kron(Lambda_Y,eye(n));
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);

sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
% d/dLambda_Y [vec(Y*Lambda_Y')]
% vec(Y*Lambda_Y') = (I ⊗ Y) vec(Lambda_Y')
% vec(Lambda_Y') = K vec(Lambda_Y)
% したがって、d/dLambda_Y [vec(Y*Lambda_Y')] = (I ⊗ Y) K
% 実際のコードでは K (I ⊗ Y) の形式を使用
sol.Lambda_Y = kron(speye(n),Y)*implicit.helper.commutation(n,n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G14_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end