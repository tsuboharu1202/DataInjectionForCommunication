function sol = G14_merge(n,m,T,Lambda_Y,Y)
% G14 は 相補性条件における Lambda_Y * Y = 0
% 制約: Y >= 0 (数値的安定性のため tolerance*eye を使うが、本質的には >= 0)
% 相補性: Lambda_Y * Y = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = kron(Lambda_Y,eye(n));
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);


sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
% 修正: d/dLambda_Y [Lambda_Y * Y] = Y なので、kron(speye(n), Y) を使う
sol.Lambda_Y = kron(speye(n),Y)*implicit.helper.commutation(n,n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G14_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end