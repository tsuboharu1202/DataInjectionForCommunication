function sol = G7_merge(n,m,T,Lambda3,F3)
% G7 は 相補性条件における F3Λ3 = 0
n2_2 = (n+m)*(n+m);
eye = speye(n+m);
sub_matrix = kron(Lambda3,eye);

sol = struct();
sol.dL = sub_matrix*implicit.helper.dF3_dL(n,m);
sol.dY = sub_matrix*implicit.helper.dF3_dY(n,m);
sol.dAlpha = sparse(n2_2,1);
sol.dBeta = sparse(n2_2,1);
sol.dtDelta = sparse(n2_2,1);

sol.dLambda1 = sparse(n2_2,(3*n+m)*(3*n+m));
sol.dLambda3 = kron(eye,F3)*implicit.helper.commutation(n+m,n+m);
sol.Lambda_Alpha = sparse(n2_2,1);
sol.Lambda_Beta = sparse(n2_2,1);
sol.Lambda_tDelta = sparse(n2_2,1);
sol.Data = sparse(n2_2,T*(2*n+m));

sol.G7_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta];
end