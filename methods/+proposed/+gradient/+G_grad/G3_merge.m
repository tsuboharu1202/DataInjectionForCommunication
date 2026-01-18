function sol = G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B)
% G3 は 相補性条件 をdeltaで変微分したもの
n1 = 2*n+2*m;
sol = struct();
Left = kron(Lambda',speye(n1));
dF_dL = proposed.gradient.helper.dF_dL(n,m,T,X,Z,U);
sol.dL = Left*dF_dL;

sol.dLambda = kron(speye(n1),F1)*core.helper.commutation(n1,n1);

% dF/dDelta（スカラーなので vec(∂F/∂δ) を作る）
dF_dDelta = proposed.gradient.helper.dF_dDelta(n,m,B); % (n1^2 x 1)
sol.dtDelta = Left * dF_dDelta;  % ★ここが本質


sol.dLambda_P = sparse(n1*n1,n*n);
sol.Data = Left*proposed.gradient.helper.dF_dD(n,m,T,L);
sol.dLambda3 = sparse(n1*n1, n*n);
sol.G3_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P, sol.dLambda3];

end