function sol = G5_merge(n,m,X,Lambda_P,L)
% G5 は Lambda_P の相補性条件
n1 = 2*n+2*m;
sol = struct();
sol.dL = kron(Lambda_P,X);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = kron(speye(n),X*L)*implicit.helper.commutation(n,n);

Ex = [sparse(n,n), speye(n),sparse(n,m)];
sol.Data = kron(Lambda_P*L',Ex);
sol.G5_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P];
end