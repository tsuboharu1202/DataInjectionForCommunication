function sol = G5_merge(n,m,T,X,Lambda_P,L)
% G5 は Lambda_P の相補性条件
n1 = 2*n+2*m;
sol = struct();
sol.dL = 0.5*kron(Lambda_P,X) + 0.5*kron(Lambda_P*X,speye(n))*core.helper.commutation(T,n);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = kron(speye(n),X*L)*core.helper.commutation(n,n);

Ex = [sparse(n,n), speye(n),sparse(n,m)];
sol.Data = 0.5*kron(Lambda_P*L',Ex) + 0.5*kron(Lambda_P*Ex,L')*core.helper.commutation(2*n+m,T);
sol.dLambda3 = sparse(n*n, n*n);
sol.G5_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P, sol.dLambda3];
end