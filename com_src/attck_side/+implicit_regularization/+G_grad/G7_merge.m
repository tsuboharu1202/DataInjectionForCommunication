function sol = G7_merge(n,m,T,X,L)
% G7 は XL = (XL)'
n1 = 2*n+2*m;
sol = struct();
sol.dL = kron(speye(n),X) - kron(X,speye(n))*implicit.helper.commutation(T,n);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = sparse(n*n,n*n);
C2nmT = implicit.helper.commutation(2*n+m,T);
Ex = [sparse(n,n), speye(n),sparse(n,m)];
sol.Data = kron(L',Ex)- kron(Ex,L')*C2nmT;
sol.G7_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P];
end