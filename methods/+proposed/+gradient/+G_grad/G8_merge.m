function sol = G8_merge(n, m, T)
sol = struct();

Knn = core.helper.commutation(n,n);
r = n*n; % rows

sol.dtDelta   = sparse(r,1);
sol.dL        = sparse(r,T*n);

n1 = 2*n+2*m;
sol.dLambda   = sparse(r,n1*n1);
sol.dLambda_P = sparse(r,n*n);

sol.dLambda3  = speye(r) + Knn;

sol.Data      = sparse(r, T*(2*n+m));

sol.G8_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
end
