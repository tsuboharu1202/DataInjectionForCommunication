function F1 = F1_of_D(n,m,T,B,L,delta,D)
% D = [Z; X; U] in R^{(2n+m) x T}

Ez = [speye(n), sparse(n,n+m)];
Ex = [sparse(n,n), speye(n), sparse(n,m)];
Eu = [sparse(m,2*n), speye(m)];

Z = Ez*D;   % n x T
X = Ex*D;   % n x T
U = Eu*D;   % m x T

Xm_L_Sym = (X*L + (X*L)')/2;
F1 = [Xm_L_Sym, (Z*L)', zeros(n,m), (U*L)';
    Z*L,      Xm_L_Sym, delta*B,  zeros(n,m);
    zeros(m,n), delta*B', eye(m), zeros(m,m);
    U*L, zeros(m,n), zeros(m,m), eye(m)];
end
