function dF_dL = dF_dL(n,m,T,X,Z,U)
E1 = [speye(n);
    sparse(n+2*m,n)];
E2 = [sparse(n,n);
    speye(n);
    sparse(2*m,n)];
X1 = [X;
    sparse(2*m+n,T)];
X2 = [sparse(n,T);
    X;
    sparse(2*m,T)];
Z2 = [sparse(n,T);
    Z;
    sparse(2*m,T)];
U4 = [sparse(2*n+m,T);
    U];
dF_dL = kron(E1,X1)+ kron(E2,X2) + kron(E2,Z2) +...
    kron(E1,Z2) +kron(Z2,E1)*implicit.helper.commutation(T,n) +...
    kron(E1,U4) +kron(U4,E1)*implicit.helper.commutation(T,n);

end