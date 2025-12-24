function dF1_dY = dF1_dY(n,m)
E1 = [speye(n);
    zeros(2*n+m,n)];

E2 = [zeros(n,n);
    speye(n);
    zeros(n+m,n)];

E3 = [zeros(n,n);
    zeros(n,n);
    speye(n);
    zeros(m,n)];
dF1_dY = kron(E1,E1) + kron(E2,E2) + kron(E3,E2) + kron(E2,E3)*implicit.helper.commutation(n,n);