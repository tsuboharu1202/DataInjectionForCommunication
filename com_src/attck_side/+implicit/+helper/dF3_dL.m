function dF3_dL = dF3_dL(n,m)


E1 = [speye(n);
    zeros(m,n)];


E2 = [zeros(n,m);
    speye(m)];


dF3_dL = kron(E2,E1) + kron(E1,E2)*implicit.helper.commutation(m,n);