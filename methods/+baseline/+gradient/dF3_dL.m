function dF3_dL = dF3_dL(n,m)

E1 = [speye(n);
    zeros(m,n)];


E2 = [zeros(n,m);
    speye(m)];

% 誤りの可能性あり
% dF3_dL = kron(E2,E1) + kron(E1,E2)*gradient.commutation(m,n);
dF3_dL = kron(E1, E2) + kron(E2, E1) * gradient.commutation(m,n);