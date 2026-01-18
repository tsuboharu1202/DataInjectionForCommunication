% function dF2_dD = dF2_dD(n,m,B,T,G,Phi)
% Ez = [speye(n), zeros(n,n), zeros(n,m)];
% Ex = [zeros(n,n), speye(n), zeros(n,m)];
% Eu = [zeros(m,n), zeros(m,n), speye(m)];

% n1 = 3*n+m;
% n2 = n+T;

% Eleft = [sparse(n,T);
%     speye(T)];

% Eright = [Ez - B*Eu;
%     -Ex;
%     sparse(n+m,2*n+m)];

% GPhi = G*Phi;


% dG_dD = kron(Eleft,Eright);

% dF2_dD = kron(GPhi,speye(n1))*dG_dD + kron(speye(n1), GPhi)*implicit.helper.commutation(n1,n2)*dG_dD;
% end