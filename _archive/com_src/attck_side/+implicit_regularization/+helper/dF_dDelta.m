% function v = dF_dDelta(n,m,B)
% n1 = 2*n+2*m;

% dF = sparse(n1,n1);

% r2 = (n+1):(2*n);
% c2 = (n+1):(2*n);
% r3 = (2*n+1):(2*n+m);
% c3 = (2*n+1):(2*n+m);

% % δB block: (row2, col3)
% dF(r2, c3) = B;

% % δB^T block: (row3, col2)
% dF(r3, c2) = B';

% v = dF(:); % vec(∂F/∂δ)
% end
