% function dpinvGamma_dGamma = dpinvGamma_dGamma(n,m,T,Gamma)

% Cnm_T = core.helper.commutation(n+m,T);
% GGT = Gamma*Gamma';
% GGT_inv = GGT^(-1);

% term1 = kron(GGT_inv',speye(T))*Cnm_T;
% term2 = -kron(speye(n+m),Gamma')*kron(GGT_inv',GGT_inv)*(kron(Gamma,speye(m+n)) + kron(speye(m+n),Gamma))*Cnm_T;
% dpinvGamma_dGamma = term1 + term2;
% end

function dpinvGamma_dGamma = dpinvGamma_dGamma(n,m,T,Gamma)
p = n + m;    % rows of Gamma
q = T;        % cols of Gamma

Kpq = core.helper.commutation(p,q);   % vec(Gamma') = Kpq * vec(Gamma)

G   = Gamma*Gamma';                        % p x p
Ginv = G \ speye(p);                       % safer than G^(-1)

% term1: vec(dGamma' * Ginv) = (Ginv' ⊗ I_q) vec(dGamma')
term1 = kron(Ginv', speye(q)) * Kpq;

% vec(dG) = (Gamma ⊗ I_p) vec(dGamma) + (I_p ⊗ Gamma) vec(dGamma')
%         = [(Gamma ⊗ I_p) + (I_p ⊗ Gamma)Kpq] vec(dGamma)
dG_dGamma = kron(Gamma, speye(p)) + kron(speye(p), Gamma) * Kpq;   % (p^2 x pq)

% term2: -vec(Gamma' Ginv (dG) Ginv)
% vec(A (dG) B) = (B' ⊗ A) vec(dG)
A = Gamma' * Ginv;                  % q x p
term2 = -kron(Ginv', A) * dG_dGamma; % (pq x pq)

dpinvGamma_dGamma = term1 + term2;
end
