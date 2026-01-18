% function Pi = Pi_of_D(n,m,T,D)
% % Pi_of_D : Pi(D) = I_T - pinv(Gamma)*Gamma
% % D is (2n+m) x T, ordered as [Z; X; U]
% %
% % Gamma = [U; X]  ( (m+n) x T )

% Z = D(1:n, :); %#ok<NASGU>
% X = D(n+1:2*n, :);
% U = D(2*n+1:2*n+m, :);

% Gamma = [U; X];  % (m+n) x T

% Pi = eye(T) - pinv(Gamma) * Gamma; % (T x T)
% end
