function test_dF_dD_FD(n,m,T,B,L,delta,D)

fprintf('\n=== dF_dD 有限差分検証 ===\n\n');

F1_base = helper.F1_of_D(n,m,T,B,L,delta,D);
vecF1_base = F1_base(:);

J = helper.dF_dD(n,m,T,L); % あなたの dF_dD
fprintf('dF_dD size = %dx%d (期待: %dx%d)\n', size(J,1), size(J,2), ...
    (2*n+2*m)^2, (2*n+m)*T);

nTrials = 10;
eps_list = [1e-6, 1e-7, 1e-8];

for trial = 1:nTrials
    E = randn(size(D));           % (2n+m) x T
    vecE = E(:);
    
    fprintf('Trial %d/%d:\n', trial, nTrials);
    for e = 1:numel(eps_list)
        eps = eps_list(e);
        
        F1_p = helper.F1_of_D(n,m,T,B,L,delta, D + eps*E);
        lhs = (F1_p(:) - vecF1_base)/eps;   % FD
        
        rhs = J * vecE;                    % analytic
        
        relerr = norm(lhs-rhs)/max(1,norm(lhs));
        fprintf('  eps=%.1e: relerr=%.3e\n', eps, relerr);
    end
end

fprintf('\n=== 検証完了 ===\n');
end
