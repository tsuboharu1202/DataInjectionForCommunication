function diag_rank(H)
s = svd(full(H),'econ');
fprintf("sigma_max = %.3e\n", s(1));
fprintf("sigma_min = %.3e\n", s(end));
fprintf("tail (last 10):\n"); disp(s(end-9:end));
fprintf("cond_est = %.3e\n", s(1)/max(s(end),realmin));

% 2通りの tol で数値ランクを見る
tol1 = max(size(H))*eps(s(1));
tol2 = 1e-12*s(1);  % こっちは少し強め
r1 = sum(s > tol1);
r2 = sum(s > tol2);
fprintf("rank(tol=default)=%d, rank(tol=1e-12*s1)=%d\n", r1, r2);
end
