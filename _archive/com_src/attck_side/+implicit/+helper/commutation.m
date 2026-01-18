% function C = commutation(p, q)
% % commutation: commutation matrix C_{p,q} such that vec(A^T) = C vec(A)
% % for A in R^{p×q}. Size: (pq)×(pq).

% C = sparse(p*q, p*q);
% % For each (i,j) in A, maps to (j,i) in A^T.
% % vec(A) stacks columns: index = i + (j-1)*p
% % vec(A^T) index = j + (i-1)*q
% for i = 1:p
%     for j = 1:q
%         row = j + (i-1)*q;
%         col = i + (j-1)*p;
%         C(row, col) = 1;
%     end
% end
% end





