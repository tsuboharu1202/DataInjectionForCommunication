% function [X,Z] = simulate_openloop(A,B,U,x0)
%     % X: n×T,  Z: n×T  (Zは1ステップ先)
%     [n, T] = deal(size(A,1), size(U,2));
%     % 初期値はx0 = 0に基準値を設定
%     if nargin<4, x0 = randn(n,1); end
%     X = zeros(n,T); Z = zeros(n,T);
%     x = x0;
%         for t = 1:T
%             X(:,t) = x;
%             x = A*x + B*U(:,t);
%             Z(:,t) = x;
%         end
% end
