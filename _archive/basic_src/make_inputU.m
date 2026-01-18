% function U = make_inputU(m)
% % U = m× crg.Const.SAMPLE_COUNT
% % 行列の生成。それぞれの要素は標準正規分布から作成
%     % U = randn(m, cfg.Const.SAMPLE_COUNT);
%     T = cfg.Const.SAMPLE_COUNT;
%     U = zeros(m,T); U(:,1) = randn(m,1);
%     for t = 2:T
%         st = 0.3;
%         U(:,t) = U(:,t-1) + st*randn(m,1);
%     end
% end