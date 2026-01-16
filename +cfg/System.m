classdef System
    % System - デフォルトシステム設定
    %
    % 使用例:
    %   A = cfg.System.A;
    %   B = cfg.System.B;
    %   [n, m] = cfg.System.getDimensions();
    
    properties (Constant)
        % デフォルトシステム行列（3次元1入力）
        A = [-0.192, -0.936, -0.814;
            -0.918,  0.729, -0.724;
            -0.412, -0.135, -0.516]
        
        B = [-0.554; 0.735; 0.528]
        
        % システム次元
        n = 3  % 状態次元
        m = 1  % 入力次元
    end
    
    methods (Static)
        function [n, m] = getDimensions()
            % システム次元を取得
            n = cfg.System.n;
            m = cfg.System.m;
        end
        
        function Phi = getDefaultPhi(n, T, coef)
            % デフォルトのPhi行列を取得
            %   coef: Phi11の係数（デフォルト: 1e-1）
            if nargin < 3
                coef = 1e-1;
            end
            Phi.Phi11 = coef * eye(n);
            Phi.Phi12 = zeros(n, T);
            Phi.Phi22 = -eye(T);
        end
        
        function rho = calcRhoRough(A)
            % 理論的下界 rho_rough を計算
            %   rho_rough = 1 / prod(|不安定固有値|)
            eigA = eig(A);
            gamma_prod = 1;
            for i = 1:numel(eigA)
                if abs(eigA(i)) > 1
                    gamma_prod = gamma_prod * abs(eigA(i));
                end
            end
            rho = 1 / gamma_prod;
        end
    end
end
