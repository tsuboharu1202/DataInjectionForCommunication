classdef Const
    properties (Constant)
        %　全体jのステップ数
        SAMPLE_COUNT = 20
        
        % 数値微分の既定ステップ（相対スケール）
        FD_STEP = eps^(1/3)
        
        % 攻撃側制約
        ATTACKER_UPPERLIMIT = 0.005
        
        
        EPSILON = 0.02
        NOISE_BOUND = 10^(-4)
        
        % 乱数生成
        SEED = []
    end
end
