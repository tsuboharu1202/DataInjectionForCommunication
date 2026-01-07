classdef Const
    properties (Constant)
        %　全体jのステップ数
        SAMPLE_COUNT = 20
        
        % 数値微分の既定ステップ（相対スケール）
        FD_STEP = eps^(1/3)
        
        % 攻撃側制約
        ATTACKER_UPPERLIMIT = 1e-2
        
        
        EPSILON = 1e-4
        
        % 乱数生成
        SEED = []
        
        % IDGSM関連の定数
        IDGSM_ALPHA = 1e-3  % ステップサイズ
        MAX_ITERATION = 40  % 最大反復回数
    end
end
