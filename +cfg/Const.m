classdef Const
    properties (Constant)
        %　全体jのステップ数
        SAMPLE_COUNT = 20
        
        % 数値微分の既定ステップ（相対スケール）
        FD_STEP = eps^(1/3)
        
        % 攻撃側制約
        ATTACKER_UPPERLIMIT = 1e-4
        
        
        EPSILON = 1e-3
        NOISE_BOUND = 1e-4
        
        % 乱数生成
        SEED = []
        
        % IDGSM関連の定数
        IDGSM_ALPHA = 1e-5  % ステップサイズ
        IDGSM_STAGNATION_STEPS = 10  % 停滞判定のステップ数
        IDGSM_RHO_CHANGE_THRESHOLD = 1e-6  % 停滞判定の閾値（deltaの変化量）
        IDGSM_RANDOM_NOISE_SCALE = 0.1  % 局所最適解回避のランダムノイズスケール
        IDGSM_ESCAPE_LOCAL_MIN = true  % 局所最適解回避を有効にするか
        MAX_ITERATION = 20  % 最大反復回数
    end
end
