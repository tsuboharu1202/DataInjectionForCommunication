classdef AttackType
    % AttackType: 攻撃手法の種類
    %
    % 攻撃方向（deltaを大きくする/小さくする）は opts.direction で指定してください。
    %   'positive': deltaを大きくする
    %   'negative': deltaを小さくする（デフォルト）
    
    properties (Constant)
        % 直接法（1回の勾配計算）
        DIRECT_DGSM_DELTA = "DIRECT_DGSM_DELTA"
        
        % 反復法（複数回の勾配計算）
        IMPLICIT_IDGSM_DELTA = "IMPLICIT_IDGSM_DELTA"
    end
end
