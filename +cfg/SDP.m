classdef SDP
    % SDP - SDP関連パラメータ設定
    %
    % 注意: gamma は実験ごとに明示的に指定してください（デフォルト値は提供しません）
    
    properties (Constant)
        % 保守性の判定マージン（1.01 = 1%マージン）
        TOLERANCE_MARGIN = 1.01
        
        % SDPソルバー設定
        SOLVER_VERBOSE = 0  % 0: 非表示, 1: 表示
    end
end
