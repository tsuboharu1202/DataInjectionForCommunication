classdef SystemData
    % SystemData: システムデータを保持するクラス
    %
    % プロパティ:
    %   A, B: システム行列
    %   X, Z, U: 軌道データ（X: 状態、Z: 次状態、U: 入力）
    %   Phi11, Phi12, Phi22: S-procedureの重み行列
    
    properties (SetAccess = private)
        A double = []
        B double = []
        X double = []
        Z double = []
        U double = []
        
        Phi11 double = []
        Phi12 double = []
        Phi22 double = []
    end
    
    methods
        function obj = SystemData(A,B,X,Z,U,Phi11,Phi12,Phi22)
            arguments
                A double = []
                B double = []
                X double = []
                Z double = []
                U double = []
                
                Phi11 double = []
                Phi12 double = []
                Phi22 double = []
            end
            obj.A = A; obj.B = B;
            obj.X = X; obj.Z = Z; obj.U = U;
            obj.Phi11 = Phi11; obj.Phi12 = Phi12; obj.Phi22 = Phi22;
            obj = obj.validate();
        end
        
        % =======================
        % 更新メソッド（イミュータブル風）
        % =======================
        function newObj = changeAB(obj, A,B)
            newObj = datasim.SystemData(A,B,obj.X,obj.Z,obj.U,obj.Phi11,obj.Phi12,obj.Phi22);
        end
        
        function newObj = withX(obj, newX)
            newObj = obj;
            newObj.X = newX;
            newObj = newObj.validate();
        end
        
        function newObj = withU(obj, newU)
            newObj = obj;
            newObj.U = newU;
            newObj = newObj.validate();
        end
        
        function newObj = withZ(obj, newZ)
            newObj = obj;
            newObj.Z = newZ;
            newObj = newObj.validate();
        end
        
        function newObj = withXZU(obj,newX,newZ,newU)
            newObj = obj;
            newObj.X = newX;
            newObj.Z = newZ;
            newObj.U = newU;
            newObj = newObj.validate();
        end
        
        function newObj = withDefaultPhi(obj, phi11_coef)
            % withDefaultPhi: デフォルトのPhi行列を設定
            %   phi11_coef: Phi11 = phi11_coef * eye(n) の係数（デフォルト: 0.1）
            if nargin < 2 || isempty(phi11_coef)
                phi11_coef = 0.1;
            end
            n = size(obj.A, 1);
            T = size(obj.X, 2);
            newObj = obj;
            newObj.Phi11 = phi11_coef * eye(n);
            newObj.Phi12 = zeros(n, T);
            newObj.Phi22 = -eye(T);
        end
        
        % =======================
        % SDP関連メソッド
        % =======================
        function [sol, K, Y, L, diagnostics] = solve_sdp_on_data(obj)
            % solve_sdp_on_data: baseline SDPを解く
            [sol, K, Y, L, diagnostics] = baseline.solve_sdp(obj);
        end
        
        function [spectral, rho] = calc_spectral_rho(obj)
            % calc_spectral_rho: スペクトル半径とrhoを計算
            [sol, K, ~, ~, ~] = baseline.solve_sdp(obj);
            spectral = max(abs(eig(obj.A + obj.B*K)));
            rho = sol.rho;
        end
        
        function rho_rough = rough_quantization_lower_bound(obj)
            % rough_quantization_lower_bound: 理論的な量子化下界を計算
            %   rho_rough = 1 / prod(|不安定固有値|)
            rho_rough = cfg.System.calcRhoRough(obj.A);
        end
    end
    
    methods (Static)
        function obj = create(A, B, X, Z, U, phi11_coef)
            % create: システムデータを簡単に生成するファクトリーメソッド
            %   Phi行列は自動生成される
            %
            % 使用例:
            %   data = datasim.SystemData.create(A, B, X, Z, U);
            %   data = datasim.SystemData.create(A, B, X, Z, U, 0.1);  % Phi11係数指定
            
            if nargin < 6 || isempty(phi11_coef)
                phi11_coef = 0.1;
            end
            
            n = size(A, 1);
            T = size(X, 2);
            Phi11 = phi11_coef * eye(n);
            Phi12 = zeros(n, T);
            Phi22 = -eye(T);
            
            obj = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
        end
    end
    
    methods (Access = private)
        function obj = validate(obj)
            % validate: データの整合性をチェック
            
            % A は正方行列
            if ~isempty(obj.A)
                if size(obj.A,1) ~= size(obj.A,2)
                    error('SystemData:InvalidA', 'A は正方行列である必要があります。');
                end
            end
            
            % 次元整合性チェック
            if ~isempty(obj.A) && ~isempty(obj.B)
                n = size(obj.A,1);
                
                if size(obj.B,1) ~= n
                    error('SystemData:InvalidB', 'B の行数は A の次元と一致する必要があります。');
                end
                if ~isempty(obj.X) && size(obj.X,1) ~= n
                    error('SystemData:InvalidX', 'X の行数は状態次元と一致する必要があります。');
                end
                if ~isempty(obj.Z) && size(obj.Z,1) ~= n
                    error('SystemData:InvalidZ', 'Z の行数は状態次元と一致する必要があります。');
                end
                if ~isempty(obj.U) && size(obj.U,1) ~= size(obj.B,2)
                    error('SystemData:InvalidU', 'U の行数は入力次元と一致する必要があります。');
                end
                
                if ~isempty(obj.X) && ~isempty(obj.U)
                    if size(obj.X,2) ~= size(obj.U,2)
                        error('SystemData:DimensionMismatch', 'X と U のサンプル数が一致しません。');
                    end
                end
                if ~isempty(obj.Z) && ~isempty(obj.U)
                    if size(obj.Z,2) ~= size(obj.U,2)
                        error('SystemData:DimensionMismatch', 'Z と U のサンプル数が一致しません。');
                    end
                end
            end
        end
    end
end
