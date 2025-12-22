classdef SystemData
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
            obj = obj.validate();   % 生成時チェック
        end
        
        % =======================
        % 不変っぽい更新API
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
        
        function [sol, K, Y, L, diagnostics] = solve_sdp_on_data(obj)
            [sol, K, Y, L, diagnostics] = original_thesis.solve_sdp(obj);
        end
        
        function [spectral, rho] = calc_spectral_rho(obj)
            [sol, K, ~, ~, ~] = original_thesis.solve_sdp(obj);
            spectral = max(abs(eig(obj.A + obj.B*K)));
            rho = sol.rho;
        end
        
        % 元のシステム行列Aの不安定な固有値の絶対値の積を計算し、逆数を取る
        % これが最も荒い量子化下界
        function rho_rough = rough_quantization_lower_bound(obj)
            eigA = eig(obj.A);
            unstable_eigA = abs(eigA(abs(eigA) > 1));
            rho_rough = prod(unstable_eigA);
            rho_rough = 1/rho_rough;
        end
        
    end
    
    
    methods (Access = private)
        function obj = validate(obj)
            % 正方・対称性
            if ~isempty(obj.A)
                assert(size(obj.A,1)==size(obj.A,2), 'A must be square.');
            end
            
            % 次元整合
            if ~isempty(obj.A) && ~isempty(obj.B)
                n = size(obj.A,1);
                assert(size(obj.B,1)==n, 'B rows must match A.');
                if ~isempty(obj.X), assert(size(obj.X,1)==n, 'X rows must match A.'); end
                if ~isempty(obj.Z), assert(size(obj.Z,1)==n, 'Z rows must match A.'); end
                if ~isempty(obj.U), assert(size(obj.U,1)==size(obj.B,2), 'U rows must match input dim.'); end
                
                if ~isempty(obj.X) && ~isempty(obj.U)
                    assert(size(obj.X,2) == size(obj.U,2), 'X cols must match U cols.');
                end
                if ~isempty(obj.Z) && ~isempty(obj.U)
                    assert(size(obj.Z,2) == size(obj.U,2), 'Z cols must match U cols.');
                end
            end
        end
        
    end
end
