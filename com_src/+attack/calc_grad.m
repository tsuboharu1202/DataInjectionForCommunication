function [gradX_pi, gradZ_pi, gradU_pi] = calc_grad(sd,opts)

    h = cfg.Const.FD_STEP;

    % ---- 基準の支配固有値 λ0（複素） ----

    % -------- X の Πλ 勾配 --------
    if isempty(sd.X)
        gradX_pi = [];
    else
        gradX_pi = zeros(size(sd.X));
        for k = 1:numel(sd.X)
            Xp = sd.X; Xp(k) = Xp(k) + h;
            Xm = sd.X; Xm(k) = Xm(k) - h;
            gp = calc_rho(sd.withX(Xp));
            gm = calc_rho(sd.withX(Xm));
            g  = (gp - gm)/(2*h);                         % g: dλ/dX_k  (complex)
            gradX_pi(k) =  g; % Πλ(g)
        end
    end

    % -------- Z の Πλ 勾配 --------
    if isempty(sd.Z)
        gradZ_pi = [];
    else
        gradZ_pi = zeros(size(sd.Z));
        for k = 1:numel(sd.Z)
            Zp = sd.Z; Zp(k) = Zp(k) + h;
            Zm = sd.Z; Zm(k) = Zm(k) - h;
            gp = calc_rho(sd.withZ(Zp));
            gm = calc_rho(sd.withZ(Zm));
            g  = (gp - gm)/(2*h);
            gradZ_pi(k) = g;
        end
    end

    % -------- U の Πλ 勾配 --------
    if isempty(sd.U)
        gradU_pi = [];
    else
        gradU_pi = zeros(size(sd.U));
        for k = 1:numel(sd.U)
            Up = sd.U; Up(k) = Up(k) + h;
            Um = sd.U; Um(k) = Um(k) - h;
            gp = calc_rho(sd.withU(Up));
            gm = calc_rho(sd.withU(Um));
            g  = (gp - gm)/(2*h);
            gradU_pi(k) = g;
        end
    end

    % ===== ローカル関数 =====
    function rho = calc_rho(sysd)
        % SDP → K → 支配固有値 λ（絶対値最大のもの）
        [sol, K, Y, L, diagnostics] = original_thesis.solve_sdp(sysd);
        % tDelta = sol.tDelta;                 % 複素の場合あり
        Ac = sysd.A + sysd.B*K;
        lambda = eig(Ac);
        rho = max(abs(lambda));
    end

end
