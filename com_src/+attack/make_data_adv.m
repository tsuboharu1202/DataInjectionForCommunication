function [X_adv, Z_adv, U_adv] = make_data_adv(sd, X_grad,Z_grad, U_grad)
X_adv = d_adv_add(sd.X, X_grad);
Z_adv = d_adv_add(sd.Z, Z_grad);
U_adv = d_adv_add(sd.U, U_grad);

    function d_adv = d_adv_add(d_ori, d_grad)
        
        eps_att = cfg.Const.ATTACKER_UPPERLIMIT;   % スカラー
        % 要素単位で： d_adv = d_ori + ε * |d_ori| * sign(real(d_grad))
        % 勾配が複素数の場合でも、実部の符号を使う
        d_adv = d_ori + eps_att* abs(d_ori) .* sign(real(d_grad));
        
    end
end