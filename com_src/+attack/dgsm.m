function [X_adv, Z_adv, U_adv] = dgsm(sd)
    
    [gradX_pi, gradZ_pi, gradU_pi] = attack.calc_grad(sd);
    [X_adv,Z_adv, U_adv] = attack.make_data_adv(sd,gradX_pi,gradZ_pi,gradU_pi);

end
