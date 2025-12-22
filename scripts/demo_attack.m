
function [delta_adv,data_adv, K_adv] = demo_attack(data_ori)% U生成→データ生成→SDP→表示の最小デモ
    % --------------
    [X_adv, Z_adv, U_adv] = attack.dgsm(data_ori);

    data_adv = datasim.SystemData(A,B,X_adv,Z_adv,U_adv,H,E1,E2,Phi11,Phi12,Phi22);

    [sol_adv, K_adv, ~, ~, ~] = solve_quantized_sdp(data_adv);
    % 結果の軽い表示
    delta_adv = (sol_adv.tDelta)^(-1/2);
end