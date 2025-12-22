function data_ori = make_sd(epsilon, noise_bound)
    % 1) 連続の種 or 再現性
    rng();

    % 2) システム＆重み
    [n,m,T] = deal(3,1,cfg.Const.SAMPLE_COUNT);
    [A,B] = datasim.make_lti();

    % 3) 入力とデータ取得
    V = make_inputU(m);
    [X,Z,U] = datasim.simulate_openloop_stable(A,B,V);
    % U = V;
    % [X,Z] = datasim.simulate_openloop(A,B,U);
    visualize.plot_data(X);

    H  = epsilon*[1;1;1]; 
    E1 = [1,1,1];       
    E2 = 1;       

    Phi11 = noise_bound*cfg.Const.SAMPLE_COUNT * eye(n);
    Phi12 = zeros(n,T);
    Phi22 = -eye(T);    

    % === ここから下はあなたの元コードのまま ===
    data_ori = datasim.SystemData(A,B,X,Z,U,H,E1,E2,Phi11,Phi12,Phi22);