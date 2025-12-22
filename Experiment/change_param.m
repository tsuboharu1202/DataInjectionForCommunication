function new_sd = change_param(sd, epsilon, bound)
    temp_sd = make_sd(epsilon, bound);
    new_sd = temp_sd.changeAB(sd.A,sd.B);
end