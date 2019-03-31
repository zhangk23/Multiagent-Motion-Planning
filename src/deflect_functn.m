function [ Fo_n, lambda ] = deflect_functn(r, ro, r_g)

    delta_r = r - ro;
    
    lambda = 1; % no need to assign lambda, this is only a stakeholder to work with previously written get_si_dot
    
    
     Fo_x = delta_r(1) / norm(delta_r);
     Fo_y = delta_r(2) / norm(delta_r);

     Fo_n = [Fo_x; Fo_y];

end




