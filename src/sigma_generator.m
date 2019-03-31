function [ sig ] = sigma_generator(r, r_o, rho, rho_e, rho_o, f)

    dm = 0.8;
    dr = 1;
    dc = 1.25;
    
    dij = norm(r-r_o);
   
    a = -(2/(dr-dc)^3);
    b = 3*(dr+dc)/((dr-dc)^3);
    c = -6*dr*dc/(dr-dc)^3;
    d = (dc^2)*(3*dr-dc)/((dr-dc)^3);
    
    if dij >= dm && dij < dr
        sig = 1;
    elseif dij >= dr && dij < dc
        sig = (a * (dij^3)) + (b * (dij^2)) + (c * dij) + d;    
    else
        sig = 0;
    end
%     
 end



    
