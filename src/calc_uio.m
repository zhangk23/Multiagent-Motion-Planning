function [uio] = calc_uio(agentA_i, agentB_o)

%     disp('-------calc_uio---------')
%     disp('-------------------')
    
    dm = 0.8;
    dc = 1.25;
    
    if isempty(agentB_o.r_old)
        obs_u = 0;
    else
        obs_u = norm(agentB_o.r - agentB_o.r_old)/0.5;
    end  
    
    d_io = norm(agentA_i.r - agentB_o.r);
    r_io = agentB_o.r - agentA_i.r;
    
    u_s = obs_u * dc / (r_io(1)*cos(agentA_i.psi)+ r_io(2)*sin(agentA_i.psi));

    uio = agentA_i.u_c * (d_io - dm)/(dc - dm) + u_s+(dc - d_io)/(dc- dm);

end