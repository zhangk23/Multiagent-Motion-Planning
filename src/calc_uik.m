function [uik] = calc_uik(agentA_i, agentA_k)
    % calculates u_i|k, the safe velocity if agent i w.r.t a class-A agent
    
%     disp('-----calc_uik-------')
    dm = 0.8;
    de = 1;
    
    d_ik = norm(agentA_i.r - agentA_k.r);
    r_ki = agentA_i.r - agentA_k.r;
    
%     disp('------u_e------')
%     disp(agentA_i.u_e)
%     disp('-----r_ki------')
%     disp(r_ki)
%     disp('-----agentA_k.u-------')
%     disp(agentA_k.u)
%     disp('-----i_psi-------')
%     disp(agentA_k.psi)
%     disp(agentA_i.psi)
    
    u_s = (agentA_k.u) * (r_ki(1)*cos(agentA_k.psi)+r_ki(2)*sin(agentA_k.psi)) / (r_ki(1)*cos(agentA_i.psi)+r_ki(2)*sin(agentA_i.psi));

    uik = agentA_i.u_e * (d_ik - dm)/(de - dm) + agentA_i.epsilon * u_s+(de - d_ik)/(de - dm);
    
%     disp('-----u_s-------')
%     disp(u_s)
%     disp('-----uik-------')
%     disp(uik)
%     disp('---------------')
    
end