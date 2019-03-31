function [roc] = calc_obs_rate(agent_i, agent_o)

    if isempty(agent_o.r_old) == 1
        obs_u = 0;
    else
        obs_u = norm(agent_o.r - agent_o.r_old)/0.5;
    end
    
%     disp('---obs linear speed------')
%     disp(obs_u)
%     disp('--------------------------')

    r_ji = agent_i.r - agent_o.r;
    roc = obs_u * (r_ji(1)*cos(agent_i.psi) + r_ji(2)*sin(agent_i.psi)) - obs_u * (r_ji(1)*cos(agent_o.psi) + r_ji(2)*sin(agent_o.psi));

end
