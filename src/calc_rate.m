function [roc] = calc_rate(agent_i, agent_j)
    r_ji = agent_i.r - agent_j.r;
    roc = agent_i.u * (r_ji(1)*cos(agent_i.psi) + r_ji(2)*sin(agent_i.psi)) - agent_j.u * (r_ji(1)*cos(agent_j.psi) + r_ji(2)*sin(agent_j.psi));
end