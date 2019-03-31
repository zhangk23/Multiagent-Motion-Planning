function [dij] = worst_dij(agent, neighborsA)
    rate_previous = 0;
    worst_agent_index = 1;
    for i = 1:length(neighborsA)
        rate = calc_rate(agent, neighborsA(i));
        if rate > rate_previous
            worst_agent_index = i;
        end
        rate_previous = rate;
    end
    if isempty(neighborsA)
        dij = 1.5;
    else
        dij = norm(agent.r - neighborsA(worst_agent_index).r);
    end 
end