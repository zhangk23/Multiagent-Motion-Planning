function [dij] = worst_dio(agent, neighborsB)  
    rate_previous = 0;
    worst_agent_index = 1;
    for i = 1:length(neighborsB)
        rate = calc_obs_rate(agent, neighborsB(i));
        if rate > rate_previous
            worst_agent_index = i;
        end
        rate_previous = rate;
    end
    if isempty(neighborsB)
        dij = 1.5;
    else
        dij = norm(agent.r - neighborsB(worst_agent_index).r);
    end 
end