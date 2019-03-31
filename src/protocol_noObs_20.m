clear
close all

% generate colors
colors = zeros(20,3);
for c = 1:20
    colors(c,:) = [rand,rand,rand];
end

colors(9,:) = [0.8500, 0.3250, 0.0980];
colors(10,:) = [0.75, 0, 0.75];
colors(13,:) = [0.6350, 0.0780, 0.1840]	;
colors(15,:) = [0.75, 0.75, 0];


% define some constant distances
dm = 0.8;
de = 0.9;
dc = 0.98;
Rc = 1.25;

% create class A agents
A1 = ClassA;
A2 = ClassA;
A3 = ClassA;
A4 = ClassA;
A5 = ClassA;
A6 = ClassA;
A7 = ClassA;
A8 = ClassA;
A9 = ClassA;
A10 = ClassA;
A11 = ClassA;
A12 = ClassA;
A13 = ClassA;
A14 = ClassA;
A15 = ClassA;
A16 = ClassA;
A17 = ClassA;
A18 = ClassA;
A19 = ClassA;
A20 = ClassA;

% create class B agents
B1 = ClassB;
B2 = ClassB;
B3 = ClassB;
B4 = ClassB;

% set agent locations 
A1.rg = [-4,4];
A1.r = [-6,-3.5];
A1.theta = 0;

A2.rg = [-4,2];
A2.r = [1,-3];
% A2.r = [1,-6];
A2.theta = 0;

A3.rg = [-4,0];
A3.r = [2,3];
A3.theta = 0;

A4.rg = [-4,-2];
A4.r = [4,-3];
A4.theta = 0;

A5.rg = [-4,-4];
A5.r = [3,3];
% A5.r = [5,5];
A5.theta = 0;

A6.rg = [-2,4];
A6.r = [0,4];
A6.theta = 0;

A7.rg = [-2,2];
A7.r = [-10,4];
A7.theta = 0;

A8.rg = [-2,0];
A8.r = [-5,0];
A8.theta = 0;

A9.rg = [-2,-2];
A9.r = [2,0];
A9.theta = 0;

A10.rg = [-2,-4];
A10.r = [6,3];
A10.theta = 0;

A11.rg = [0,4];
A11.r = [0,6];
A11.theta = 0;

A12.rg = [0,2];
A12.r = [1,-2];
A12.theta = 0;

A13.rg = [0,0];
A13.r = [-5,-6];
A13.theta = 0;

A14.rg = [0,-2];
A14.r = [-4,6];
A14.theta = 0;

A15.rg = [0,-4];
A15.r = [0,2];
A15.theta = 0;

A16.rg = [2,4];
A16.r = [-10,-4];
A16.theta = 0;

A17.rg = [2,2];
A17.r = [0,8];
A17.theta = 0;

A18.rg = [2,0];
A18.r = [-10,0];
A18.theta = 0;

A19.rg = [2,-2];
A19.r = [-1,-8];
A19.theta = 0;

A20.rg = [2,-4];
A20.r = [-8,4];
A20.theta = 0;

% --------------------------------------------------------------
% B1.r = [-2,-3];
B1.r = [-2,-6];
B2.r = [-6,3];
B3.r = [-1,3];
B4.r = [4,0];
obstacles = [];

agents = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20];

count_agent = 0;
count = 0;
dt = 0.5; % time step for integration (seconds)
T = 10000*dt; % measurement time step

xlim([-15,10]);
ylim([-10,10]);
hold on
% plot goal positions of agents

plot(A1.rg(1),A1.rg(2),'rs')
plot(A2.rg(1),A2.rg(2),'bs')
plot(A3.rg(1),A3.rg(2),'gs')
plot(A4.rg(1),A4.rg(2),'ys')
plot(A5.rg(1),A5.rg(2),'cs')
plot(A6.rg(1),A6.rg(2),'s','Color',colors(1,:))
plot(A7.rg(1),A7.rg(2),'s','Color',colors(2,:))
plot(A8.rg(1),A8.rg(2),'s','Color',colors(3,:))
plot(A9.rg(1),A9.rg(2),'s','Color',colors(4,:))
plot(A10.rg(1),A10.rg(2),'s','Color',colors(5,:))
plot(A11.rg(1),A11.rg(2),'s','Color',colors(6,:))
plot(A12.rg(1),A12.rg(2),'s','Color',colors(7,:))
plot(A13.rg(1),A13.rg(2),'s','Color',colors(8,:))
plot(A14.rg(1),A14.rg(2),'s','Color',colors(9,:))
plot(A15.rg(1),A15.rg(2),'s','Color',colors(10,:))
plot(A16.rg(1),A16.rg(2),'s','Color',colors(11,:))
plot(A17.rg(1),A17.rg(2),'s','Color',colors(12,:))
plot(A18.rg(1),A18.rg(2),'s','Color',colors(13,:))
plot(A19.rg(1),A19.rg(2),'s','Color',colors(14,:))
plot(A20.rg(1),A20.rg(2),'s','Color',colors(15,:))


% plot obstaclestable positions
% B1 = plot(B1.r(1), B1.r(2),'kx');
% B2 = plot(B2.r(1), B2.r(2),'kx');
% plot(B3.r(1), B3.r(2),'kx')
% plot(B4.r(1), B4.r(2),'kx')

% legend([goal2 goal3 goal5 goal7 goal9 goal12 goal14 goal15 goal18 goal20 B1 B2],{'(-4,2)','(-4,0)','(-4,-4)','(-2,2)','(-2,-2)','(0,2)','(0,-2)','(0,-4)','(2,0)','(2,-4)','(-7,-0.5)','(1,4)'},'AutoUpdate','off');

% plot start positions of agents
plot(A1.r(1), A1.r(2),'ro')
plot(A2.r(1), A2.r(2),'bo')
plot(A3.r(1), A3.r(2),'go')
plot(A4.r(1), A4.r(2),'yo')
plot(A5.r(1), A5.r(2),'co')
plot(A6.r(1),A6.r(2),'o','Color',colors(1,:))
plot(A7.r(1),A7.r(2),'o','Color',colors(2,:))
plot(A8.r(1),A8.r(2),'o','Color',colors(3,:))
plot(A9.r(1),A9.r(2),'o','Color',colors(4,:))
plot(A10.r(1),A10.r(2),'o','Color',colors(5,:))
plot(A11.r(1),A11.r(2),'o','Color',colors(6,:))
plot(A12.r(1),A12.r(2),'o','Color',colors(7,:))
plot(A13.r(1),A13.r(2),'o','Color',colors(8,:))
plot(A14.r(1),A14.r(2),'o','Color',colors(9,:))
plot(A15.r(1),A15.r(2),'o','Color',colors(10,:))
plot(A16.r(1),A16.r(2),'o','Color',colors(11,:))
plot(A17.r(1),A17.r(2),'o','Color',colors(12,:))
plot(A18.r(1),A18.r(2),'o','Color',colors(13,:))
plot(A19.r(1),A19.r(2),'o','Color',colors(14,:))
plot(A20.r(1),A20.r(2),'o','Color',colors(15,:))

vidfile = VideoWriter('20_noObs.mp4','MPEG-4');
open(vidfile);
count_agent_11 = 0;
count_agent_21 = 0;
count_agent_31 = 0;
count_agent_41 = 0;
for tau = dt : dt : T
    data_pos = [];
    for i = 1:length(agents)
        % --------------- update neighborsA and neighborsB -------------
        agents(i).neighborsA = [];
        agents(i).neighborsB = [];
        other_agents = agents;
        other_agents(i) = [];

        for j = 1:length(other_agents)
            if norm(agents(i).r - other_agents(j).r) <= Rc
                agents(i).neighborsA = horzcat(agents(i).neighborsA,other_agents(j));
            end
        end

        for k = 1:length(obstacles)
            if norm(agents(i).r - obstacles(k).r) <= Rc
                agents(i).neighborsB = horzcat(agents(i).neighborsB,obstacles(k));
            end
        end

        % start protocol to update agent linear velocity u
        % -------------------- Scenario 1 ----------------------
        % note that for all four scenarios, the second sub-scenario is
        % commented out (refer to report)
        if isempty(agents(i).neighborsB)==1 && agents(i).mu == 0
            disp('scenario 1 ...')
            if agents(i).dik <= de && agents(i).dik >= dm
                disp('1---1')
                count_agent_11 = count_agent_11 + 1;
                agents(i).u = max(0,agents(i).min_uik_A);
                disp(agents(i).min_uik_A);
%                     disp(agents(i).u);
%             elseif agents(i).dik < dc && agents(i).dik > de
%                 disp('1---2')
%                 agents(i).u =  agents(i).u_e;

%             elseif  agents(i).dik >= dc
            elseif agents(i).dik > de
                 disp('1---3')
                 agents(i).u =  agents(i).u_c;
            else
                disp('1----unknown')
                disp(agents(i).dik)
            end
        % -------------------- Scenario 2 ----------------------
        elseif isempty(agents(i).neighborsB)==0 && agents(i).mu == 0
            disp('scenario 2 ...')
            if agents(i).dio <= de && agents(i).dio >= dm
                disp('2---1')
                count_agent_21 = count_agent_21 + 1;
                agents(i).u = agents(i).min_uio;

%             elseif agents(i).dio < dc && agents(i).dio > de
%                 disp('2---2')
%                 agents(i).u =  agents(i).u_e;

            elseif  agents(i).dio > de
                disp('2---3')
                agents(i).u =  agents(i).u_c;
            else
                disp('2----unknown')
                disp(agents(i).dio)
            end

        % -------------------- Scenario 3 ----------------------
        elseif isempty(agents(i).neighborsB)==1 && agents(i).mu == 1
            disp('scenario 3 ...')
            if agents(i).dik <= de && agents(i).dik >= dm
                disp('3---1')
                count_agent_31 = count_agent_31 + 1;
                agents(i).u = agents(i).min_uik_A;

%             elseif agents(i).dik < dc && agents(i).dik > de
%                  disp('3---2')
%                  agents(i).u =  agents(i).u_e;

            elseif  agents(i).dik > de
                 disp('3---3')
                 agents(i).u =  agents(i).u_c;
            else
                disp('3----unknown')
                disp(agents(i).dik)
            end

        % -------------------- Scenario 4 ----------------------
        elseif isempty(agents(i).neighborsB)==0 && agents(i).mu==1
            disp('scenario 4 ...')
            if agents(i).dik <= de && agents(i).dik >= dm
                disp('4---1')
                count_agent_41 = count_agent_41 + 1;
                agents(i).u = agents(i).min_uik_AB;
                disp(agents(i).u);

%             elseif agents(i).dik < dc && agents(i).dik > de
%                 disp('4---2')
%                 agents(i).u =  agents(i).u_e;

            elseif  agents(i).dik > de
                disp('4---3')
                agents(i).u =  agents(i).u_c;
            else
                disp('4----unknown')
                disp(agents(i).dik)
            end
        end

        % update agent angular velocity w
        agents(i).w = -agents(i).kw * (agents(i).theta - agents(i).psi) + agents(i).psi_dot;

        % update agent position based on new u, w and theta
        ct = cos(agents(i).theta);
        st = sin(agents(i).theta);
        [A] = [ct,0 ; st,0 ; 0,1] * [agents(i).u ; agents(i).w];

        agents(i).r(1) = agents(i).r(1) + A(1) * dt;
        agents(i).r(2) = agents(i).r(2) + A(2) * dt;
        agents(i).theta = wrapToPi(agents(i).theta + A(3)*dt);
        
        % plot agent positions
        if mod(count_agent,20) == 0
        	plot(agents(i).r(1),agents(i).r(2),'ro');
        elseif mod(count_agent,20) == 1
            plot(agents(i).r(1),agents(i).r(2),'bo');
        elseif mod(count_agent,20) == 2
            plot(agents(i).r(1),agents(i).r(2),'go');
        elseif mod(count_agent,20) == 3
            plot(agents(i).r(1),agents(i).r(2),'yo');
        elseif mod(count_agent,20) == 4
            plot(agents(i).r(1),agents(i).r(2),'co');
        elseif mod(count_agent,20) == 5
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(1,:));
        elseif mod(count_agent,20) == 6
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(2,:));
        elseif mod(count_agent,20) == 7
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(3,:));
        elseif mod(count_agent,20) == 8
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(4,:));
        elseif mod(count_agent,20) == 9
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(5,:));
        elseif mod(count_agent,20) == 10
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(6,:));
        elseif mod(count_agent,20) == 11
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(7,:));
        elseif mod(count_agent,20) == 12
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(8,:));
        elseif mod(count_agent,20) == 13
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(9,:));
        elseif mod(count_agent,20) == 14
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(10,:));
        elseif mod(count_agent,20) == 15
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(11,:));
        elseif mod(count_agent,20) == 16
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(12,:));
        elseif mod(count_agent,20) == 17
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(13,:));
        elseif mod(count_agent,20) == 18
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(14,:));
        elseif mod(count_agent,20) == 19
            plot(agents(i).r(1),agents(i).r(2),'o','Color',colors(15,:));
        end
        
        drawnow()
        count_agent = count_agent + 1;
        data_pos = vertcat(data_pos, agents(i).r);
   
    end
        
    % record agent data every 40 iterations (each row represents one agent's position)
    if mod(count,40)==0
        filename = sprintf('noObs_iteration%d.xls',count);
        xlswrite(filename, data_pos);
    end
   
    count = count + 1;
    M(count) = getframe;
    writeVideo(vidfile, M(count));
end
close(vidfile);


