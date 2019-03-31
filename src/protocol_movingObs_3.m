clear
close all
beep off

% generate colors
colors = zeros(20,3);
for c = 1:20
    colors(c,:) = [rand,rand,rand];
end

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

A3.rg = [-4,0];
A3.r = [2,3];
A3.theta = 0;

A4.rg = [-4,-2];
A4.r = [4,-3];
A4.theta = 0;

A7.rg = [-2,2];
A7.r = [-10,4];
A7.theta = 0;

A8.rg = [-2,0];
A8.r = [-5,0];
A8.theta = 0;

A9.rg = [-2,-2];
A9.r = [2,0];
A9.theta = 0;

A12.rg = [0,2];
A12.r = [1,-2];
A12.theta = 0;

A15.rg = [0,-4];
A15.r = [0,2];
A15.theta = 0;


% --------------------------------------------------------------
B1.r = [-1,-1];
% B1.theta = -pi/2;
B2.r = [0,0];
B3.r = [3,-3]; 
B3.theta = pi;
B4.r = [-2,3];
B4.theta = -pi/2;

B1.u = 0.06;
B2.u = 0.08;
B3.u = 0.1;
B4.u = 0.2;

obstacles = [B1,B2,B3,B4]; % to run case without obstacles, assign obstacles=[] empty

agents = [A9,A12,A15];

count = 0;
dt = 0.5; % time step for integration (seconds)
T = 10000*dt; % measurement time step

xlim([-5,5]);
ylim([-5,5]);
hold on
% plot goal positions of agents
plot(A9.rg(1),A9.rg(2),'rs')
plot(A12.rg(1),A12.rg(2),'bs')
plot(A15.rg(1),A15.rg(2),'gs')

% plot start positions of agents
plot(A9.r(1),A9.r(2),'ro')
plot(A12.r(1),A12.r(2),'bo')
plot(A15.r(1),A15.r(2),'go')

% plot obstaclestable positions
plot(B1.r(1), B1.r(2),'kx')
plot(B2.r(1), B2.r(2),'kx')
plot(B3.r(1), B3.r(2),'kx')
plot(B4.r(1), B4.r(2),'kx')

vidfile = VideoWriter('3_movingObs.mp4','MPEG-4');
open(vidfile);

for tau = dt : dt : T
    
    % run ClassA agents
    for i = 1:length(agents) 
        % --------------- update neighborsA (and neighborsB) -------------
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
                agents(i).u = max(0, agents(i).min_uik_A);

%             elseif agents(i).dik < dc && agents(i).dik > de
%                 disp('1---2')
%                 agents(i).u =  agents(i).u_e;

%             elseif  agents(i).dik >= dc
            elseif agents(i).dik > de
                 disp('1---3')
                 agents(i).u =  agents(i).u_c;
            else 
                disp(agents(i).dik)
                agents(i).u = 0;
            end
        % -------------------- Scenario 2 ----------------------
        elseif isempty(agents(i).neighborsB)==0 && agents(i).mu == 0
            disp('scenario 2 ...')
            if agents(i).dio <= de && agents(i).dio >= dm
                disp('2---1')
                agents(i).u = agents(i).min_uio;
% 
%                 elseif agreents(i).dio < dc && agents(i).dio > de
%                     disp('2---2')
%                     agents(i).u =  agents(i).u_e;

            elseif  agents(i).dio > de
                disp('2---3')
                agents(i).u =  agents(i).u_c;
            end

        % -------------------- Scenario 3 ----------------------
        elseif isempty(agents(i).neighborsB)==1 && agents(i).mu == 1
            disp('scenario 3 ...')
            if agents(i).dik <= de && agents(i).dik >= dm
                disp('3---1')
                agents(i).u = agents(i).min_uik_A;

            elseif agents(i).dik < dc && agents(i).dik > de
                 disp('3---2')
                 agents(i).u =  agents(i).u_e;

            elseif  agents(i).dik > dc
                 disp('3---3')
                 agents(i).u =  agents(i).u_c;
            end

        % -------------------- Scenario 4 ----------------------
        elseif isempty(agents(i).neighborsB)==0 && agents(i).mu==1
            disp('scenario 4 ...')
            if agents(i).dik <= de && agents(i).dik >= dm
                disp('4---1')
                agents(i).u = agents(i).min_uik_AB;

            elseif agents(i).dik < dc && agents(i).dik > de
                disp('4---2')
                agents(i).u =  agents(i).u_e;

            elseif  agents(i).dik > dc
                disp('4---3')
                agents(i).u =  agents(i).u_c;
            end
        end

        % update agent angular velocity w
        agents(i).w = -agents(i).kw * (agents(i).theta - agents(i).psi) + agents(i).psi_dot;

        ct = cos(agents(i).theta);
        st = sin(agents(i).theta);
        [A] = [ct,0 ; st,0 ; 0,1] * [agents(i).u ; agents(i).w];

        agents(i).r(1) = agents(i).r(1) + A(1) * dt;
        agents(i).r(2) = agents(i).r(2) + A(2) * dt;
        agents(i).theta = wrapToPi(agents(i).theta + A(3)*dt);    
             
        % plot agent positions
        if mod(count,3) == 0
        	plot(agents(i).r(1),agents(i).r(2),'ro');
        elseif mod(count,3) == 1
            plot(agents(i).r(1),agents(i).r(2),'bo');
        elseif mod(count,3) == 2
            plot(agents(i).r(1),agents(i).r(2),'go');
        end
        
        drawnow()
        count = count + 1;

        
        M(i) = getframe;
        writeVideo(vidfile, M(i));

    end
    
    % update/plot ClassB positions
    for obs = 1:length(obstacles)
        obstacles(obs).r_old = obstacles(obs).r;

        obs_w_max = pi/4;
        obs_w_min = -pi/4;
        obs_w_rand = (obs_w_max - obs_w_min)*rand + obs_w_min;
        obstacles(obs).w = obs_w_rand;
        
        obs_u_max = 0.1;  %u_o
        obs_u_min = 0;
        obs_u_rand = (obs_u_max - obs_u_min)*rand + obs_u_min;
        obstacles(obs).u = obs_u_rand;

        ct = cos(obstacles(obs).theta);
        st = sin(obstacles(obs).theta);
        [A_obs] = [ct,0 ; st,0 ; 0,1] * [obstacles(obs).u ; obstacles(obs).w];
        
        obstacles(obs).r(1) = obstacles(obs).r(1) + A_obs(1) * dt;
        obstacles(obs).r(2) = obstacles(obs).r(2) + A_obs(2) * dt;
        obstacles(obs).theta = wrapToPi(obstacles(obs).theta + A_obs(3)*dt);
        
        plot(obstacles(obs).r(1),obstacles(obs).r(2),'kx');
    end
end
close(vidfile);

