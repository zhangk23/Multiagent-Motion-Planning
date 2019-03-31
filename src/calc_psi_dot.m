function [psi, psi_dot] = calc_psi_dot(r, rg, theta, u, obs_r)

    %object radius
    rho = 0.4; %predetermined
    %obstacle radius
    rho_o = 0.4; %predetermined too
    %minimum distance between object and obstacle
    rho_e = 0.8; %predetermined as well
    f = 0.1;

    numberOfObstacles = size(obs_r,1);
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    %gets sigma for all  obs_r with respect to r and store in array
    sig_array = zeros([numberOfObstacles,1]);
    for i = 1 : numberOfObstacles
        r_o = obs_r(i,:);
        sig = sigma_generator(r, r_o, rho, rho_e, rho_o, f);
        sig_array(i) = sig;
    end
    
    %calculates the attractive fields at every point r (by calling the attract fxn)
    sig_finold = 1;
    sig_final = 1;
    for i = 1 : numberOfObstacles
        sig_final = sig_finold*(1-sig_array(i)); % !!!!!!!!!!!!!!!!!!!!!!
        sig_finold = sig_final;
    end
    [Fg_n] = attract_functn(r, rg);
                
    Fg_star = Fg_n * sig_final;
    %disp(Fg_n)
    %calculate the linear velocity of the object

    %call deflect_function to get Fo_n for each obstacle (send in arguments r & r_o)
    %this is a for loop that will generate an array for lambda for each obstacle
    Fo_starold = [0;0];
    Fo_star = [0;0];
    lam_array = zeros([numberOfObstacles,1]);

    for i = 1 : numberOfObstacles
        r_o = obs_r(i,:);
        %allow this function to output lambda as another variable, for use at
        %later to find si_dot
        [Fo_n , lam] = deflect_functn( r, r_o, rg);
        lam_array(i) = lam;
        %get Fo_star
        Fo_star = Fo_starold + (Fo_n * (sig_array(i))); % !!!!!!!!!!!!!!!!!!!!!!
        Fo_starold = Fo_star;
    end

    %add the both attractive and deflective forces
    
    F_star = Fg_star + Fo_star;
    Fx_star = F_star(1);
    Fy_star = F_star(2);

    %get the si_dot to calculate ang vel, w (the attractive half of F_star)
    lambda = 2;
    r_att = r - rg;
    get_si_d = get_si_dot (r_att, theta, lambda, Fx_star, Fy_star, u); 
    si_dot1 = sig_final * get_si_d;
    
    %get the si_dot to calculate ang vel, w (the deflective half of F_star)
    si_dot2old = 0;
    si_dot2 = 0;
    %get ss_dot2       
    for i = 1 : numberOfObstacles
        lam = lam_array(i);
        r_def = r - obs_r(i,:);
        si_dot_ = get_si_dot (r_def, theta, lam, Fx_star, Fy_star, u);
        si_dot2 = si_dot2old + (si_dot_ * (1-sig_array(i)));
        si_dot2old = si_dot2;
    end    

    % get si & si_dot
    psi = (atan2(Fy_star,Fx_star));
%     psi_dot = get_si_dot (r_att, theta, lambda, Fx_star, Fy_star, u); 
    psi_dot = (si_dot1 + si_dot2);
    
end





