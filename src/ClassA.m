classdef ClassA
    properties
       dm = 0.8;
       de = 0.9;
       dc = 0.98;
       epsilon = 0.5;

       u = 0 % linear velocity
       w = 0 % angular velocity
       r % position (x,y)
       theta % orientation
       rg % goal position (xg, yg)
       ku = 0.075
       kw = 2.5
       obs_r
       
       neighborsA % array containing class A neighbors
       neighborsB % array containing class B neighbors
       
       psi % oritentaion of potential field arctan(Fyi_star/Fxi_star)
       psi_dot % TO DO ...
       u_c %chcek
       u_e % check
       mu % check
       dik % check 
       dio % check
       min_uio % check
       min_uik_A % check
       min_uik_AB
    end
    
    methods
%         function u = get.u(obj)
%             u = obj.ku*tanh(norm(obj.r - obj.rg));
%         end
        function u_c = get.u_c(obj)
            u_c = obj.ku*tanh(norm(obj.r - obj.rg));
        end
        
        function u_e = get.u_e(obj)
            u_e = obj.u_c;
        end
        
        function mu = get.mu(obj)
           mu = 0;
           for i = 1: length(obj.neighborsA)
               if isempty(obj.neighborsA(i).neighborsB) == 0
                   mu = 1;
               end
           end
        end

        function min_uik_A = get.min_uik_A(obj)
            arr_uik = [];
            for i=1:length(obj.neighborsA)
                uik = calc_uik(obj, obj.neighborsA(i));
                arr_uik = horzcat(arr_uik,uik);
            end
            min_uik_A = min(arr_uik);
        end
        
        function min_uik_AB = get.min_uik_AB(obj)
            arr_u = [];
            for i=1:length(obj.neighborsA)
                uik = calc_uik(obj, obj.neighborsA(i));
                arr_u = horzcat(arr_u,uik);
            end
            for j=1:length(obj.neighborsB)
                uio = calc_uio(obj, obj.neighborsB(j));
                arr_u = horzcat(arr_u,uio);
            end
            min_uik_AB = min(arr_u);
        end

        function min_uio = get.min_uio(obj)
            arr_uio = [];
            for i=1:length(obj.neighborsB)
                uio = calc_uio(obj, obj.neighborsB(i));
                arr_uio = horzcat(arr_uio,uio);
            end
            min_uio = min(arr_uio);
        end
        
        function dik = get.dik(obj)
            dik = worst_dij(obj, obj.neighborsA);
        end
        
        function dio = get.dio(obj)
            dio = worst_dio(obj, obj.neighborsB);
        end
                
        function obs_r = get.obs_r(obj)
            obs_r = [];
            for i = 1: length(obj.neighborsA)
            	obs_r = vertcat(obs_r, obj.neighborsA(i).r);
            end
            for j = 1: length(obj.neighborsB)
            	obs_r = vertcat(obs_r, obj.neighborsB(j).r);
            end
        end
       
        function psi = get.psi(obj)
           [one,two] = calc_psi_dot(obj.r, obj.rg, obj.theta, obj.u, obj.obs_r);
           psi = one;
        end
        
        function psi_dot = get.psi_dot(obj)
           [one,two] = calc_psi_dot(obj.r, obj.rg, obj.theta, obj.u, obj.obs_r);
           psi_dot = two;
        end

    end
end