function [ Fg_n ] = attract_functn(r , r_g)
%This code takes in r_g is an imput goal location and computes the
%attractive field

delta_r = (r - r_g);
x = delta_r(1);
y = delta_r(2);

%calculate phi for each obstacle
%makes use of only the middle line of the field. Possibly the fastest path
%to destination ---------- commment out code below for original
%case
xo =  r_g(1)-r(1);
yo =  r_g(2)-r(2);
phi = (atan2(-yo, -xo) + (pi));
P = [cos(phi), sin(phi)];
Px = P(1);
Py = P(2);

%for original path (the full attractive field)--------comment out for
%middle line case
% Px = 1;
% Py = 0;

%initialize variables
lambda = 2;


Fg_x =((lambda-1)*Px*(x^2)) + (lambda*Py*x*y) - (Px*(y^2));
Fg_y =((lambda-1)*Py*(y^2)) + (lambda*Px*x*y) - (Py*(x^2));

% !!!!!!!!!!!!!!!!!!!!!!!!!!!
% Fg_x = x^2 - y^2;
% Fg_y = 2*x*y;

Fg = [Fg_x;Fg_y];

%normalixe the forces
norm_Fg = norm(Fg);

%check if the goal is reached
if r == r_g %all(delta_r(:))== 0;
    Fg_n = [0;0];   
else
    Fg_n = Fg/norm_Fg;
    
end
end

