function [si_dot] = get_si_dot (r, tetha, lambda, Fx_val, Fy_val, u)
syms x y;
% r_rg = [0;0];
% rrr=r-r_rg;
% dt = .475;

Px = 1;%exp(rrr(1) * -dt);
Py = 0;%exp(rrr(2) * -dt);
Fx = ((lambda-1)*Px*(x^2)) + (lambda*Py*x*y) - (Px*(y^2));
Fy = ((lambda-1)*Py*(y^2)) + (lambda*Px*x*y) - (Py*(x^2));

%differntiate the forces, F_star
dfx_dx = diff (Fx,x);
dfx_dy = diff (Fx,y);
dfy_dx = diff (Fy,x);
dfy_dy = diff (Fy,y);
ct = cos(tetha);
st = sin(tetha);


a = r(1);
b = r(2);

dfx_dx = subs(dfx_dx,x,a);
dfx_dy = subs(dfx_dy,y,b);
dfy_dx = subs(dfy_dx,y,b);
dfy_dy = subs(dfy_dy,x,a);

%compute si_dot
si_dot = ((((dfy_dx*ct) + (dfy_dy*st))* Fx_val) - (((dfx_dx*ct) + (dfx_dy*st))*Fy_val)) * u;
si_dot = double(si_dot);
end
