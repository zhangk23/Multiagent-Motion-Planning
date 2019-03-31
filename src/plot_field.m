[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);

lamda = 2;
Px = 0;
Py = 1;
% c = (x.^2 + y.^2)^.5;
% u = ((lamda-1)*Px*(x.^2)) + (lamda*Py*x.*y) - (Px*(y.^2)) ;
% v = ((lamda-1)*Py*(y.^2)) + (lamda*Px*x.*y) - (Py*(x.^2)) ;

u = x.^2 - y.^2;
v = 2.*x.*y;

% u = -y.^2;
% v = x.*y;

% u = x;
% v = y;

figure
quiver(x,y,u,v,'LineWidth',1)

startx = -1:0.1:1;
starty = ones(size(startx));
% streamline(x,y,u,v,startx,starty)
axis equal
axis([-1 1 -1 1])
