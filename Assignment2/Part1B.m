clear
clc

nx = 150;                           % Number of points in the x-direction
ny = 100;							% Number of points in the y-direction
W = 1;                              % Width of box in metres (m)
L = 1.5;							% Length of box in metres (m)
x = linspace(0,L,nx);               % Construct linspace for x
y = linspace(0,W,ny);				% Construct linspace for y

G = zeros((nx*ny));                 % Construct G square matrix, nx*ny on a side
Vx = 1;                            	% Set boundary condition for x=0 and x=1
Vy = 0;								% Set boundary condition for y=0 and y=1
B = zeros((nx*ny),1);				% Construct B matrix

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;				% Determine index for diagonal
        if i == 1
            G(n,n) = 1;                	% Set boundary condition for x=0
            B(n) = Vx;
        elseif i == nx
            G(n,n) = 1;					% Set boundary condition for x=1
            B(n) = Vx;
        elseif j == 1					% Set boundary condition for y=0
            G(n,n) = 1;
            B(n) = Vy;
        elseif j == ny					% Set boundary condition for y=1
            G(n,n) = 1;
            B(n) = Vy;
        else							% Set values for non boundary elements
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;
        end
    end
end

V1 = G\B;								% Solve for potentials
V = reshape(V1,[ny,nx]);				% Convert generated vector into a matrix

Analytical_V = zeros(ny,nx);			% Create matrix for analytical solution
b = L/2;								% Adjust W for analytical solution
x2 = linspace(-b,b,nx);					% Create linspace for adjusted width

for i = 1:nx							% Perform analytical solution
    for j = 1:ny
        for n = 1:2:151
            Analytical_V(j,i) = Analytical_V(j,i) + (4*Vx*cosh(n*pi*x2(i)/W)*sin(n*pi*y(j)/W)/(n*pi*cosh(n*pi*b/W)));
        end
    end
end

[X,Y] = meshgrid(x,y);
[X2,Y2] = meshgrid(x2,y);

figure(1)								% Plot numerical solution
mesh(X,Y,V)
axis([0 L 0 W 0 Vx])
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Potential (V)')
title('Potential Plot using Numerical Techniques')
title(c,'Potential (V)')
print('Part1B_1','-dpng')

figure(2)								% Plot analytical solution
mesh(X2,Y2,Analytical_V)
axis([-b b 0 W 0 1.2])
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Potential (V)')
title('Potential Plot using the Analytical Solution')
title(c,'Potential (V)')
print('Part1B_2','-dpng')