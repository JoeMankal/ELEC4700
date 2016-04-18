clear
clc

nx = 150;                           % Number of points in the x-direction
ny = 100;							% Numver of points in the y-direction 
W = 1;                              % Width of box in metres (m)
L = 1.5;							% Length of the box in metres (m)
x = linspace(0,L,nx);               % Construct linspace for x
y = linspace(0,W,ny);				% Construct linspace for y
G = zeros(nx*ny);                   % Construct G as a square matrix, nx*ny on a side
Vo = 1;                             % Set boundary condition for x=0
B = zeros((nx*ny),1);				% Construct B matrix

for i = 1:nx
    for j = 1:ny
        n = j + ((i-1)*ny);				% Determine index for diagonal
        if i == 1
            G(n,n) = 1;                	% Boundary condition at x=0
            B(n) = Vo;					% Set value of boundary condition
        elseif i == nx
            G(n,n) = 1;                 % Another boundary colorbar at x = 1
        else
            G(n,n) = -2;				% Set values for non boundary elements
            G(n,n-ny) = 1;
            G(n,n+ny) = 1;
        end
    end
end

V1 = G\B;								% Solve for potential
V = reshape(V1,[ny,nx]);				% Convert the generated vector into a matrix

[X,Y] = meshgrid(x,y);
mesh(X,Y,V)								% Generate 3-D plot of potential
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Potential (V)')
title('Potential Plot')
c = colorbar;
title(c,'Potential (V)')
view([0 0])
print('Part1A','-dpng')