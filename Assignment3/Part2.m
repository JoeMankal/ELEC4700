clear
clc

sigma_outside = 1;					% Set conductance for outside of the boxes
sigma_inside = 0.01;				% Set conductance for inside of the boxes
L = 200e-9;							% Set length of space in metres (m)
W = 100e-9;							% Set width of space in metres (m)
Lb = L*0.35;						% Set length of box
Wb = W*0.35;						% Set width of box
nx = 150;							% Number of points in x direction
ny = 100;							% Number of points in y direction
x = linspace(0,L,nx);               % Construct linspace for x
y = linspace(0,W,ny);				% Construct linspace for y
dx = x(2) - x(1);
dy = y(2) - y(1);
sigma = zeros(ny,nx);               % Construct matrix for conductance
G = zeros((nx*ny));                 % Construct G square matrix, nx*ny on a side
Vx = 0.1;                           % Set boundary condition for x=0
B = zeros((nx*ny),1);				% Construct B matrix

for i = 1:nx						% Set conductance for each point whether it is within a box or not
    for j = 1:ny
        if (y(j) <= Wb) && (x(i) >= Lb) && (x(i) <= L-Lb) ||...
                (y(j) >= W-Wb) && (x(i) >= Lb) && (x(i) <= L-Lb)
            sigma(j,i) = sigma_inside;
        else 
            sigma(j,i) = sigma_outside;
        end   
    end
end

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;											% Determine index for diagonal
        if i == 1													% Set boundary condition for x = 0		
            G(n,n) = 1;
            B(n) = Vx;
        elseif i == nx												% Set boundary condition for x = 1
            G(n,n) = 1;
        elseif j == 1												% Set matrix values at boundary of y = 0
            G(n,n) = -0.5*((3 * sigma(j,i)) + sigma(j,i-1) +...
                sigma(j,i+1) + sigma(j+1,i));
            G(n,n+1) = 0.5*(sigma(j,i) + sigma(j+1,i));
            G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
            G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
        elseif j == ny												% Set matrix values at boundary of y = 1
            G(n,n) = -0.5*((3 * sigma(j,i)) + sigma(j,i-1) +...
                sigma(j,i+1) + sigma(j-1,i));
            G(n,n-1) = 0.5*(sigma(j,i) + sigma(j-1,i));
            G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
            G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
        else 														% Set matrix values for non boundary elements
            G(n,n) = -0.5*((4 * sigma(j,i)) + sigma(j-1,i) +...
                sigma(j+1,i) + sigma(j,i+1) + sigma(j,i-1));
            G(n,n+1) = 0.5*(sigma(j,i) + sigma(j+1,i));
            G(n,n-1) = 0.5*(sigma(j,i) + sigma(j-1,i));
            G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
            G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
        end
    end
end

V1 = G\B;							% Solve for potentials
V = reshape(V1,[ny,nx]);			% Convert generated vector into a matrix

[Ex1,Ey1] = gradient(V,dx,dy);			% Determine Electric field in the x and y directions from the gradient of the potential
Ex = -Ex1;							% Electric field in the negative gradient of the potential
Ey = -Ey1;

[X,Y] = meshgrid(x,y);

figure(1)							% Generate potential plot
mesh(X,Y,V)
axis([0 L 0 W 0 Vx])
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Potential (V)')
title('Potential Plot')
title(c,'Potential (V)')
print('Part2A','-dpng')

figure(2)							% Generate plot of electrif field as vectors
quiver(X,Y,Ex,Ey)
axis([0 L 0 W])
xlabel('X (m)')
ylabel('Y (m)')
title('Plot of Electric Field Vectors')
print('Part2B','-dpng')