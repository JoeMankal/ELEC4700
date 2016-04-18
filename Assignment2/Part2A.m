clear
clc

sigma_outside = 1;					% Set conductance for outside of the boxes
sigma_inside = 0.01;				% Set conductance for inside of the boxes
L = 1.5;							% Set length of space in metres (m)
W = 1;								% Set width of space in metres (m)
Lb = L*0.35;						% Set length of box
Wb = W*0.35;						% Set width of box
nx = 150;							% Number of points in x direction
ny = 100;							% Number of points in y direction
x = linspace(0,L,nx);               % Construct linspace for x
y = linspace(0,W,ny);				% Construct linspace for y
sigma = zeros(ny,nx);               % Construct matrix for conductance
G = zeros((nx*ny));                 % Construct G square matrix, nx*ny on a side
Vx = 1;                             % Set boundary condition for x=0
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

[Ex1,Ey1] = gradient(V);			% Determine Electric field in the x and y directions from the gradient of the potential
Ex = -Ex1;							% Electric field in the negative gradient of the potential
Ey = -Ey1;

Jx = sigma.*Ex;						% Electric field times the conductance gives the current density
Jy = sigma.*Ey;
J = sqrt(Jx.^2 + Jy.^2);

I_in = sum(Ex(:,1));				% Current on the left edge of the area is the sum of the current densities for x = 0
I_out = sum(Ex(:,nx));				% Current on the right edge of the are is the sum of the current densities for x = 1
I = (I_in + I_out)/2				% Averaging two current values

[X,Y] = meshgrid(x,y);

figure(1)							% Generate conductance plot
mesh(X,Y,sigma)
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Conductance (S)')
title('Conductance Plot')
title(c,'Conductance (S)')
view([0 90])
print('Part2A_1','-dpng')

figure(2)							% Generate potential plot
mesh(X,Y,V)
axis([0 L 0 W 0 Vx])
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Potential (V)')
title('Potential Plot')
title(c,'Potential (V)')
view([0 90])
print('Part2A_2','-dpng')

figure(3)							% Generate plot of electric field in the x-direction
mesh(X,Y,Ex)
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Electric Field in the x-direction (V/m)')
title('Plot of Electric Field in the x-direction')
title(c,'Electric Field Strength (V/m)')
view([0 90])
print('Part2A_3','-dpng')

figure(4)							% Generate plot of electric field in the direction
mesh(X,Y,Ey)
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Electric Field in the y-direction (V/m)')
title('Plot of Electric Field in the y-direction')
title(c,'Electric Field Strength (V/m)')
view([0 90])
print('Part2A_4','-dpng')

figure(5)							% Generate plot of electrif field as vectors
quiver(X,Y,Ex,Ey)
axis([0 L 0 W])
xlabel('X (m)')
ylabel('Y (m)')
title('Plot of Electric Field Vectors')
print('Part2A_5','-dpng')

figure(6)							% Generate plot of current density
mesh(X,Y,J)
c = colorbar;
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Current Density (A/m)')
title('Plot of Current Density')
ylabel(c,'Current Density (A/m)')
view([0 90])
print('Part2A_6','-dpng')