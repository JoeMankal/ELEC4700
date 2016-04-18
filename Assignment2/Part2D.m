clear
clc

sigma_outside = 1;							% Conductance for the area outside of the boxes
L = 1.5;									% Set length of space in metres (m)
W = 1;										% Set width of space in metres (m)
Lb = L*0.35;								% Set length of boxes
Wb = W*0.35;								% Set width of boxes
Vx = 1;                                 	% Set boundary condition for x=0
nx = 150;									% Set number of points in the x-direction
ny = 100;									% Set number of points in the y-direction
sigma_var = 0.001:0.009:1;                  % Set varying conductances 
I = zeros(length(sigma_var),1);             % Create matrix for currents of different conductance values 
x = linspace(0,L,nx);   		            % Construct linspace for x
y = linspace(0,W,ny);						% Construct linspace for y

for l = 1:length(sigma_var)  				% Loop for different conductance values of the boxes
    sigma_inside = sigma_var(l);     		% Set conductance value of the boxes
    sigma = zeros(ny,nx);               	% Construct m,atrix for conductance
    G = zeros((nx*ny));                 	% Construct G square matrix, nx*ny on a side                      
    B = zeros((nx*ny),1);					% Construct B matrix

    for i = 1:nx							% Set conductance for each point whether it is within the box or not
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
            n = j + (i-1)*ny;										% Determine the index for diagonal
            if i == 1												% Set boundary condition for x = 0
                G(n,n) = 1;
                B(n) = Vx;
            elseif i == nx											% Set boundary condition for x = 1
                G(n,n) = 1;
            elseif j == 1											% Set matrix values at boundary of y = 0
                G(n,n) = -0.5*((3 * sigma(j,i)) + sigma(j,i-1) +...
                    sigma(j,i+1) + sigma(j+1,i));
                G(n,n+1) = 0.5*(sigma(j,i) + sigma(j+1,i));
                G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
                G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
            elseif j == ny											% Set matrix values at boundary of y = 1
                G(n,n) = -0.5*((3 * sigma(j,i)) + sigma(j,i-1) +...
                    sigma(j,i+1) + sigma(j-1,i));
                G(n,n-1) = 0.5*(sigma(j,i) + sigma(j-1,i));
                G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
                G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
            else 													% Set matrix values for non boundary elements 
                G(n,n) = -0.5*((4 * sigma(j,i)) + sigma(j-1,i) +...
                    sigma(j+1,i) + sigma(j,i+1) + sigma(j,i-1));
                G(n,n+1) = 0.5*(sigma(j,i) + sigma(j+1,i));
                G(n,n-1) = 0.5*(sigma(j,i) + sigma(j-1,i));
                G(n,n-ny) = 0.5*(sigma(j,i) + sigma(j,i-1));
                G(n,n+ny) = 0.5*(sigma(j,i) + sigma(j,i+1));
            end
        end
    end
    V1 = G\B;									% Solve for potentials
    V = reshape(V1,[ny,nx]);					% Convert generated vectors into a matrix

    [Ex1,Ey1] = gradient(V);					% Gradient of the potential to determine electric field
    Ex = -Ex1;									% Negate values since the E field is the negative gradient of V
    Ey = -Ey1;

    Jx = sigma.*Ex;								% Electric field times the conductance gives the current density
    Jy = sigma.*Ey;

    I_in = sum(Ex(:,1));						% Current on the left edge of the area is the sum of the current density at x = 0
    I_out = sum(Ex(:,nx));						% Current on the right edge of the area is the sum of the current desnity at x = 1
    I(l) = (I_in + I_out)/2;            		% Averaging the two current values
end

figure(1)										% Generate plot for current against conductance of the boxes
plot(sigma_var,I)
grid on
xlabel('Conductance of box (S)')
ylabel('Current (A)')
title('Current vs. Conductance of Box')
print('Part2D','-dpng')