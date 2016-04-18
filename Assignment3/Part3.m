clear
clc

C.q_0 = 1.60217653e-19;                 % electron charge
C.m_0 = 9.10938215e-31;                 % electron mass
C.kb = 1.3806504e-23;                   % Boltzmann constant

m_n = 0.26 * C.m_0;                     % Effective mass of electrons
X_length = 200e-9;                      % Size of x region in m 
Y_length = 100e-9;                      % Size of y region in m
T = 300;                                % Tempurature in Kelvin
taugh_mn = 0.2e-12;                     % Mean time between collisions
Delta_t = 1e-15;                        % Time step
nElectrons = 1000;                      % Number of electrons
nSims = 2000;                           % Number of simulations
v_g = sqrt(2 * C.kb * T / m_n);         % Thermal velocity
electrons = zeros(nElectrons,6);        % Setup matrix for electrons
colours = zeros(nElectrons,3);          % Setup matrix for colours
X = zeros(1,2);                         % Setup array of two x-values to plot line segment of electron tregectory
Y = zeros(1,2);                         % Setup array of two y-values to plot line segment of electron tregectory
P_scat = 1 - exp(-Delta_t / taugh_mn);
Step_Size = 5e-9;                       % Set size of boxes for density and termperature map
Density = zeros(int64((Y_length/Step_Size)),int64((X_length/Step_Size)));	% Create matrix to store density values
sigma_outside = 1;                      % Set conductance for outside of the boxes
sigma_inside = 0.01;                    % Set conductance for inside of the boxes
Lb = X_length*0.35;                     % Set length of box
Wb = Y_length*0.35;                     % Set width of box
nx = 150;                               % Number of points in x direction
ny = 100;                               % Number of points in y direction
x = linspace(0,X_length,nx);            % Construct linspace for x
y = linspace(0,Y_length,ny);            % Construct linspace for y
dx = x(2) - x(1);
dy = y(2) - y(1);
sigma = zeros(ny,nx);                   % Construct matrix for conductance
G = zeros((nx*ny));                     % Construct G square matrix, nx*ny on a side
Vx = 0.1;                               % Set boundary condition for x=0
B = zeros((nx*ny),1);                   % Construct B matrix
Box1_specular = true;					% Set whether the bottom box is specular or diffuse
Box2_specular = true;					% Set whether the top box is specular or diffuse
Correct_dir = false;					% Create boolean to determine if re-thermalized electrons have the correct direction
InBox = true;							% Create boolean to determine if newly created electrons are within one of the boxes
x1 = 0.8e-7;							% Set first x boundary of both boxes
x2 = 1.2e-7;							% Set second x boundary of both boxes
y11 = 0;								% Set first y boundary of the bottom box
y12 = 0.4e-7;							% Set second y boundary for the bottom box
y21 = 0.6e-7;							% Set first y boundary for the top box
y22 = 1e-7;								% Set second y boundary for the top box
Box_X = [x1 x2 x2 x1 x1];				% Set array for x values of corners of both boxes
Box_Y1 = [y11 y11 y12 y12 y11];			% Set array for y values of corners of bottom box
Box_Y2 = [y21 y21 y22 y22 y21];			% Set array for y values of corners of top box

for i = 1:nx                            % Set conductance for each point whether it is within a box or not
    for j = 1:ny
        if (y(j) <= Wb) && (x(i) >= Lb) && (x(i) <= X_length-Lb) ||...
                (y(j) >= Y_length - Wb) && (x(i) >= Lb) &&...
                (x(i) <= X_length - Lb)
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

V1 = G\B;                               % Solve for potentials
V = reshape(V1,[ny,nx]);                % Convert generated vector into a matrix
[Ex1,Ey1] = gradient(V,dx,dy);          % Determine Electric field in the x and y directions from the gradient of the potential
Ex = -Ex1;                              % Electric field in the negative gradient of the potential
Ey = -Ey1;
Fx = C.q_0 * Ex;                        % Force felt by electrons in the x direction
Fy = C.q_0 * Ey;                        % Force on electron in the y-direction
a_x = Fx / m_n;                         % Acceleration in x-direction
a_y = Fy / m_n;                         % Acceleration in y-direction

for i = 1:nElectrons
    colours(i,:) = rand(1,3);                                               % Determine random colour for electron
    electrons(i,1) = X_length * rand;                                       % Assign electron a random x position
    electrons(i,3) = Y_length * rand;                                       % Assign electron a random y position
    InBox = true;                                                           % Set boolean 'InBox' to true to be able to check 
	while InBox == true                                                     % if the elctron is inside one of the boxes
        if (electrons(i,3) <= y12 || electrons(i,3) >= y21) &&...           % Determines if the electron is within one of the boxes
                (electrons(i,1) >= x1 && electrons(i,1) <= x2)	
            electrons(i,1) = X_length * rand;                               % Re-initialize x position
            electrons(i,3) = Y_length * rand;                               % Re-initialize y position
        else
            InBox = false;                                                  % Set to false if electron is not in a box			
        end
    end
    electrons(i,5) = v_g * randn;                                           % Determine velocity in the x direction at random for MaxB Dist.
    electrons(i,6) = v_g * randn;                                           % Determine velocity in the y direction at random for MaxB Dist.
end

figure(1)                                                       % Create Figure1 to display electron tregectories and temperature 
plot(Box_X,Box_Y1,'k-',Box_X,Box_Y2,'k-')
axis([0 X_length 0 Y_length]);
title('Trajectory of 15 Electrons at 300k')
xlabel('X (m)')
ylabel('Y (m)')
hold on

for i = 2:nSims                             % For loop for total number of simulated points
    for n = 1:nElectrons                    % For loop over all electrons
        
        Correct_dir = false;
        
        if P_scat > rand													% Determine if the electron will scatter
                electrons(n,5) = v_g * randn;								% Rethermalize the electron should it scatter
                electrons(n,6) = v_g * randn;
        end
        
        acc_x_index = round(electrons(n,1)/dx);
        acc_y_index = round(electrons(n,3)/dy);
        if acc_x_index < 1
            acc_x_index = 1;
        end
        if acc_y_index < 1
            acc_y_index = 1;
        end
        
        electrons(n,5) = electrons(n,5) + (a_x(acc_y_index,acc_x_index) * Delta_t);
        electrons(n,6) = electrons(n,6) + (a_y(acc_y_index,acc_x_index) * Delta_t);
        electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);		% Update x location of electron
        electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);		% Update y location of electron
        
		previous_x = electrons(n,1) - (electrons(n,5) * Delta_t);			% Store previous x location of electron
		previous_y = electrons(n,3) - (electrons(n,6) * Delta_t);			% Store previous y location of electron
		
        if electrons(n,2) >= x1 && electrons(n,2) <= x2 &&...				% Determine if an electron collides with the bottom box when it is specular
                electrons(n,4) <= y12 && Box1_specular == true
            if previous_x < x1                                              % Adjust electron tregectory accordingly
                electrons(n,5) = electrons(n,5) * -1;
                electrons(n,2) = x1;
            elseif previous_x > x2
                electrons(n,5) = electrons(n,5) * -1;
                electrons(n,2) = x2;
            elseif previous_y > y12
                electrons(n,6) = electrons(n,6) * -1;
                electrons(n,4) = y12;
            end
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);	% Update x location of electron
            electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);	% Update y location of electron
        elseif electrons(n,2) > x1 && electrons(n,2) < x2 &&...				% Determine if an electron collides with the bottom box when it is diffuse
                electrons(n,4) < y12 && Box1_specular == false			
            electrons(n,5) = v_g * randn;									% Rethermalize electron
            electrons(n,6) = v_g * randn;
            if previous_x < x1												% Adjust electron location to edge of box to ensure electrons trajectory,
                electrons(n,1) = x1;										% is outside of the box on next iteration
            elseif previous_x > x2
                electrons(n,1) = x2;
            elseif previous_y > y12
                electrons(n,3) = y12;
            end
            while Correct_dir == false				
                if ((previous_x < x1) && (electrons(n,5) >= 0)) ||...		% Determine if the new direction of the re-thermalized electron is correct,
                        ((previous_x > x2) && (electrons(n,5) <= 0)) ||...	% that is, the electron will not hit into the box on the next iteration
                        ((previous_y > y12) && (electrons(n,6) <= 0))
                    electrons(n,5) = v_g * randn;							% Rethermalize the electron if the direction is not correct
                    electrons(n,6) = v_g * randn;
                else
                    Correct_dir = true;
                end
            end
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);	% Update x location of electron
            electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);	% Update y location of electron
        elseif electrons(n,2) >= x1 && electrons(n,2) <= x2 &&...           % Determine if an electron collides with the top box when it is specular
                electrons(n,4) >= y21 && Box2_specular == true
            if previous_x < x1                                              % Adjust electron tregectory accordingly
                electrons(n,5) = electrons(n,5) * -1;
                electrons(n,1) = x1;
            elseif previous_x > x2
                electrons(n,1) = electrons(n,5) * -1;
                electrons(n,2) = x2;
            elseif previous_y < y21
                electrons(n,6) = electrons(n,6) * -1;
                electrons(n,3) = y21;
            end
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);	% Update x location of electron
            electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);	% Update y location of electron
        elseif electrons(n,2) >= x1 && electrons(n,2) <= x2 &&...			% Determine if an electron collides with the top box when it is diffuse
                electrons(n,4) >= y21 && Box2_specular == false
            electrons(n,5) = v_g * randn;									% Rethermalize electron
            electrons(n,6) = v_g * randn;
            if previous_x < x1												% Adjust electron location to edge of box to ensure electrons trajectory,
                electrons(n,1) = x1;										% is outside of the box on next iteration
            elseif previous_x > x2
                electrons(n,1) = x2;
            elseif previous_y < y21
                electrons(n,3) = y21;
            end
            while Correct_dir == false
                if ((previous_x < x1) && (electrons(n,5) >= 0)) ||...		% Determine if the new direction of the re-thermalized electron is correct,
                        ((previous_x > x2) && (electrons(n,5) <= 0)) ||...	% that is, the electron will not hit into the box on the next iteration
                        ((previous_y < y21) && (electrons(n,6) >= 0))
                    electrons(n,5) = v_g * randn;							% Rethermalize the electron if the direction is not correct
                    electrons(n,6) = v_g * randn;
                else
                    Correct_dir = true;
                end
            end
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);	% Update x location of electron
            electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);	% Update y location of electron
        end
            
        if electrons(n,2) > X_length										% Determine if electron will exceed right x boundary
            electrons(n,1) = 0;												% Adjust electron trajectory if boundary is passed
            electrons(n,2) = electrons(n,5) * Delta_t;	
        elseif electrons(n,2) < 0											% Determine if electron will exceed left x boundary
            electrons(n,1) = X_length;										% Adjust electron trajectory if boundary is passed
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);	
        end
        
        if electrons(n,4) >= Y_length || electrons(n,4) <= 0				% Determine if electron exceeds top or bottom y boundary
            electrons(n,6) = -1 * electrons(n,6);							% Adjust electron trajectory if boundary is passed
        end
        
        if n < 15
            X = [electrons(n,1) electrons(n,2)];
            Y = [electrons(n,3) electrons(n,4)];
            plot(X,Y,'color',colours(n,:))                                      % Update plot with new trajectory
        end
        electrons(n,1) = electrons(n,2);
        electrons(n,3) = electrons(n,4);
    end
    %pause(1e-19)                                            				% Pause to allow plot to update
end
print('Part3_A','-dpng')														% Save plot as an image to add into report

for n = 1:(X_length/Step_Size)												% Loop over defined step size in x for the density/tempurature map
    for m = 1:(Y_length/Step_Size)											% Loop over defined step size in y for the density/tempurature map
        for i = 1:nElectrons												% Loop over all electrons
            if (electrons(i,1) > (Step_Size*(n-1))) &&...					% Determine if the electron is within the defined square
                    (electrons(i,1) < (Step_Size*n)) &&...
                    (electrons(i,3) > (Step_Size*(m-1))) &&...
                    (electrons(i,3) < (Step_Size*m))
                Density(m,n) = Density(m,n) + 1;							% Increment density
            end
        end
    end
end

figure(2)
imagesc(Density)															% Create map for the desnity of electrons
title('Density Map of Electrons')
xlabel('X (5X10^{-9} m)')
ylabel('Y (5X10^{-9} m)')
set(gca,'Ydir','Normal')
c = colorbar;
title(c,'Number of Electrons')
print('Part3_B','-dpng')													% Save density map as an image to add into report