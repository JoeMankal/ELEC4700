clear
clc

C.q_0 = 1.60217653e-19;                 % electron charge
C.hb = 1.054571596e-34;                 % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;                 % electron mass
C.kb = 1.3806504e-23;                   % Boltzmann constant
C.eps_0 = 8.854187817e-12;              % vacuum permittivity
C.mu_0 = 1.2566370614e-6;               % vacuum permeability
C.c = 299792458;                        % speed of light
C.g = 9.80665;                          % metres per s²

m_n = 0.26 * C.m_0;                     % Effective mass of electrons
X_length = 200e-9;                      % Size of x region in m 
Y_length = 100e-9;                      % Size of y region in m
T = 300;                                % Tempurature in Kelvin
taugh_mn = 0.2e-12;                     % Mean time between collisions
Delta_t = 1e-15;                        % Time step
nElectrons = 30000;                     % Number of electrons
nSims = 2000;                           % Number of simulations
v_g = sqrt(2 * C.kb * T / m_n);         % Thermal velocity
electrons = zeros(nElectrons,6);        % Setup matrix for electrons
colours = zeros(nElectrons,3);          % Setup matrix for colours
Temp_avg = zeros(2,1);                  % Setup matrix for average temperature of electrons
X = zeros(1,2);                         % Setup array of two x-values to plot line segment of electron tregectory
Y = zeros(1,2);                         % Setup array of two y-values to plot line segment of electron tregectory
E_xfield = 0.1 / X_length;              % Constant electric field in the x direction
E_yfield = 0;                           % Electric field in the y-direction
Fx = C.q_0 * E_xfield;                  % Force felt by electrons in the x direction
Fy = C.q_0 * E_yfield;                  % Force on electron in the y-direction
a_x = Fx / m_n;                         % Acceleration in x-direction
a_y = Fy / m_n;                         % Acceleration in y-direction
P_scat = 1 - exp(-Delta_t / taugh_mn);	% Determine scattering probability
e_con = (1e15) * (100^2);               % Electron concentration in m^-2
I_x = zeros(1,nSims);
Step_Size = 5e-9;                       % Set size of boxes for density and termperature map
Density = zeros(int64((Y_length/Step_Size)),int64((X_length/Step_Size)));	% Create matrix to store density values
Temps = zeros(int64((Y_length/Step_Size)),int64((X_length/Step_Size)));		% Create matrix to store tempurature values

for i = 1:nElectrons
    colours(i,:) = rand(1,3);                                   % Determine random colour for electron
    electrons(i,1) = X_length * rand;                           % Assign electron a random x position
    electrons(i,3) = Y_length * rand;                           % Assign electron a random y position
    direction = 2*pi*rand;                                      % Determine random angle for direction of electron
    electrons(i,5) = v_g * cos(direction);                      % Determine corresponding velocity in the x direction
    electrons(i,6) = v_g * sin(direction);                      % Determine corresponding velocity in the y direction
    Temp_avg(2) = Temp_avg(2) + ((electrons(i,5)^2 +...
        electrons(i,6)^2) * m_n / (2 * C.kb));					% Sum temperature of all electrons
end

V_dx = mean(electrons(:,5));                % Calculate average drift velocity of electrons for initial setup
I_x(1) = Y_length * X_length * C.q_0 * e_con * V_dx;  

figure(1)                                   % Create Figure1 to display electron tregectories and temperature 
subplot(2,1,1)                              % Setup first subplot for electron tregectories
axis([0 X_length 0 Y_length]);
title('Trajectory of 15 Electrons at 300k')
xlabel('X (m)')
ylabel('Y (m)')
hold on

Temp_avg(2) = Temp_avg(2)/nElectrons;       % Calculate average temperature
subplot(2,1,2)                              % Setup second subplot for temperature
plot(Temp_avg(1),Temp_avg(2),'b.')
hold on
axis([0 (Delta_t*nSims) 0 1000])
title('Average Temperature of System over Time')
xlabel('Time (s)')
ylabel('Temperature (K)')

for i = 2:nSims                             % For loop for total number of simulated points
    Temp_avg(2) = 0;						% Set average temperature variable to 0 
    subplot(2,1,1)							% Select first subplot
    for n = 1:nElectrons                    % For loop over all electrons
        
        if P_scat > rand													% Determine if the electron will scatter
            electrons(n,5) = v_g * randn;									% Rethermalize the electron should it scatter
            electrons(n,6) = v_g * randn;
        end
        
        electrons(n,5) = electrons(n,5) + (a_x * Delta_t);
        electrons(n,6) = electrons(n,6) + (a_y * Delta_t);
        electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);       % Update x location of electron
        electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);       % Update y location of electron
        
        if electrons(n,2) > X_length                                        % Determine if electron will exceed right x boundary
            electrons(n,1) = 0;                                             % Adjust electron trajectory if boundary is passed
            electrons(n,2) = electrons(n,5) * Delta_t;
        elseif electrons(n,2) < 0                                           % Determine if electron will exceed left x boundary
            electrons(n,1) = X_length;                                      % Adjust electron trajectory if boundary is passed
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);
        end
        
        if electrons(n,4) >= Y_length || electrons(n,4) <= 0                % Determine if electron exceeds top or bottom y boundary
            electrons(n,6) = -1 * electrons(n,6);                           % Adjust electron trajectory if boundary is passed
        end
        
        if n < 15 
            X = [electrons(n,1) electrons(n,2)];
            Y = [electrons(n,3) electrons(n,4)];
            plot(X,Y,'color',colours(n,:))                                  % Update plot with new trajectory
        end
        Temp_avg(2) = Temp_avg(2) + ((electrons(n,5)^2 +...                 % Sum temperature of electrons throughout loop
            electrons(n,6)^2) * m_n / (2 * C.kb));
        electrons(n,1) = electrons(n,2);
        electrons(n,3) = electrons(n,4);
    end
    V_dx = mean(electrons(:,5));                                            % Calculate average drift velocity of electrons
    I_x(i) = X_length * C.q_0 * e_con * V_dx;
    Temp_avg(2) = Temp_avg(2,1)/nElectrons;                                 % Determine average temperature of electrons
    Temp_avg(1) = Temp_avg(1,1) + Delta_t;
    subplot(2,1,2)                                                          % Select second subplot
    plot(Temp_avg(1),Temp_avg(2),'b.')  									% Update temperature plot
    %pause(1e-19)                                            				% Pause to allow plot to update
end
print('Part1_A','-dpng')                                                    % Save plot as an image to add into report

Time = (0:1:nSims-1) * Delta_t;
figure(2)                                                                   % Create Figure2 to display electron tregectories and temperature 
plot(Time,I_x)
title('Current')
xlabel('Time (s)')
ylabel('Current (A)')
print('Part1_B','-dpng')

for n = 1:(X_length/Step_Size)												% Loop over defined step size in x for the density/tempurature map
    for m = 1:(Y_length/Step_Size)											% Loop over defined step size in y for the density/tempurature map
        for i = 1:nElectrons												% Loop over all electrons
            if (electrons(i,1) > (Step_Size*(n-1))) &&...					% Determine if the electron is within the defined square
                    (electrons(i,1) < (Step_Size*n)) &&...
                    (electrons(i,3) > (Step_Size*(m-1))) &&...
                    (electrons(i,3) < (Step_Size*m))
                Density(m,n) = Density(m,n) + 1;							% Increment density
                Temps(m,n) = Temps(m,n) + ((electrons(i,5)^2 +...			% Add temperature to sum
                    electrons(i,6)^2) * m_n / (2 * C.kb));
            end
        end
        Temps(m,n) = Temps(m,n) / Density(m,n);								% Take average of the temperature
    end
end

figure(3)
imagesc(Density)															% Create map for the desnity of electrons
title('Density Map of Electrons')
xlabel('X (5X10^{-9} m)')
ylabel('Y (5X10^{-9} m)')
set(gca,'Ydir','Normal')
c = colorbar;
title(c,'Number of Electrons')
print('Part1_C','-dpng')													% Save density map as an image to add into report

figure(4)
imagesc(Temps)																% Create map for the temperature of electrons
title('Temperature Map of Electrons')
xlabel('X (5X10^{-9} m)')
ylabel('Y (5X10^{-9} m)')
set(gca,'Ydir','Normal')
c = colorbar;
title(c,'Temperature (K)')
print('Part1_D','-dpng')													% Save temperature map as an image to add into report