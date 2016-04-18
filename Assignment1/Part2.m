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
nElectrons = 100;                       % Number of electrons
nSims = 1000;                           % Number of simulations
v_g = sqrt(C.kb * T / m_n);             % Thermal velocity

electrons = zeros(nElectrons,6);        % Setup matrix for electrons
colours = zeros(nElectrons,3);			% Setup matrix for colours
Temp_avg = zeros(2,1);					% Setup matrix for average temperature of electrons
X = zeros(1,2);							% Setup array of two x-values to plot line segment of electron tregectory
Y = zeros(1,2);							% Setup array of two y-values to plot line segment of electron tregectory
vel_dist = zeros(nElectrons,1);			% Setup array for velocity distribution
P_scat = 1 - exp(-Delta_t / taugh_mn);	% Determine scattering probability

for i = 1:nElectrons
    colours(i,:) = rand(1,3);                                   % Determine random colour for electron
    electrons(i,1) = X_length * rand;   						% Assign electron a random x position
    electrons(i,3) = Y_length * rand;							% Assign electron a random y position
    electrons(i,5) = v_g * randn;                      			% Determine velocity in the x direction at random for MaxB Dist.
    electrons(i,6) = v_g * randn;                      			% Determine velocity in the y direction at random for MaxB Dist.
    vel_dist(i) = sqrt(electrons(i,5)^2 + electrons(i,6)^2);	% Add velocity of electron to array for ditribution
    Temp_avg(2,1) = Temp_avg(2,1) + ((electrons(i,5)^2 +...
        electrons(i,6)^2) * m_n / (2 * C.kb));					% Sum temperature of all electrons
end

figure(1)							% Create Figure1 to display electron tregectories and temperature 
subplot(2,1,1)                      % Setup first subplot for electron tregectories
axis([0 X_length 0 Y_length]);
title('Trajectory of 100 Electrons at 300k')
xlabel('X (m)')
ylabel('Y (m)')
hold on

Temp_avg(2,1) = Temp_avg(2,1)/nElectrons;   % Calculate average temperature
subplot(2,1,2)                              % Setup second subplot for temperature
plot(Temp_avg(1,1),Temp_avg(2,1),'b.')
hold on
axis([0 (Delta_t*nSims) 0 600])
title('Average Temperature of System over Time')
xlabel('Time (s)')
ylabel('Temperature (K)')

for i = 2:nSims                             								% For loop for total number of simulated points
    Temp_avg(2,1) = 0;														% Set average temperature variable to 0 
    subplot(2,1,1)															% Select first subplot
    for n = 1:nElectrons                    								% For loop over all electrons
        
        if P_scat > rand													% Determine if the electron will scatter
            electrons(n,5) = v_g * randn;									% Rethermalize the electron should it scatter
            electrons(n,6) = v_g * randn;
        end
        
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
        
        X = [electrons(n,1) electrons(n,2)];
        Y = [electrons(n,3) electrons(n,4)];
        plot(X,Y,'color',colours(n,:))                                      % Update plot with new trajectory
        Temp_avg(2,1) = Temp_avg(2,1) + ((electrons(n,5)^2 +...             % Sum temperature of electrons throughout loop
            electrons(n,6)^2) * m_n / (2 * C.kb));
        electrons(n,1) = electrons(n,2);
        electrons(n,3) = electrons(n,4);
    end
    Temp_avg(2,1) = Temp_avg(2,1)/nElectrons;                               % Determine average temperature of electrons
    Temp_avg(1,1) = Temp_avg(1,1) + Delta_t;
    subplot(2,1,2)                                                          % Select second subplot
    plot(Temp_avg(1,1),Temp_avg(2,1),'b.')									% Update temperature plot
    pause(1e-9) 															% Pause to allow plot to update
end

print('Part2_1','-dpng')													% Save plot as an image to add into report

figure(2)
histogram(vel_dist)															% Create histogram of velocity distribution
title('Velocity Distribution of Electrons After Initialization')
xlabel('Velocity (m/s)')
ylabel('Count')
print('Part2_2','-dpng')													% Save histogram as an image to add into report