clear
clc

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres per s²

m_n = 0.26 * C.m_0;                 % Effective mass of electrons
X_length = 200e-9;                  % Size of x region in m 
Y_length = 100e-9;                  % Size of y region in m
T = 300;                            % Tempurature in Kelvin
taugh_mn = 0.2e-12;                 % Mean time between collisions

Delta_t = 1.5e-15;                  % Time step
nElectrons = 10;                    % Number of electrons
nSims = 1000;                       % Number of simulations
v_g = sqrt(C.kb * T / m_n);         % Thermal velocity
electrons = zeros(nElectrons,6);
Temps = 0;
Temp_avg = zeros(nSims,2);
vel_dist = zeros(nElectrons,1);
P_scat = 1 - exp(-Delta_t / taugh_mn);

for i = 1:nElectrons
    [electrons(i,1), electrons(i,2)] = deal(X_length * rand);
    [electrons(i,3), electrons(i,4)] = deal(Y_length * rand);
    electrons(i,5) = v_g * randn;
    electrons(i,6) = v_g * randn;
    vel_dist(i) = sqrt(electrons(i,5)^2 + electrons(i,6)^2);
    Temps = Temps + (electrons(i,5)^2 + electrons(i,6)^2) * m_n / (2 * C.kb);
end

Temp_avg(1,2) = Temps/nElectrons; 

% figure
% histogram(vel_dist,50)
% mean(vel_dist);

for i = 2:nSims
    tic
    for n = 1:nElectrons
        %Temps = 0;
        if P_scat > rand
            electrons(n,5) = v_g * randn;
            electrons(n,6) = v_g * randn;
        end
        
        if electrons(n,1) >= X_length
            electrons(n,1) = 0;
            electrons(n,2) = electrons(n,5) * Delta_t;
        elseif electrons(n,1) <= 0
            electrons(n,1) = X_length;
            electrons(n,2) = X_length + (electrons(n,5) * Delta_t);
        else
            electrons(n,2) = electrons(n,1) + (electrons(n,5) * Delta_t);
        end
        
        if electrons(n,3) >= Y_length || electrons(n,3) <= 0
            electrons(n,6) = -1 * electrons(n,6);
        end
        
        electrons(n,4) = electrons(n,3) + (electrons(n,6) * Delta_t);
        %subplot(2,1,1)
        X = [electrons(n,1) electrons(n,2)];
        Y = [electrons(n,3) electrons(n,4)];
        plot(X,Y)
        axis([0 X_length 0 Y_length]);
        xlabel('X (m)')
        ylabel('Y (m)')
        hold on
        %Temps = Temps + (electrons(n,5)^2 + electrons(n,6)^2) * m_n / (2 * C.kb);
        electrons(n,1) = electrons(n,2);
        electrons(n,3) = electrons(n,4);
    end
    toc
    %Temp_avg(i,2) = Temps/nElectrons;
    %Temp_avg(i,1) = Delta_t * (i-1);
%     hold off
%     subplot(2,1,2)
%     plot(Temp_avg(1:i,1),Temp_avg(1:i,2),'b')
%     axis([0 (Delta_t*nSims) 0 inf])
%     xlabel('Time (s)')
%     ylabel('Temperature (K)')
    pause(1e-9)
end

plot(Temp_avg(1:i,1),Temp_avg(1:i,2),'b')