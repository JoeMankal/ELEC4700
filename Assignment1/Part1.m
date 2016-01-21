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

Delta_t = 1.5e-15;                  % Time step
nElectrons = 100;                    % Number of electrons
nSims = 1000;                       % Number of simulations
v_g = sqrt(2 * C.kb * T / m_n);     % Thermal velocity
electrons = zeros(nElectrons,4,nSims);
Temps = zeros(nElectrons,1);
Temp_avg = zeros(nSims,2);

for i = 1:nElectrons
    electrons(i,1,1) = X_length * rand;
    electrons(i,2,1) = Y_length * rand;
    direction = 2*pi*rand;
    electrons(i,3,:) = v_g * cos(direction);
    electrons(i,4,:) = v_g * sin(direction);
    Temps(i) = ((electrons(i,3,1)^2 + electrons(i,4,1)^2) * m_n /...
        (2 * C.kb));
end

Temp_avg(1,2) = mean(Temps);

for i = 2:nSims
    for n = 1:nElectrons
        if electrons(n,1,i-1) >= X_length
            electrons(n,1,i) = (electrons(n,3,i) * Delta_t);
            electrons(n,1,i-1) = NaN;
        elseif electrons(n,1,i-1) <= 0
            electrons(n,1,i) = X_length + (electrons(n,3,i) * Delta_t);
            electrons(n,1,i-1) = NaN;
        else
            electrons(n,1,i) = electrons(n,1,i-1) + (electrons(n,3,i) * Delta_t);
        end
        if electrons(n,2,i-1) >= Y_length || electrons(n,2,i-1) <= 0
            electrons(n,4,:) = electrons(n,4,i) * -1;
        end
        electrons(n,2,i) = electrons(n,2,i-1) + (electrons(n,4,i) * Delta_t);
        subplot(2,1,1)
        plot(squeeze(electrons(n,1,1:i)),squeeze(electrons(n,2,1:i)))
        axis([0 X_length 0 Y_length]);
        xlabel('X (m)')
        ylabel('Y (m)')
        hold on
        Temps(n) =(electrons(n,3,i)^2 + electrons(n,4,i)^2) *...
            m_n / (2 * C.kb);
    end
    Temp_avg(i,2) = mean(Temps);
    Temp_avg(i,1) = Delta_t * (i-1);
    hold off
    subplot(2,1,2)
    plot(Temp_avg(1:i,1),Temp_avg(1:i,2),'b')
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    pause(0.01)
end