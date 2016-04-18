clear
clc

kb = 1.38064852e-23;
T = 300;
R = 10;
C = 10e-6;
B = 1 / (2*pi*R*C);

v_RMS = sqrt(4*kb*T*R*B);