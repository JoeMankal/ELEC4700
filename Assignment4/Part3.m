clear
clc

R = 20;                            
C = 10e-6;                          
Vi = 1;

Delta_t = 1e-6;
t = 0:Delta_t:3e-3;
Vo_n1 = zeros(1,length(t));
Vo_n1(1) = 0; 
for i = 1:(length(t) - 1)
    i_noise = (Vi/(5*R))*randn;
    Vo_n1(i+1) = (Delta_t * (((Vi - Vo_n1(i)) / (R*C)) + (i_noise/C))) + Vo_n1(i);
end
figure(1)
plot(t,Vo_n1)
title('Output of Low Pass Filter with a Current Noise Source')
xlabel('Time (s)')
ylabel('V_o(t) (V)')
grid on
print('Part1_A','-dpng')

Vo_n2 = zeros(1,length(t));
Vo_n2(1) = 0; 
for i = 1:(length(t) - 1)
    i_noise = (1/(5*R))*randn;
    Vo_n2(i+1) = (Delta_t * ((-Vo_n2(i) / (R*C)) + (i_noise/C))) + Vo_n2(i);
end
V_RMS = sqrt(mean(Vo_n2.^2));

Fs = 1/Delta_t;
t = (0:Delta_t:(10 - Delta_t))';
N = size(t,1);
Y = fftshift(fft(Vo_n2,N));
dF = Fs/N;
f = (-Fs/2):dF:(Fs/2-dF);
figure(2)
plot(f,abs(Y)/N)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response with an Input Voltage of 0V')
print('Part1_B','-dpng')

R1 = 10;
C1 = 10e-6;
Vo_n3 = zeros(1,length(t));
Vo_n3(1) = 0; 
for i = 1:(length(t) - 1)
    i_noise = (1/(5*R))*randn;
    Vo_n3(i+1) = (Delta_t * ((-Vo_n3(i) / (R1*C1)) + (i_noise/C1))) + Vo_n3(i);
end
V_RMS2 = sqrt(mean(Vo_n3.^2));
Y2 = fftshift(fft(Vo_n3,N));
figure(3)
plot(f,abs(Y2)/N)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response with R = 10\Omega and C = 10\muF')
print('Part1_C','-dpng')

R2 = 50;
C2 = 10e-6;
Vo_n4 = zeros(1,length(t));
Vo_n4(1) = 0; 
for i = 1:(length(t) - 1)
    i_noise = (1/(5*R))*randn;
    Vo_n4(i+1) = (Delta_t * ((-Vo_n4(i) / (R2*C2)) + (i_noise/C2))) + Vo_n4(i);
end
V_RMS3 = sqrt(mean(Vo_n4.^2));
Y3 = fftshift(fft(Vo_n4,N));
figure(4)
plot(f,abs(Y3)/N)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response with R = 50\Omega and C = 10\muF')
print('Part1_D','-dpng')

R3 = 10;
C3 = 20e-6;
Vo_n5 = zeros(1,length(t));
Vo_n5(1) = 0; 
for i = 1:(length(t) - 1)
    i_noise = (1/(5*R))*randn;
    Vo_n5(i+1) = (Delta_t * ((-Vo_n5(i) / (R3*C3)) + (i_noise/C3))) + Vo_n5(i);
end
V_RMS4 = sqrt(mean(Vo_n5.^2));
Y4 = fftshift(fft(Vo_n5,N));
figure(5)
plot(f,abs(Y4)/N)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response with R = 10\Omega and C = 20\muF')
print('Part1_E','-dpng')