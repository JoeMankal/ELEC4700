clear
clc

R = 20;                             % Set resistance of resistance in ohms
C = 10e-6;                          % Set capacitance of capacitor in farads
Vi = 1;

Transfer_func = @(w) abs((1 + (2*pi*i*w*R*C)).^(-1));

f = logspace(0,5,1000000);
Gain = Transfer_func(f);

figure(1)
loglog(f,Gain)
xlabel('Frequency (Hz)')
ylabel('log(|H(jw)|)')
title('Frequency Response of a Low-Pass Filter')
grid on
print('Part1_A','-dpng')

Vo = @(t) Vi * (1 - exp(-t/(R*C)));
Delta_t1 = 1e-6;
Delta_t2 = 1e-5;
Delta_t3 = 1e-4;
t1 = 0:Delta_t1:2e-3;
t2 = 0:Delta_t2:2e-3;
t3 = 0:Delta_t3:2e-3;
Vo_n1 = zeros(1,length(t1));
Vo_n2 = zeros(1,length(t2));
Vo_n3 = zeros(1,length(t3));
[Vo_n1(1),Vo_n2(1),Vo_n3(1)] = deal(0); 
for i = 1:(length(t1) - 1)
    Vo_n1(i+1) = (Delta_t1 * (Vi - Vo_n1(i)) / (R*C)) + Vo_n1(i);
end
for i = 1:(length(t2) - 1)
    Vo_n2(i+1) = (Delta_t2 * (Vi - Vo_n2(i)) / (R*C)) + Vo_n2(i);
end
for i = 1:(length(t3) - 1)
    Vo_n3(i+1) = (Delta_t3 * (Vi - Vo_n3(i)) / (R*C)) + Vo_n3(i);
end
figure(2)
fplot(Vo,[0 2e-3])
title('Analytical Solution for Output of Low Pass Filter')
xlabel('Time (s)')
ylabel('V_o(t) (V)')
grid on
print('Part1_B','-dpng')
figure(3)
plot(t1,Vo_n1,t2,Vo_n2,t3,Vo_n3)
title('Discrete Solution for Output of Low Pass Filter')
legend('\Deltat = 1 \mus','\Deltat = 10 \mus','\Deltat = 100 \mus')
xlabel('Time (s)')
ylabel('V_o(t) (V)')
grid on
print('Part1_C','-dpng')

freq_1 = 1e3;
freq_2 = 1e4;
freq_3 = 1e5;
E_1 = @(t) sin(2*pi*freq_1*t);
E_2 = @(t) sin(2*pi*freq_2*t);
E_3 = @(t) sin(2*pi*freq_3*t);
Delta_t1 = 1/(100*freq_1);
Delta_t2 = 1/(100*freq_2);
Delta_t3 = 1/(100*freq_3);
t1 = 0:Delta_t1:5e-3;
t2 = 0:Delta_t2:5e-3;
t3 = 0:Delta_t3:5e-3;
E1 = E_1(t1);
E2 = E_2(t2);
E3 = E_3(t3);
Vo2_n1 = zeros(1,length(t1));
Vo2_n2 = zeros(1,length(t2));
Vo2_n3 = zeros(1,length(t3));
[Vo2_n1(1),Vo2_n2(1),Vo2_n3(1)] = deal(0); 
for i = 1:(length(t1) - 1)
    Vo2_n1(i+1) = (Delta_t1 * (E1(i) - Vo2_n1(i)) / (R*C)) + Vo2_n1(i);
end
for i = 1:(length(t2) - 1)
    Vo2_n2(i+1) = (Delta_t2 * (E2(i) - Vo2_n2(i)) / (R*C)) + Vo2_n2(i);
end
for i = 1:(length(t3) - 1)
    Vo2_n3(i+1) = (Delta_t3 * (E3(i) - Vo2_n3(i)) / (R*C)) + Vo2_n3(i);
end

figure(4)
plot(t1,Vo2_n1,t2,Vo2_n2,t3,Vo2_n3)
title('Response of Low Pass Filter for Various Input Frequncies')
legend('f = 1 kHz','f = 10 kHz','f = 100 kHz')
xlabel('Time (s)')
ylabel('V_o(t) (V)')
grid on
print('Part1_D','-dpng')

freq1 = 1e3;
Fs1 = 5*freq1;
Delta_t4 = 1/Fs1;
t4 = (0:Delta_t4:(1000 - Delta_t4))';
N1 = size(t4,1);
Vi1 = @(t) sin(2*pi*freq1*t);
E4 = Vi1(t4);
Vout1 = zeros(1,length(t4));
Vout1(1) = 0;

for i = 1:(length(t4) - 1)
    Vout1(i+1) = (Delta_t4 * (E4(i) - Vout1(i)) / (R*C)) + Vout1(i);
end

Y = fftshift(fft(Vout1,N1));
dF1 = Fs1/N1;
f1 = (-Fs1/2):dF1:((Fs1/2)-dF1);
figure(5)
plot(f1,abs(Y)/N1)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response of sin(2\pi10^{3}t)')
print('Part1_E','-dpng')

freq2 = 1e4;
Fs2 = 4*freq2;
Delta_t5 = 1/Fs2;
t5 = (0:Delta_t5:(100 - Delta_t5))';
N2 = size(t5,1);
Vi2 = @(t) sin(2*pi*freq2*t);
E5 = Vi2(t5);
Vout2 = zeros(1,length(t5));
Vout2(1) = 0;

for i = 1:(length(t5) - 1)
    Vout2(i+1) = (Delta_t5 * (E5(i) - Vout2(i)) / (R*C)) + Vout2(i);
end

Y2 = fftshift(fft(Vout2,N2));
dF2 = Fs2/N2;
f2 = (-Fs2/2):dF2:(Fs2/2-dF2);
figure(6)
plot(f2,abs(Y2)/N2)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response of sin(2\pi10^{4}t)')
print('Part1_F','-dpng')

freq3 = 1e6;
Fs3 = 4*freq3;
Delta_t6 = 1/Fs3;
t6 = (0:Delta_t6:(1 - Delta_t6))';
N3 = size(t6,1);
Vi3 = @(t) sin(2*pi*freq3*t);
E6 = Vi3(t6);
Vout3 = zeros(1,length(t6));
Vout3(1) = 0;

for i = 1:(length(t6) - 1)
    Vout3(i+1) = (Delta_t6 * (E6(i) - Vout3(i)) / (R*C)) + Vout3(i);
end

Y3 = fftshift(fft(Vout3,N3));
dF3 = Fs3/N3;
f3 = (-Fs3/2):dF3:(Fs3/2-dF3);
figure(7)
plot(f3,abs(Y3)/N3)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Response of sin(2\pi10^{6}t)')
print('Part1_G','-dpng')