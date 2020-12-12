clear;
close all;
clc;

%% data 
k_p = 1.41e5; %N/m
m_p = 0.128*0.385; %kg
A_p = 0.0375*0.385; %m^2
r_p = 32; %Nm/kg/s

m_h = 0.000804; %kg
A_h = 0.00785; %m^2
R_h = 30; %N/m

V = 0.0172; % m^3

c = 343; %[m/s]
rho = 1.2; %[kg/m^3]
%% question 1

omega = 2*pi*(0:0.1:500); % [rad/s] %cambiato questo in rad/s per usarlo dopo nella forma analitica

M_p = m_p/A_p^2;
M_h = m_h/A_h^2;
C_p = A_p^2/k_p;
C_v = V/rho/c^2;

%Analytic bridge impedance
Zb1=1i*omega.*M_h+R_h;
Zb2=(1./(1i*(C_v*omega)));
Zeq=((1./Zb1 + 1./Zb2).^-1);
Zb=Zeq+r_p+1i*(M_p*omega-(1./(C_p*omega)));

figure (1)
subplot 211
plot (omega/2/pi, db(abs(Zb)));
xlabel('Frequency [Hz]');
ylabel('|Z_b| [kg/s]');
grid on;
xlim([0 500]);

subplot 212
plot (omega/2/pi, angle(Zb)*180/pi);
xlabel('Frequency [Hz]');
ylabel('[deg]')
grid on;
xlim([0 500]);

%% question 2
L0 = 0.64516;
Fs = 44100;% [Hz]
Nfft = Fs;% Number of fft points
LA1 = 220;% [Hz]

Ndel = round(Fs/LA1); %Number of delay samples for delay line
omega = linspace(0, 2*pi*(Fs/2),Nfft);

Zb1=1i*omega.*M_h+R_h;
Zb2=(1./(1i*(C_v*omega)));
Zeq=((1./Zb1 + 1./Zb2).^-1);
Zb=Zeq+r_p+1i*(M_p*omega-(1./(C_p*omega)));
Z = Zb./(max(abs(Zb(50:length(Zb))))); %perchè z a 0 + infinito


                        %excitation point
% schema L1-L2 Nut ========== E1-E2 ========== Bridge R1-R2

beta = 50/100; %plucking position
alfa = 1; % damping factor

N_right = round(Ndel*beta);
N_left = round(Ndel*(1-beta));

[H_e1r1, ~] = freqz([zeros(1, N_right-1) 1], alfa, Nfft, Fs);
[H_e2e1, ~] = freqz([zeros(1, 2*N_left-1) 1], -alfa, Nfft, Fs);
[H_e2r1, ~] = freqz([zeros(1, (2*N_left)+N_right-1) 1], -alfa, Nfft, Fs);
[H_r2e2, F] = freqz([zeros(1, N_right-1) 1], alfa, Nfft, Fs);

H_loop = -1 .* H_r2e2 .* H_e2e1 .* H_e1r1;

S = 1./(1-H_loop);

figure(3)
plot(F, db(abs(S)));
xlim([0 500])

[integrator, ~] = freqz(1,[1 -0.999], Nfft, Fs);

figure(4)
plot(F, db(abs(integrator)));
xlim([0 500])
%%
Heb = (1+H_e2r1).*H_e1r1.*S.*Z'.*integrator;
Heb_n = Heb./max(abs(Heb));

figure(5)
plot(F, abs(Heb_n));
for i=1:4
 xline(LA1*i);
end
legend('Normalized', 'Not normalized', 'F0')
xlim([50 500]);

Lmax = 3e-3;
x2 = Lmax/(beta*L0) + Lmax/(L0*(1-beta));
x1=[x2 zeros(1,length(F)-length(x2))];

%%
y2=(fft(x1).*Heb_n');
cut=20;
y2_n=[zeros(1,cut) y2(cut+1:end)];

figure(6)
plot(F, (abs(y2_n)));
xlim([0 500]);
%%
sigy2=ifft(y2_n);
plot(1./F, sigy2)

Decay=0:1/(length(y2_n)):1-1/(length(y2_n));
PassaBass=(abs(y2_n)./max(abs(y2_n).*(1-Decay)));
signal=ifft(PassaBass).*(1-Decay);
plot(signal);
n=1;

figure(1000)
plot(n*omega/2/pi,(abs(y2_n))); 
xlim([18*n 1000])
for i=1:4
 xline(LA1*i);
end
%legend('Normalized', 'Not normalized', 'F0')

soundsc(abs(signal), Fs);

figure(1001)
t=2000*pi./omega;
plot(t,(real(ifft(PassaBass))))
xlim([0 2]);











