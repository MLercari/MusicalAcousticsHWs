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
%LA = 80;% [Hz]
LA1=440;
Ndel = round(Fs/LA1); %Number of delay samples for delay line
omega = linspace(0, pi*(Fs),Nfft);

Zb1=1i*omega.*M_h+R_h;
Zb2=(1./(1i*(C_v*omega)));
Zeq=((1./Zb1 + 1./Zb2).^-1);
Zb=Zeq+r_p+1i*(M_p*omega-(1./(C_p*omega)));

%H_ = dfilt.delay(latency)

beta = 50/100; %plucking position
%alfa = 1; % damping factor

N_right = round(Ndel*beta);
N_left = round(Ndel*(1-beta));

H_e1r1 = dfilt.delay(N_right);
H_e2e1 = dfilt.delay(2*N_left);
H_e2r1 = dfilt.delay(2*N_left+N_right);
H_r2e2 = dfilt.delay(N_right);
x=rand(1,Ndel);
x=[x zeros(1,Fs-Ndel)];
y=filter(H_e2r1,x);
H_e2r1FFT=fft(y)./fft(x);
y=filter(-1,1,x);
y=filter(H_r2e2,y);
y=filter(H_e2e1,y);
y=filter(H_e1r1,y);
H_loop = fft(y)./fft(x);
%H_loop = -1 .* H_r2e2 .* H_e2e1 .* H_e1r1;
S = 1./(1-H_loop);
Snorm = S./max(abs(S));
plot(1:length(x),(abs(S)));
xlim([1 500])

[integrator, ~] = freqz(1,[1 -1], Nfft, Fs);

%Heb = (H_e2r1).*H_e1r1.*S.*Z'.*integrator+H_e1r1.*S.*Z'.*integrator;
%% test segnale
%Devo tagliare Zb per evitare che mandi tutto a caso
ZbCut=[zeros(1,25) Zb(26:end)];
integrator1=integrator';
integratorCut=[zeros(1,25) integrator1(26:end)];
integratorNorm=integratorCut./(max(abs(integratorCut)));
x2=ones(1,1); %impulso in ingresso (si può cambiare)
x1=[x2 zeros(1,Fs-length(x2))];
y1=filter(H_e1r1,x1);
y2=(fft(y1).*Snorm);
numb=40;
Znorm=[zeros(1,numb) Zb(numb+1:end)]./(max(Zb(numb:22000)));
y3=(y2.*Znorm);
y4=y3.*(integratorNorm); %H_e1r1.*S.*Z'.*integrator;
y5=y4.*H_e2r1FFT;
y6=(y4+y5);
y6Norm=[zeros(1,round(LA1/2)) y6(LA1/2:Fs/2)./max(abs(y6(LA1/2:Fs/2)))];
Decay=0:1/(length(y6Norm)):1-1/(length(y6Norm));
%il decay dovrebbe essere calcolato così
%=================================
%ro=1.2; %kg/m3 densità aria
%nu=1.81e-5; %kg/ms viscosità dinamica aria a 15°
%diam = 4e-4; %m diametro corda
%fc=1:500;
%R_L=2*pi*nu+2*pi*diam*sqrt(pi*nu*ro.*fc);
%Q_friction=2*pi*mu.*fc./R_L;
%Q_VeTE=0;
%Q_0=Q_friction+Q_VeTE;
%Decay = -pi.*fc./Q_0;
%plot(fc,Decay);
%============================
signal=ifft(abs(y6Norm).*(1-Decay)).*(1-Decay);
%plot(signal);
n=1;

figure(100)
plot(n*omega(1:Fs/2+1)/2/pi,(abs(y6Norm))); 
xlim([18*n 10000])
for i=1:4
 xline(LA1*i);
end
%legend('Normalized', 'Not normalized', 'F0')
soundsc(abs(signal))


