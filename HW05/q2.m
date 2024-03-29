clc, clear all, close all;
%% Question2
L = 0.645;   % [m] Length of the string
fE2 = 82;   % [Hz]
mu = 5.78e-3;    % [kg m^-1] linear density
T_str = 4*(L^2)*mu*fE2^2;   % [N] tension
c = sqrt(T_str/mu);
beta = 1/5;

Fs = 44100; %sampling frequency
T_E2E1 = (2*L*(1-beta))/c;
T_E1R1 = (beta*L)/c;
T_S = (2*L)/c;

N_E = round(T_E2E1*Fs);
N_E1R1 = round(T_E1R1*Fs);
N_S = round(T_S*Fs);

%% Digital filter design with Simulink

%filter from Simulink
signalLen = 5; %duration of simulation

%load simulation
open_system(['ques2.slx'], 'loadonly');

%perform simulation
simulation = sim(['ques2.slx'], signalLen);

input = simulation.input.data;
output = simulation.filter.data;

%bridge impedance from q1
% Parameters of the bridge 
k_p = 1.41e5; % [N/m] top plate stiffness
m_p = 0.128*0.385;    % [kg] top plate effective mass
A_p = 0.0375*0.385;   % [m^2] top plate area
R_p = 32; % [N s m^-1] top plate resistance

m_h = 8.04e-4;   % [kg]  air piston mass
A_h = 7.85e-3;   % [m^2] air piston area
R_h = 30;    % [Ns/m^5] sound hole resistance

V = 0.0172;  %[m^3] air cavity volume
rho = 1.2;   % [kg m^-3] air density
c = 343;    % [m/s] sound velocity in air

L_p = m_p/(A_p^2);  % plate inertance
L_h = m_h/(A_h^2);  % sound hole inertance
C_p = A_p^2/k_p;    % plate compliance
C_v = V/(rho*c^2);  % cavity compliance

%plot of the bridge impedance

f = 1:Fs/length(input):Fs-(1/length(input)) +1;  %[Hz] range of frequency in Hz
w = 2*pi*f; %[rad/s] range of angular frequency 

Z_p = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1); %[Kg/m^4s] plate impedance
Z_v = 1./(1i.*w*C_v); %[Ks/m^4 s] cavity impedance
Z_h = 1i.*w*L_h + R_h; %[Ks/m^4 s] hole impedance

Z_b = Z_p + ( Z_v .* Z_h ) ./ (Z_v + Z_h ); %[Ks/m^4 s] bridge impedance

%filter without bridge impedance
H_1 = (abs(fft(output) ./ fft(input)))';

%filter with bridge impedance
Z_b(1) = 1;
maxZ = maxk(abs(Z_b), 2);
H = H_1.*(Z_b./(maxZ(2)));


figure(1)
semilogy(f, abs(H)./(max(abs(H))));
xlim([0 500])
ylim([0 1])
grid on

%% analytical q2

%fare una formula analitica di H, parametrica rispetto a beta (definire i
%delay time rispetto a beta e rispetto alla nota (il valore di c cambia a
%seconda della frequenza (possibile? non renderebbe il filtro non lineare?)
%provare a plottare diversi filtri, magari fissando c e cambiando beta
%(punto di applicazione) o fissando beta e cambiando c (nota suonata)
%ricordati di plottare anche la fase dei filtri

% analytical filter in Laplace domain

H_E = 0.5*(1 + 0.99*exp(-T_E2E1*1i.*w));

H_E1R1 = 0.99*exp(-T_E1R1*1i.*w);

%H_ma = 0.5.*((1 - exp(-2*(1/Fs).*1i.*w))./(1 - exp(-(1/Fs)*1i.*w))); 

S = 1./(1 - 0.99*exp(-T_S.*1i.*w));

H_EB = 2.*H_E.*(Z_b./(maxZ(2))).*H_E1R1.*S.*(1./(1i.*w));

H_EB(1) = 0;
figure(2)
semilogy(f, abs(H_EB)./(abs(max(H_EB))))
ylim([0 1])
xlim([0 500])
grid on

%%  q3

t = linspace(0, signalLen,  Fs*signalLen + 1); %time vector

ai = zeros(1, length(t)); %impulsive signal = displacement at 0.003 [m]
ai(1) = 10;

Ai = fft(ai, Fs*signalLen + 1); %fft of input signal

Fo = Ai.*H_EB; %system frequency response 

Fo(1) = 1; %remove singolarities

fo = ifft(Fo, Fs*signalLen + 1, 'nonsymmetric'); %force at the bridge = system time response

figure(3)
plot(t, abs(fo));



% filename = 'force.wav';
% audiowrite(filename, 100000.*abs(fo),Fs);
