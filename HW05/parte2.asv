clc, clear all, close all;

%% Parameters
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

%% Question2

L = 0.64516;   % [m] Length of the string
fE2 = 82;   % [Hz]
mu = 5.78e-3;    % [kg m^-1] linear density
T_str = 4*(L^2)*mu*fE2^2;   % [N] tension
c = sqrt(T_str/mu);
beta = 1/5;


Fs = 20*c/(beta*L); %sampling frequency
T_E2E1 = (2*L*(1-beta))/c;
T_E1R1 = (beta*L)/c;
T_S = (2*L)/c;

duration = 2.5;   % [s]
t = 0:1/Fs:duration;

N_E = round(T_E2E1*Fs);
N_E1R1 = round(T_E1R1*Fs);
N_S = round(T_S*Fs);
%{
syms z;
H_E = 0.5*(1 - 0.99*(z^(-N_E)));
H_E1R1 = 0.99*(z^(-N_E1R1));
H_loop = 0.5*(1 + z^(-1))*z^(-N_S);
H_S = 1/(1 - H_loop);
integr = 2/(1-0.9*z^(-1));
%}

[H_E, ~, ~] = freqz([0.5; zeros(N_E - 1, 1); 0.5*0.8],1, length(t), Fs);
[H_E1R1, ~, ~] = freqz([zeros(N_E1R1, 1); 0.8],1, length(t), Fs);
%[H_S, ~, ~] = freqz(1, [zeros(N_S, 1); -0.5; -0.5], length(t), Fs);
movAvg = freqz([0.5 0.5], 1, length(t), Fs);
loopDel = freqz([zeros(N_S,1); 1], 1, length(t), Fs);
H_S = 1./(1 - 0.999.*movAvg.*loopDel);
[integr, f] = freqz(2, [1; 0.9], length(t), Fs);

%f = w/(2*pi);


%f = 0:1/Fs:; %[Hz] range of frequency in Hz
w = 2*pi*f; %[rad/s] range of angular frequency 

Z_p = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1); %[Kg/m^4s] plate impedance
Z_v = 1./(1i.*w*C_v); %[Ks/m^4 s] cavity impedance
Z_h = 1i.*w*L_h + R_h; %[Ks/m^4 s] hole impedance

Z = Z_p + ( Z_v .* Z_h ) ./ (Z_v + Z_h ); %[Ks/m^4 s] bridge impedance
Z = 2*Z./(1i.*w);
Z = Z.*(A_p^2);

Z(1) = 0;
%Z = Z./(max(abs(Z(50:length(Z))))); %perchè z a 0 + infinito


H1 = H_E.*H_E1R1.*H_S.*Z;


%% quest 3
disp = 0.003;
Teq = T_str*(disp/(beta*L) + (disp/(1 - beta)*L));

a0 = zeros(length(t), 1);
%a0(1) = (Teq*20)/(beta*L*mu);
%a0(1) = disp*Fs/c;  %sarà giusto?
a0(1) = (c^2).*(disp/(beta*L) + (disp/(1 - beta)*L));

imp = fft(a0, length(t));

out = ifft(imp.*H1, length(t));


figure(69);
subplot 211;
semilogy(f, abs(H1));
xlim([0 500]);
xlabel('f [Hz]');
ylabel('|H_{EB}| [kg]');
grid on;

subplot 212;
plot(f, angle(H1));
xlim([0 500]);
xlabel('f [Hz]');
ylabel('\angle H_{EB}');
grid on;

%{
figure(420);
plot(t, real(out));
%}


filter = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',10,'PassbandRipple',0.2, ...
         'SampleRate',Fs);
     
[filt, ~] = freqz(filter, length(H1));
out2 = ifft(imp.*H1.*filt, length(t));

figure(666)
plot(t, real(out2));
hold on;
plot(t, zeros(1,length(t)), 'LineStyle', '--', 'Color', 'k');
grid on;
xlabel("t [s]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex');
ylabel("F [N]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')


figure(777)
semilogy(f, abs(fft(out2, length(f))));
ylim([ 0 1e4]);
xlim([0 600]);

%% displacement function
x = 0:c/Fs:L;
%disp
y_x = (disp/(beta*L)).*x.*(heaviside(x) - heaviside(x - beta*L)) +( disp - (disp/(L - beta*L)).*(x - beta*L)).*heaviside(x - beta*L);
plot(x , y_x);
xlim([ 0 L])
ylim([0 disp])
%vel
v_x = diff(y_x)./diff(x);
x(end) = [];
plot(x , v_x);
%curv
a_x = diff(v_x)./diff(x);
x(end) = [];
plot(x , a_x);
%% Simulink (non funzionante)

%{

%load simulation
open_system(['ques2.slx'], 'loadonly');

%perform simulation
simulation = sim(['ques2.slx'], duration);

input = simulation.input.data;
output = simulation.filter.data;

%filter without bridge impedance
Hs = (abs(fft(output) ./ fft(input)))';

Hsz = Hs'.*Z;

figure(22)

subplot 211;
semilogy(f, abs(Hsz));
xlim([0 500]);
xlabel('f [Hz]');
ylabel('|H_{EB}| [kg]');
grid on;

subplot 212;
plot(f, angle(Hsz));
xlim([0 500]);
xlabel('f [Hz]');
ylabel('\angle H_{EB}');
grid on;

out3 = ifft(imp.*Hsz.*filt, length(t));

figure(667)
plot(t, real(out3));
hold on;
plot(t, zeros(1,length(t)), 'LineStyle', '--', 'Color', 'k');
grid on;
xlabel("t [s]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex');
ylabel("F [N]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')

%}
