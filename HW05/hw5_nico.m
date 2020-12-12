clear;
close all;
clc;
%%
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

M_p = m_p/A_p^2;
M_h = m_h/A_h^2;
C_p = A_p^2/k_p;
C_v = V/rho/c^2;
R_p = 32; %Nm/kg/s
R_v = 0; %kg/s

%% 
dur = 4;
beta = 1/2;
L0 = 0.65;
%A = 0.36e-6;
E = 5e9;
T0 = 82;
d = 3e-3;
mu = 7.2e-3;
Fs = 100000; %usare questa per gli spettri 
%Fs = 88100; %usare questa per la risposta nel tempo

E4 = 329.63;
B = 246.94;
G = 196.00;
D = 146.83;
A = 110;
E1 = 82.41;

notes = [E1, A, D, G, B, E4];
H = [];
outputs = [];
t = 0:1/Fs:dur-1/Fs;
f = 0:Fs/length(t):Fs-1/length(t);
w = 2*pi*f;

Z_eq = zeros(1, length(t));
for jj = 1:length(t)
    Z_eq(jj) = (1i*w(jj)*M_p + 1/(1i*w(jj)*C_p) + R_p + ...
        (1/(1i*w(jj)*C_v) + R_v )*((1i*w(jj)*M_h + R_h)/(1i*w(jj)*M_h + R_h ...
        + 1/(1i*w(jj)*C_v) + R_v)));%*A_p;
end
Y_eq = 1./Z_eq;
Z = Y_eq/max(abs(Y_eq));
Z(1)=0;
Z = Z_eq/max(abs(Z_eq(2:end)));

% figure(1000)
% plot(w/2/pi, abs(Z));
% xlim([0 500]);


for ii = 1:length(notes)
    N = round(Fs/notes(ii));
    N_right = round(N*beta);
    N_left = round(N*(1-beta));
    
    x = linspace(0, L0, N);
    
    H_e1r1 = freqz([zeros(N_right,1);1],1,length(t));
    H_e2r1 = freqz([zeros(N_left + N, 1);-1],1, length(t));
    %H_e2e1 = freqz([zeros(2*N_left, 1);-1],1, length(t));
    %H_r2e2 = H_e1r1;
    H_loop = freqz([zeros(2*N, 1);1],1, length(t));
    H_loop(1) = 1.01;
    %H_loop = R_b*H_r2e2.*H_e2e1.*H_e1r1;

    R_b = -1;

    integral = freqz([1], [1,-1], length(t));

   %H_eb = 1/2*((1 + H_e2r1).*H_e1r1.*Z'.*integral.*(1-R_b))./(1-H_loop);
    H_eb = 1/2*((1 + H_e2r1).*H_e1r1.*Z'.*integral)./(1-H_loop);
    H_eb = abs(H_eb).*exp(1i*angle(H_eb)) - abs(H_eb).*exp(-1i*angle(H_eb));
    H(ii,:) = H_eb;
    
    Lmax = 3e-3;
    deltaamp = Lmax/(beta*L0) + Lmax/(L0*(1-beta));
    acc = zeros(size(t));
    acc(1) = -deltaamp;
    spec_acc = fft(acc);
    spec_out = spec_acc.*H_eb';
    output = ifft(spec_out).*exp(-2*t);
    outputs(ii,:) = output;

    figure(5)
    subplot (6,1,ii)
    semilogy(f, abs(H_eb), 'Linewidth', 1.1);
    xlim([0 500]);
    ylabel('|H_{eb}|'); xlabel('Frequency [Hz]');
    figure(6)
    subplot (6,1,ii)
    plot(f, angle(H_eb)*180/pi, 'LineWidth', 1.1);
    xlim([0 500]);
    ylabel('angle H_{eb}'); xlabel('Frequency [Hz]');
    figure(7)
    subplot(6,1,ii)
    semilogy(f, abs(spec_out));
    xlim([0 500]);
    figure(8)
    subplot(6,1,ii)
    plot(t, abs(output));

end

%% Es3 Plucking position va a 1/5!!!
Lmax = 3e-3;
deltaamp = Lmax/(beta*L0) + Lmax/(L0*(1-beta));

acc = zeros(size(t));
acc(1) = -deltaamp;
spec_acc = fft(acc);
spec_out = spec_acc.*H_eb';
output = ifft(spec_out).*exp(-2*t);

figure(100)
plot(t, abs(output));
