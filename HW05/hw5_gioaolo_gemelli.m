clear
close all
clc

%%
dur = 4;
beta = 1/3; %position coefficient
L0=0.64516; %lunghezza corda
%A = 0.36e-6; %boh
E = 5e9; %modulo young corda nylon
Eacciaio = 1.875e11; %Pa modulo young corda acciaio https://www.alloywire.it/products/stainless-steel-302/
d = 3e-3; %initial displacement of string
mu = 7.2e-3; %nylon densità lineare
muacciaio = sqrt(Eacciaio/E)*mu; %densità lineare acciaio
Fs = 44100; %frequenza di campionamento
FC = 100; %Hz frequenza corda
D = 1/FC/2; %Delay time
%?????????????????????????????????????????
N = round(Fs/FC); %Delay in sample totale
N_right = round(N*beta); %delay in sample pluck ponte
N_left = round(N*(1-beta)); %delay in sample pluck nut

%%TOP PLATE
kp=1.41*100000; %stiffness top plate N/m
mplate=0.128; %real mass in kg
corrMFactor=0.385; %fattore di riduzione
mp=mplate*corrMFactor; %modal mass
Ap=0.0375*0.385; %top plate area m2
rp=32; %top plate resistance Ns/kg
% SOUNDHOLE
mh=0.000804; %kg of the air piston
Ah=0.00785; %area of the air piston
R_h = 30; %air resistance N/m
%AIR CAVITY
V=0.0172; %volume m3
R_v=0;
maxf=500;
f=50:0.1:maxf;
omega=2*pi.*f;

%% QUESTION 1
M_p=mp/(Ap^2);
M_h=mh/(Ah^2);
C_p=Ap^2/kp;
C_v=V/1.2/(340^2);
R_p=rp;
%%
t = 0:1/Fs:dur-1/Fs;
x = linspace(0, L0, N);
f = 0:Fs/length(t):Fs-1/length(t);
w = 2*pi*f;
%%
H_e1r1 = freqz([zeros(N_right,1);1],1,length(t));
H_e2r1 = freqz([zeros(N_left+N,1);-1],1, length(t));
%H_e2e1 = freqz([zeros(2*N_left, 1);-1],1, length(t));
%H_r2e2 = H_e1r1;
H_loop = freqz([zeros(2*N, 1);1],1,length(t));

%%

R_b = -1;

integral = freqz([1], [0,1], length(t));

Z_eq = zeros(1, length(t));
for jj = 1:length(t)
    Z_eq(jj) = (1i*w(jj)*M_p + 1/(1i*w(jj)*C_p) + R_p + ...
        (1/(1i*w(jj)*C_v) + R_v )*((1i*w(jj)*M_h + R_h)/(1i*w(jj)*M_h + R_h ...
        + 1/(1i*w(jj)*C_v) + R_v)));%*A_p;
end

Y_eq = 1./Z_eq;
Z = Y_eq/max(abs(Y_eq));
Z = Z_eq/max(abs(Z_eq(2:end)));


%H_eb = 1/2*((1 + H_e2r1).*H_e1r1.*Z'.*integral.*(1-R_b))./(1-H_loop);
H_eb = 0.5*((1 + H_e2r1).*H_e1r1.*Z'.*integral)./(1-H_loop);
H_eb = abs(H_eb).*exp(1i*angle(H_eb)) - abs(H_eb).*exp(1i*(-angle(H_eb)));


figure(5)
subplot 211
semilogy(f, abs(H_eb), 'Linewidth', 1.1);
xlim([0 500]);
ylabel('|H_{eb}|'); xlabel('Frequency [Hz]');
subplot 212
plot(f, angle(H_eb)*180/pi, 'LineWidth', 1.1);
xlim([0 500]);
ylabel('angle H_{eb}'); xlabel('Frequency [Hz]');