%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2020                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters

Fs = 4*44100;% sampling frequency [Hz]
T = 1/Fs;% sampling period [s]

% Fundamental note

C2 = 65.4;% fundamental note freq [Hz]

% Boundary          

zeta_l = 1e20;% left hand normalized impedance [(rad/s)/kg*m^-2*s^-1]
zeta_b = 1000;% bridge normalized impedance [(rad/s)/kg*m^-2*s^-1]

% String parameters

L = 1.92;% string length [m]
% Te = 750;% string tension [N]
MS = 35e-3;% string total mass [kg]
rho = MS/L;% linear string density [kg/m]
Te = 4*L^2*rho*C2^2;
c = sqrt(Te/rho);% wave velocity in the string [m/s]
epsilon = 7.5e-6;% string sitffness parameter
b_1 = 0.5;% air damping coefficient [s^-1]
b_2 = 6.25e-9;% string internal friction coefficient [s]
k = epsilon;

% Spatial sampling parameters

gamma = Fs/2/C2;

% Aliasing condition

F_Nyq = Fs/2;% nyquist frequency [Hz]

if F_Nyq < C2
    disp('ATTENTION! ALIASING CONDITION NOT MET!');
end

% Number of maximum spatial steps

N_max = L/(sqrt(0.5*(c^2*T^2+4*b_2*T+sqrt((c^2*T^2+4*b_2*T)^2+16*k^2*T^2))));
N_min = ((-1+(1+16*epsilon*gamma^2)^(1/2))/(8*epsilon))^(1/2);

% Integer values

% N = 521 %Saitis
N = floor(N_min);

% Spatial sampling

X = L/N;% spatial sample length [m]

% FD parameters

mu = k^2/(c^2*X^2);% 
nu = 2*b_2*T/X^2;% 
lambda = c*T/X;% Courant number

% Hammer parameters

MH = 4.9e-3;% mass of the hammer [kg]
b_H = 1e-4;% fluid damping coefficient [s^-1]
a = 0.12;% relative striking position
w = 0.2;% window length [m]
Vh0 = 2.5;% initial hammer velocity [m/s]

% Hammer contact window definition

M0 = floor((a*L)/X);% sample of contact (center of window)
wlength = floor(w/X);
center = round(wlength/2);
begin = M0-center+1;
g = zeros(N, 1);
g(begin:begin+wlength-1) = hann(wlength);

% PDE Coefficients:

a1 = (-lambda^2*mu)/(1 + b_1*T);
a2 = (lambda^2 + 4*lambda^2*mu + nu)/(1 + b_1*T);
a3 = (2 - 2*lambda^2 - 6*lambda^2*mu - 2*nu)/(1 + b_1*T);
a4 = (-1 + b_1*T + 2*nu)/(1 + b_1*T);
a5 = -nu/(1 + b_1*T);
af = (T^2/rho)/(1 + b_1*T); %force coefficient

% Bridge boundary coefficients

br_1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+zeta_b*lambda);
br_2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+zeta_b*lambda);
br_3 = (-2*lambda^2*mu)/(1+b_1*T+zeta_b*lambda);
br_4 = (-1-b_1*T+zeta_b*lambda)/(1+b_1*T+zeta_b*lambda);
br_5 = (T^2/rho)/(1+b_1*T+zeta_b*lambda);

% Left hand (hinged string end) boundary coefficients

bl_1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+zeta_l*lambda);
bl_2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+zeta_l*lambda);
bl_3 = (-2*lambda^2*mu)/(1+b_1*T+zeta_l*lambda);
bl_4 = (-1-b_1*T+zeta_l*lambda)/(1+b_1*T+zeta_l*lambda);
bl_5 = (T^2/rho)/(1+b_1*T+zeta_l*lambda);

% Hammer felt parameters

K = 4e8;% Hammer felt stiffness
p = 2.3;% hammer fel stiffness exponent

%% Computation of the FD scheme

% Initialization

dur = 8;% duration [s]
t = linspace(0, dur, Fs*dur);

y = zeros(N, length(t));
eta = zeros(1, length(t));

d1 = 2/(1+(b_H*T)/2/MH);
d2 = (-1+(b_H*T)/2/MH)/(1+(b_H*T)/2/MH);
df = (-T^2/MH)/(1+(b_H*T)/2/MH);

% Computation loop

%First three steps
%Step 1: already ok
%Step 2
eta(2) = Vh0*T;% see Chaigne et Al.
Fh = K*abs(eta(2)-y(M0,2))^p;

%step 3
%calculation of displacement in step 3
%first calculate approximate string diplacement
for m = 2:N-1
y(m, 3) = y(m-1,2)+y(m+1,2)-y(m,1)+(T^2*N*Fh*g(m))/MS;% see Chaigne pg.5
end

%then hammer displacement
eta(3) = 2*eta(2)-eta(1)-(T^2*Fh)/MH;% see Chaigne pg.5

%begin the main loop
for n = 3:length(t)-1
    if eta(n)<y(M0,n)
        Fh = 0;
    else
        Fh = K*abs(eta(n) - y(M0,n))^p;
    end
    F = Fh*g;
    
    for m = 1:N
        if (m-1) == 0
            y(m,n+1) = bl_1*y(m,n)+bl_2*y(m+1,n)+bl_3*y(m+2,n)+bl_4*y(m,n-1)+bl_5*F(m);
        elseif (m-1) == 1
            y(m,n+1) = a1*(y(m+2,n)-y(m,n)+2*y(m-1,n))+a2*(y(m+1,n)+y(m-1,n))...% a1*(y(m+2,n)-y(m,n)+2*y(m-2,n)) non è possibile m-2...
                +a3*y(m,n)+a4*y(m,n-1)+a5*(y(m+1,n-1)+y(m-1,n-1))+af*F(m);
        elseif (m-1) == N-2
            y(m,n+1) = a1*(2*y(m+1,n)-y(m,n)+2*y(m-2,n))+a2*(y(m+1,n)+y(m-1,n))...
                +a3*y(m,n)+a4*y(m,n-1)+a5*(y(m+1,n-1)+y(m-1,n-1))+af*F(m);
        elseif (m-1) == N-1
            y(m,n+1) = br_1*y(m,n)+br_2*y(m-1,n)+br_3*y(m-2,n)+br_4*y(m,n-1)+br_5*F(m);
        else
            y(m,n+1) = a1*(y(m+2,n)+y(m-2,n))+a2*(y(m+1,n)+y(m-1,n))...
                +a3*y(m,n)+a4*y(m,n-1)+a5*(y(m+1,n-1)+y(m-1,n-1))+af*F(m);
        end
    end
       eta(n+1) = d1*eta(n)+d2*eta(n-1)+df*Fh;
end

%% Plot the displacement in time
space = linspace(0,L,N);

figure(1)
for i =  1:50:2000
    sgtitle(['Time: ' num2str(t(i)) ' s'])
    plot(space,y(:,i));
    ylim([-2.5e-4 2.5e-4]);
    xlim([0 L]);
    xlabel('String length [m]')
    ylabel('Displacement [m]')
    pause(T);
end

%% Plot the synthesized signal play it and save it on the disk

spec_pos = N-M0;

y_avg = sum(y(spec_pos-5:spec_pos+6,:), 1)/12;

figure(2)
plot(t, y_avg);
xlabel('Time [s]');
ylabel('Avg displacement [m]');

% Play the sound

soundsc(y_avg, Fs);

% Save on disk

sound = decimate(y_avg,4,'fir');

sound = (sound/max(abs(sound)))*0.9;

audiowrite('C2note.wav', sound, Fs/4);


