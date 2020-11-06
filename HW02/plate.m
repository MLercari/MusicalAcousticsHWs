%% HOMEWORK 2
close all, clear all, clc;

%% PLATE PARAMETERS
h = 3e-3;   % thickness [m]
E = 69e9;   % Young's modulus [N/m^2]
nu = 0.35;  % Poisson ratio
rho = 2650; % density [kg/m^3]
Lx = 1; % [m]
Ly = 1.4;   % [m]

a = 4e-4; % amplitude [m] initial amplitude 
alpha = 1e3; % decay factor [s^-1]

c_L = sqrt( E/(rho*(1-(nu^2))) ); %longitudinal wave speed
D = E*h^3/(12*(1 - nu^2)); %bending stiffness
sigma = rho * h; %densit� superficiale
%% 
% HW2 -1 

% SS simply supported analytical computation
f = zeros(5,5);


for i=0:5
    for j=0:5
        f(i+1, j+1) = 0.453*c_L*h*( ((i+1).^2/Lx^2) + ((j+1).^2/Ly^2) );
    end
end

f_a_ss = [f(1,1) , f(1,2), f(2,1) , f(1,3), f(2,2)];

% analytical freq.
% preso dal Fletcher e dalle slides

% clamped 
f00 = 1.654*c_L*h/(Lx^2);
f_a_c = zeros(1,5);
relative = [0.75, 1.88, 1.16, 2.27, 3.66]; %il quinto modo � in realt� quadrato
for i = 1:5
    f_a_c(i) = f00*relative(i);
end

% free
f11 = h*c_L*sqrt((1-nu)/2)/(1.2^2);
f_a_f = zeros(1,5);
relative_f = [1, 1.52, 1.94, 2.71, 2.71];
for i = 1:5
    f_a_f(i) = f11*relative_f(i);
end

% dalla Nasa tab. 4.29 assumento aspect ratio a/b = 0.72
% clamped 

% a = Lx

%nasa_coeff = [27.01, 41.72, 65.5, 66.53, 79.81];
nasa_coeff = [28.05, 45.72, 67.005, 75.06, 83.34];
%             (0,0)  (0,1)  (1,0) (0,2)  (1,1)  (m,n) 

% w a^2 sqrt(rho / D) = coeff 
% w = coeff/(a^2*sqrt(rho/D)
w_nasa_c = zeros(1,5);
for i = 1:5
    w_nasa_c(i) = nasa_coeff(i)/(Lx^2*sqrt(sigma/D));
end

f_nasa_c = w_nasa_c/(2*pi);

% free tabella 4.72
% assumendo nu = 1/6 e b/a = 1.5 
nasa_coeff_f = [9.905, 9.944, 22.245, 22.373, 27.410]; 
%             (1,1)  (0,2)    (1,2)      (2,0)       (3,0)  (m,n) 

% w a^2 sqrt(rho / D) = coeff 
% w = coeff/(a^2*sqrt(rho/D)
w_nasa_f = zeros(1,5);
for i = 1:5
    w_nasa_f(i) = nasa_coeff_f(i)/(Lx^2*sqrt(sigma/D));
end

f_nasa_f = w_nasa_f/(2*pi);

%% HW-2


% difference between adjacent modes

delta_ss = zeros(1,4);

delta_f = zeros(1,4);

delta_c = zeros(1,4);


for i = 1:4
    delta_f(i) = f_nasa_f(i+1) - f_nasa_f(i);
end


for i = 1:4
    delta_ss(i) = f_a_ss(i+1) - f_a_ss(i);
end


for i = 1:4
    delta_c(i) = f_nasa_c(i+1) - f_nasa_c(i);
end

% span between last and first modes
span_f = f_nasa_f(5) - f_nasa_f(1);
span_ss = f_a_ss(5) - f_a_ss(1);
span_c = f_nasa_c(5) - f_nasa_c(1);

plot(1:4, delta_f, 1:4, delta_ss, 1:4, delta_c);
grid on
xlabel('mode')
ylabel('\Delta f')
legend('free', 'simply supported', 'clamped');

%% HW2-3

% Time history of the frequencies 

t = linspace(0, 0.01, 100);

damp = a * exp(-alpha*t); % time decay of the amplitude

f_non_lin = [f_nasa_f, f_a_ss, f_nasa_c];
f_non_lin = (1 + 0.16.*(damp./h).^2).*f_non_lin';

% %plot
% subplot(5,3,1)
% plot(t, f_non_lin(1,:)) ; %1 mode free 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('first mode free')
% 
% subplot(5,3,4)
% plot(t, f_non_lin(2,:)) ; %2 mode free
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('second mode free')
% 
% subplot(5,3,7)
% plot(t, f_non_lin(3,:)) ; %3 mode free 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('third mode free')
% 
% subplot(5,3,10)
% plot(t, f_non_lin(4,:)) ; %4 mode free 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fourth mode free')
% 
% subplot(5,3,13)
% plot(t, f_non_lin(5,:)) ; %5 mode free 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fifth mode free')
% 
% subplot(5,3,2)
% plot(t, f_non_lin(6,:)) ; %1 mode ss 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('first mode supported')
% 
% subplot(5,3,5)
% plot(t, f_non_lin(7,:)) ; %2 mode supported
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('second mode supported')
% 
% subplot(5,3,8)
% plot(t, f_non_lin(8,:)) ; %3 mode supported 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('third mode supported')
% 
% subplot(5,3,11)
% plot(t, f_non_lin(9,:)) ; %4 mode supported 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fourth mode supported')
% 
% subplot(5,3,14)
% plot(t, f_non_lin(10,:)) ; %5 mode supported 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fifth mode supported')
% 
% subplot(5,3,3)
% plot(t, f_non_lin(11,:)) ; %1 mode clamped  
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('first mode clamped ')
% 
% subplot(5,3,6)
% plot(t, f_non_lin(12,:)) ; %2 mode clamped 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('second mode clamped ')
% 
% subplot(5,3,9)
% plot(t, f_non_lin(13,:)) ; %3 mode clamped  
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('third mode clamped ')
% 
% subplot(5,3,12)
% plot(t, f_non_lin(14,:)) ; %4 mode clamped  
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fourth mode clamped ')
% 
% subplot(5,3,15)
% plot(t, f_non_lin(15,:)) ; %5 mode clamped 
% grid on;
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% legend('fifth mode clamped ')

%% HW2-4
% space definition

xSpatialNPoints = 100;
ySpatialNPoints = 100;

x = linspace(0 , Lx, xSpatialNPoints);
y = linspace(0 , Ly, ySpatialNPoints);
z_1 = zeros(length(y), length(x));

% first approach: approximation 2 closest eigen modes
m1 = 4; %relativo al lato corto
n1 = 2; %relativo al lato lungo
m2 = 3;
n2 = 4;

for i = 1:length(x)
    for j = 1:length(y)
    z_1(j,i) = sin(m1*pi*x(i)/ Lx)*sin(n1*pi*y(j)/Ly)+ sin(m2*pi*x(i)/ Lx)*sin(n2*pi*y(j)/Ly);
    end
end
Z_1 = z_1.*1i*130.8*2*pi;
%surf(x, y, abs(Z_1));

%search maxima value around a given interval (there should be 4 points)
v_a = abs(Z_1);
M = max(v_a, [], 'all');
[y1,x1]=find(v_a >= M -1e-3);
x1 = x1*(Lx/xSpatialNPoints);
y1 = y1*(Ly/ySpatialNPoints);



%% Hw2-5
%% HW2-6
