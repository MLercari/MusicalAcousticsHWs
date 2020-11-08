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


for i=0:20
    for j=0:20
        f(i+1, j+1) = 0.453*c_L*h*( ((i+1).^2/Lx^2) + ((j+1).^2/Ly^2) );
    end
end

f_a_ss = [f(1,1) , f(1,2), f(2,1) , f(1,3), f(2,2)];

% analytical freq.
% preso dal Fletcher e dalle slides

% clamped 

f00 = 1.654*c_L*h/(Lx*Ly);
f_a_c = zeros(1,5);
relative = [1, 2.04, 2.04, 3.01, 3.66]; %il quinto modo � in realt� quadrato
           %(0,0) 
for i = 1:5
    f_a_c(i) = f00*relative(i);
end

% free
f11 = h*c_L*sqrt((1-nu)/2)/(Lx*Ly); %una media delle lunghezze come L^2
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
%nasa_coeff_f = [9.905, 9.944, 22.245, 22.373, 27.410]; 
%             (1,1)  (0,2)  (1,2)    (2,0)    (3,0)  (m,n) 

% b/a = 1.5 , nu = 0.3 dalla tabella 4.75
 
% w a^2 sqrt(rho / D) = coeff 
% w = coeff/(a^2*sqrt(rho/D)
w_nasa_f = zeros(1,5);
for i = 1:5
    w_nasa_f(i) = nasa_coeff_f(i)/(Lx^2*sqrt(sigma/(E*h^3/(12*(1 - (1/6)^2)))));
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


xlabel('mode')
bar(1:4, delta_ss);
ylim([0, 16]);
grid on;
legend('supported');


subplot 133;
xlabel('mode')
bar(1:4, delta_c);
grid on;
legend('clamped');





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
    z_1(j,i) = sin(m1*pi*x(i)/ Lx)*sin(n1*pi*y(j)/Ly) + sin(m2*pi*x(i)/ Lx)*sin(n2*pi*y(j)/Ly);
    end
end
Z_1 = z_1.*1i*130.8*2*pi;
surf(x, y, abs(Z_1));

%search maxima value around a given interval (there should be 4 points)
v_a = abs(Z_1);
M = max(v_a, [], 'all');
[y1,x1]=find(v_a >= M -1e-3);
x1_burdo = x1*(Lx/xSpatialNPoints); %punti trovati con metodo semplice
y1_burdo = y1*(Ly/ySpatialNPoints);  % punti trovati con metodo semplice


%% HW2 -4 metodo pi� cazzuto

w0 = 2*pi*130.8; %input frequency
w_mn = 2*pi.*f;
Z_mn = zeros(21,21,length(x), length(y)); %Z_mn(m,n,x,y)
A_mn = 2/(sqrt(rho*h*Lx*Ly));

for i = 1:length(x)
    for j = 1:length(y)
        for m = 0:20
            for n = 0:20
                Z_mn(m+1,n+1, i, j) = A_mn*sin((m+1)*pi*x(i)/ Lx)*sin((n+1)*pi*y(j)/Ly);
            end
        end 
    end
end

  Z_11 = squeeze(Z_mn(2,2,:,:));
  surf(x,y, Z_11);

%%
Y_in = zeros(length(x), length(y));


for m = 0:20
    for n = 0:20
        Y_in = Y_in + (squeeze(Z_mn(m+1,n+1,:,:).^2)./(w_mn(m+1,n+1)^2 - w0^2 +2*1i*alpha*w0));
    end
end

Y_in = 1i*w0.*Y_in;

surf(x,y, abs(Y_in));
xlabel("x [m]")
ylabel("y [m]")
zlabel("Y [s/Kg]")

v_a = abs(Y_in);
M = max(v_a, [], 'all');
[y1,x1]=find(v_a >= M-1e-7);

x1 = x1*(Lx/xSpatialNPoints); %punti trovati con metodo cazzuto
y1 = y1*(Ly/ySpatialNPoints);  % punti trovati con metodo cazzuto
[x1,y1]


%% Hw2-5

resFreqs = [f(3,4), f(4,2)];

L_s = 1.007; % length of the string [m]
mu_s = 12e-3; % mass per unit length [Kg/m]
f_0 = 130.8; % string frequency [Hz]
w_s = f_0*2*pi; 

T = mu_s*f_0^2*L_s^2*4; %string tension [N]
Z_c = sqrt(mu_s*T); %characteristic impedance of the string [Kg/s]

x01 = int8(x1/(Lx/xSpatialNPoints)); %conversione da coordinate a indici
y01 = int8(y1/(Ly/ySpatialNPoints));  % conversione da coordinate a indici
ki = 1i*Z_c*Y_in(x01(1),y01(1))/pi; %ki
eps = linspace(-3/w_s, 3/w_s, 10000);
a = zeros(2, length(eps));
a(1,:) = ki + eps + sqrt(eps.^2 + ki^2);
a(2,:) = ki + eps - sqrt(eps.^2 + ki^2);

figure(1)
subplot(2,1,1)
plot(w_s.*eps, w_s.*real(a(1,:)), w_s.*eps, w_s.*real(a(2,:))); 
subplot(2,1,2)
plot(w_s.*eps, w_s.*imag(a(1,:)), w_s.*eps, w_s.*imag(a(2,:))); 


%% HW2-6

beta = zeros(2,length(eps));

beta(1,:) = w_s + a(1,:).*w_s; %solution of the mechanical quantities defined as complex 
beta(2,:) = w_s + a(2,:).*w_s;
tau_1 = imag(beta(1,:)).^-1;
tau_2 = imag(beta(2,:)).^-1;

figure(2)
subplot(2,1,1)
plot(w_s.*eps, (real(beta(1,:))), w_s.*eps, (real(beta(2,:))) ); %eigenfrequencies
legend("a+","a-");
subplot(2,1,2)
plot(w_s.*eps, tau_1 , w_s.*eps, tau_2 ); %decay time
legend("a+","a-");
ylim([-2 2])
