%% HOMEWORK 2
close all, clear all, clc;

%% PLATE PARAMETERS
h = 3e-3;   % thickness [m]
E = 69e9;   % Young's modulus [N/m^2]
nu = 0.35;  % Poisson ratio
rho = 2650; % density [kg/m^3]

a = 4e-4; % amplitude [m]
alpha = 1e3; % decay factor [s^-1]



Lx = 1; % [m]
Ly = 1.4;   % [m]

c_L = sqrt( E/(rho*(1-(nu^2))) );

f = zeros(5,5);

for i=0:5
    for j=0:5
        f(i+1, j+1) = 0.453*c_L*h*( ((i+1).^2/Lx^2) + ((j+1).^2/Ly^2) );
    end
end


%% Time history of the frequencies (Q3)

ff = [6.9614 ,8.0668, 16.457, 16.712, 20.660];

fss = [11.218, 22.627, 33.715	, 41.825, 45.368];

fc = [21.391 , 34.885, 51.555, 57.354, 64.806];

t = linspace(0, 0.02, 100);

damp = a * exp(-alpha*t); % time decay of the amplitude

w1 = 2*pi*ff(1)*(1 + 0.16.*(damp./h).^2);
semilogx(t, w1);
hold on;

w2 = 2*pi*ff(2)*(1 + 0.16.*(damp./h).^2);
semilogx(t, w2);
grid on;

%% analytical freq.
% preso dal Fletcher e dalle slides

% clamped 
f00 = 1.654*c_L*h/(Lx^2);
f_a_c = zeros(1,5);
relative = [0.75, 1.88, 1.16, 2.27, 3.66]; %il quinto modo è in realtà quadrato
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

%% dalla Nasa tab. 4.29 assumento aspect ratio a/b = 0.667

% clamped 

D = E*h^3/(12*(1 - nu^2));

% a = Lx
sigma = rho * h;
nasa_coeff = [27.01, 41.72, 65.5, 66.53, 79.81];
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

%% punto 4
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
Z_1 = z_1.*(1i*130.8*(2*pi) -1000);
%surf(x, y, abs(Z_1));


%search maxima value around a given interval (there should be 4 points)
v_a = abs(Z_1);
M = max(v_a, [], 'all');
[y1,x1]=find(v_a >= M-1e-6);

