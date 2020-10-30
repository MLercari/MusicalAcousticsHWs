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

