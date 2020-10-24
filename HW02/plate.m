%% HOMEWORK 2
%% PLATE PARAMETERS
h = 3e-3;   % thickness [m]
E = 69e9;   % Young's modulus [N/m^2]
nu = 0.35;  % Poisson ratio
rho = 2650; % density [kg/m^3]

Lx = 1; % [m]
Ly = 1.4;   % [m]

c_L = sqrt( E/(rho*(1-(nu^2))) );

f = zeros(5,5);

for i=0:5
    for j=0:5
        f(i+1, j+1) = 0.453*c_L*h*( ((i+1).^2/Lx^2) + ((j+1).^2/Ly^2) );
    end
end
