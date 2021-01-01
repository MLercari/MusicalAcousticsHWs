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
fs = 4*44.1e3; %[Hz] sampling frequency
ts = 1/fs;
duration = 8; %[s] signal length
tSamples = duration*fs;
t = linspace(0, duration, tSamples); %temporal vector

% Fundamental note
f1 = 52.8221; %[Hz] frequency of C2

% Boundary          

% String parameters
b1 = 0.5; %air damping coefficient
b2 = 6.25e-9; %internal friction of the string coefficient
T0 = 750;  % [N] tension at rest
L = 1.92;    % [m] length of the string
M_s = 35e-3;    % [kg] mass of the string
mu = M_s/L;   % [kg m^-1] avg linear mass density of the string
eps = 7.5e-6; % [-] string stiffness parameter assumed equal to kappa- string stiff. coefficient
c = sqrt(T0/mu);    % [m/s] string wave velocity
x0 = 0.12*L; %[m] stricking position (see d/L from Fletcher)

% Aliasing condition
if (fs < 2*f1)
    disp("WARNING: SAMPLING FREQUENCY BELOW NYQUIST");
end

% Spatial sampling parameters
X_max = sqrt( 0.5*((c*ts)^2 + 4*b2*ts + sqrt( ((c*ts)^2 + 4*b2*ts)^2 + 16*(eps*ts)^2 )) ); %stability condition (see Saitis)

% Number of maximum spatial steps
gamma = fs/(2*f1);

% Integer values
N_max = ceil(sqrt( (-1 + sqrt(1 + 16*eps*(gamma^2)))/(8*eps) ));
X_min = L/N_max;

% Spatial sampling
N = N_max - 1;
X = L/N;    %as close as possible to Xmin (see Chaigne et al.)

if (X < X_min | X > X_max )
    display("WARNING: WRONG SPATIAL SAMPLING STEP");
end

x = linspace( 0 , L , X); %spatial vector

% FD parameters
lambda = c*ts/X; %Courant parameter
nu = 2*b2*ts/(T^2); %convenient variable n. 1
mi = eps^2/(c^2*X^2); %convenient variable n.2

% Hammer parameters
Vh0 = 2.5;  % [m/s] initial hammer velocity
af = (ts^2/mu)/(1 + b1*ts); %force coefficient
bh = 1e-4; %[s^-1] fluid damping coefficient
Mh = 4.9e-3; %[Kg] hammer mass
eta = zeros(length(t)); %hammer displacement vector
Fh = zeros(length(t)); %hammer force vector

% Hammer contact window definition

%PDE Coefficients:

% ai coefficients
a1 = (-lambda^2*mi)/(1 + b1*ts);
a2 = (lambda^2 + 4*lambda^2*mi + nu)/(1 + b1*ts);
a3 = (2 - 2*lambda^2 - 6*lambda^2*mi - 2*nu)/(1 + b1*ts);
a4 = (-1 + b1*ts + 2*nu)/(1 + b1*ts);
a5 = -nu/(1 + b1*ts);

%di coefficients
d1 = 2/(1 + bh*ts/(2*Mh));
d2 = (-1 + bh*ts/(2*Mh))/(1 + bh*T/(2*Mh));
df = (-ts^2/Mh)/(1+bh*ts/(2*Mh));

% Bridge boundary coefficients
zeta_b = 1000; %normalized impedance at the bridge

%bri coefficients
br1 = (2 - 2*lambda^2*mi - 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br2 = (4*lambda^2*mi + 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br3 = (-2*lambda^2*mi)/(1 + b1*ts + zeta_b*lambda);
br4 = (-1 -b1*ts + zeta_b*lambda)/(1 + b1*ts + zeta_b*lambda);
br5 = (ts^2/mu)/(1 + b1*ts + zeta_b*lambda);

% Left hand (hinged string end) boundary coefficients
zeta_l = 1e20; %normalized impedance

%bli coefficients
bl1 = (2 - 2*lambda^2*mi - 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl2 = (4*lambda^2*mi + 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl3 = (-2*lambda^2*mi)/(1 + b1*ts + zeta_l*lambda);
bl4 = (-1 -b1*ts + zeta_l*lambda)/(1 + b1*ts + zeta_l*lambda);
bl5 = (ts^2/mu)/(1 + b1*ts + zeta_l*lambda);

% Hammer felt parameters
p = 2.3; 
a = 0.12; 
K = 4e8; 

%% Computation of the FD scheme
% Initialization
y = zeros(length(t), length(x));

% Computation loop
for n = 1:length(x)
    %boundary conditions
    
    for m = 1:length(t)
    %initial conditions
        
    end
end

%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










