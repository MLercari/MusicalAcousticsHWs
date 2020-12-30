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
fs = 4*44.1e3;
ts = 1/fs;

% Fundamental note
f1 = 52.8221; %[Hz] frequency of C2

% Boundary          

% String parameters
b1 = 0.5;
b2 = 6.25e-9;
T0 = 750;  % [N] tension at rest
L = 1.92;    % [m] length of the string
M_s = 35e-3;    % [kg] mass of the string
mu = M_s/L;   % [kg m^-1] avg linear mass density of the string
eps = 7.5e-6;
c = sqrt(T0/mu);    % [m/s] string wave velocity

% Spatial sampling parameters
X_max = sqrt( 0.5*((c*ts)^2 + 4*b2*ts + sqrt( ((c*ts)^2 + 4*b2*ts)^2 + 16*(eps*ts)^2 )) );
gamma = fs/(2*f1);

% Aliasing condition
% Number of maximum spatial steps
N_max = ceil(sqrt( (-1 + sqrt(1 + 16*eps*(gamma^2)))/(8*eps) ));
X_min = L/N_max;

% Integer values
% Spatial sampling
N = N_max - 1;
X = L/N;    %as close as possible to Xmin (see Chaigne et al.)

if (X < X_min | X > X_max )
    display("WARNING: WRONG SPATIAL SAMPLING STEP");
end

% FD parameters

% Hammer parameters
Vh0 = 2.5;  % [m/s] initial hammer velocity

% Hammer contact window definition



%PDE Coefficients:

% Bridge boundary coefficients

% Left hand (hinged string end) boundary coefficients

% Hammer felt parameters

%% Computation of the FD scheme
% Initialization

% Computation loop
%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










