clear all
close all
clc

%% H6 Question One 

%data 
rho = 1.2;
c = 343;
F4 = 349.23;    % frequency with all finger holes closed
k = (F4*2*pi)/c;
theta = deg2rad(0.75);
L = 0.45;   % [m] length of the resonator
l = 0.003;  % [m] wall thickness
S = 6.5e-5;   % [m^2] mouth cross section - arbitrary
M = (rho*l)/S;  %[kg m^-4] mouth inertance


%{
%the input impedance as a function of x1
syms x1 
Zin = (1i*rho*c)/(pi*(x1*tan(theta))^2)*(sin(k*(L + 0.6*(x1 + ...
L)*tan(theta)))*sin(atan(k*x1)))/(sin(k*(L + 0.6*(x1 + L)*tan(theta) + ...
(1/k)*atan(k*x1)))) == 0;

%imposing Zin equal to zero in order to find x1
[solx1 , param, cond] = solve(Zin, x1 , 'ReturnConditions', true);
assume(cond)
solk = solve(solx1 > 0, param);
x1value = subs(solx1, solk);
x1value = double(x1value);

% x1 = 4.7802 m
Lcorr = L + 0.6*(x1value + L)*tan(theta);
S1 = (x1value*tan(theta))^2*pi;
a1 = x1value*tan(theta);
%diameter of resonator head = 0.1252 m
d1 = a1*2; 
a2 = (x1value + L)*tan(theta);
S2 = (a2^2)*pi;
%diameter of resonator foot =  0.1369 m
d2 = a2*2;
%}

% WITH THE MOUTH IMPEDANCE
syms x1;
th1 = atan(k*x1);
r1 = x1*tan(theta);

%URGENTE PER MATTI: controllare questa formula non sembra giusta
eqn = ((rho*c)/(pi*r1^2))*sin(k*L)*sin(k*th1) + k*c*M*sin(k*(L+th1));   % condition Zp + Zm = 0



figure(33); % check zeros to set the initial condition
fplot(eqn);
hold on;
fplot(0, 'LineStyle', '--');


solQ1 = vpasolve(eqn==0, x1, 0.8);

x1 = double(solQ1);
x2 = x1 - L;

x1 = abs(x1);
x2 = abs(x2);   % take absolute values


a1 = x1*tan(theta);
a2 = x2*tan(theta);
S1 = pi*(a1)^2;
S2 = pi*(a2)^2;

DL = (M/rho)*S1;

d1 = 2*a1;
d2 = 2*a2;

x1 = abs(x1);
x2 = abs(x2);   % take absolute values

%% H6 Question Two and Three

Lcorr = L + 0.6*a2;

% FIRST APPROACH: simplify model because conical angle is very small: 
% cylindrical bore with section S1 and finger hole with area Sh = S1. 

% characteristic impedance of a cylindrical horn
Z0 = rho*c/S2; 
%end correction for this model
Delta = 0.6*a2; % 2 refers to the foot
%new frequency
G4 = 392; 
k2 = G4*2*pi/c;
%assume the finger hole has the same area as the cylinder so that l = Delta
S = S2;

syms D
assume(D > 0)
%Zcy = (1i*rho*G4*2*pi/S)*(L -  D - Delta^2/(D + 2*Delta)) == 0;
d = D + Delta^2/(D + 2*Delta);
Zcy = 1i*Z0*tan(k2*(Lcorr - d + DL)) == 0;
%imposing Zcy equal to zero in order to find D
 solD = vpasolve(Zcy, D, 0.05);
 Dvalue = double(solD);

%FOURTH APPROACH: small finger holes and cylindrical bore
Delta = 0.6*a2; 
%S = (L*0.12/29)^2*pi;
S = pi*(0.003)^2;
%the finger hole's acoustic length is not equal to Delta chosen as 0.015
%btw 15 mm was def too much for the thickness of the walls so i set it to 3 mm -Mattia

syms D4;
assume (D4 > 0 )
% Yhole = -1i*(S/(rho*c))*cot(k2*l);
% Ypipe = -1i*(S1/(rho*c))*cot(k2*(Delta + D4));
d4 = S*(D4 + Delta)^2/(S*(D4 + Delta) + S2*l);
%Zcy4 = 1i*Z0*tan(k2*(Lcorr -d4 + DL));
eqq = k2*(Lcorr -d4 + DL);
SolD4 = vpasolve(eqq==pi, D4, 0.05);
D4Value = double(SolD4);

%{
%SECOND APPROACH: conical bore impedance minus finger hole's impedance and
% cylindrical bore impedance

%here we assume the finger hole area is the same as the head area S = S1
syms D2
assume (D2 > 0 )
Zin2 = 1i*(rho*c)/(S1) * (sin(k2*L)*sin(atan(k2*x1value)))/(sin(k2*(L + (1/k2)*atan(k2*x1value)))) ...
    - 1i*2*pi*G4*rho*(D2 + (Delta^2/(D2 + 2*Delta)))/S == 0;

%imposing Zin2 equal to zero in order to find D2
SolD2 = vpasolve(Zin2, D2, [0 inf]);
D2Value = double(SolD2);
%NO SOLUTION


%THIRD APPROACH: conical shape and small finger hole

%small finger hole: here the finger hole has the radius equal to L*12/4% .
%In this way we are sure the finger hole is located within the recorder.
S = (L*0.12/10)^2*pi; 
%end correction for this model
Delta = 0.6*a2;
%the finger hole's acoustic length is not equal to Delta chosen as 0.015
l = 0.015;

syms D3
assume (D3 > 0 )
Zin3 = 1i*(rho*c)/(S1) * (sin(k2*L)*sin(atan(k2*x1value)))/(sin(k2*(L + (1/k2)*atan(k2*x1value)))) - ...
     ( -1i*(S/(rho*c))*cot(k2*l) - ... %admittance of the hole and conical end correction
    1i*(((x1value + L - D3)*tan(theta))^2*pi)/(rho*c) *(sin(k2*(Delta + D3 + ... 
    (1/k2)*atan(k2*(x1value + L - D3))))/(sin(k2*(Delta + D3))*sin(atan(k2*(x1value + L - D3))))))^(-1) == 0;

%imposing Zin3 equal to zero in order to find D2

SolD3 = vpasolve(Zin3, D3, [0 inf]);
D3Value = double(SolD3);
%VALUE WITH NO PHYSICAL MEANING


%} 
%% question 3

k3 = 440*2*pi/c;
D = D4Value(2);
syms dist
assume (dist > 0 )
Delta_prime = (S2*l*(D + Delta))/(S*(D + Delta) + S2*l);
Delta_sec = l*(dist/S2 + Delta_prime/S2)/(dist/S2 + Delta_prime/S2 + l/S);
d3 = D + dist + Delta - Delta_sec;
Z3 = 1i*Z0*tan(k3*(Lcorr - d3 + DL));
Soldist = vpasolve(Z3 == 0, dist, 0.12*(L-D));
distValue = double(Soldist);

expr = k3*(Lcorr - d3 + DL) - pi;


%% question 4
dp = 62;    %[Pa] pressure difference btw player's mouth and exit of flue channel
fc = 2e3;   %[Hz] target spectral centroid
U = sqrt((2*dp)/rho);   %[m/s] flow velocity
h = 0.3*(U/fc);     % [m] channel thickness

nu = 1.5e-5;    %[m^2/s] kinematic viscosity of air

Re = (U*h)/nu;  % Reynolds number

%% question 5
Lc = 20e-3;   %[m] channel length
delta = sqrt((nu*Lc)/U);    %[m] boundary layer thickness at the exit


%% question 6
Str = (F4*h)/U;
alfa = 0.35/h;   % [m^-1] de la Cuadra, p. 32
Pout = 10^(9/2)*2e-5; %[Pa]
Vac = U/10;
delj = Lc/sqrt(Re);
delAc = sqrt(2*nu/(F4*2*pi));
eta0 = (Vac*h*delj)/(U*delAc);

H = 0.015;  %[m] jet width = mouth width (CORRECT Q1!!)
b = 0.45*h;
syms W;
eta = eta0*exp(alfa*W);
Qin = b*H*U*(1 + tanh((eta)/b));
pippo = Qin - Pout/((1i*(rho*c)/S1)*tan(k*Lcorr));

figure(77);
fplot(abs(pippo));
hold on;
fplot(0);

solutionz = vpasolve(pippo==0, W, 0);








