clear all
close all
clc

%% H6 Question One 

%data 
rho = 1.2;
c = 343;
F4 = 349.23;
k = (F4*2*pi)/c;
theta = deg2rad(0.75);
L = 0.45;

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

S1 = (x1value*tan(theta))^2*pi;
a1 = x1value*tan(theta);
%diameter of resonator head = 0.1252 m
d1 = a1*2; 
a2 = (x1value + L)*tan(theta);
S2 = (a2^2)*pi;
%diameter of resonator foot =  0.1369 m
d2 = a2*2;

%% H6 Question Two and Three

% FIRST APPROACH: simplify model because conical angle is very small: 
% cylindrical bore with section S1 and finger hole with area Sh = S1. 

% characteristic impedance of a cylindrical horn
Z0 = rho*c/S1; 
%end correction for this model
Delta = 0.6*a1; 
%new frequency
G4 = 392; 
k2 = G4*2*pi/c;
%assume the finger hole has the same area as the cylinder so that l = Delta
S = S1;

syms D
assume(D > 0)
Zcy = (1i*rho*G4*2*pi/S)*(L -  D - Delta^2/(D + 2*Delta)) == 0;

%imposing Zcy equal to zero in order to find D
 solD = vpasolve(Zcy, D, [0 inf]);
 Dvalue = double(solD);
%D is 0.4473 m -> no physical meaning 

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

%small finger hole: here the finger hole has the radius equal to L*12/2% = 0.054 (for a tone).
%In this way we are sure the finger hole is located within the recorder.
S = (L*0.12)^2*pi; 
%end correction for this model
Delta = 0.6*a2;
%the finger hole's acoustic length is not equal to Delta but proportional
%to the finger hole area
l = Delta*(S/S1);

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

%FOURTH APPROACH: small finger holes and cylindrical bore
Delta = 0.6*a1; 
syms D4
assume (D3 > 0 )
Zcy4 = 1i*Z0*tan(k2*L) - (-1i*(S/rho*c)*cot(k2*l) -1i*(S1/(rho*c))*cot(k2*(Delta + D4)))^(-1) == 0;
SolD4 = vpasolve(Zcy4, D4, [0 inf]);
D4Value = double(SolD4);

%% question 4
dp = 62;    %[Pa] pressure difference btw player's mouth and exit of flue channel
fc = 2e3;   %[Hz] target spectral centroid
U = sqrt((2*dp)/rho);   %[m/s] flow velocity
h = 0.3*(U/fc);     % [m] channel thickness

nu = 1.5e-5;    %[m^2/s] kinematic viscosity of air

Re = (U*h)/nu;  % Reynolds number

%% question 5
Lc = 0.2;   %[m] channel length
delta = sqrt((nu*Lc)/U);    %[m] boundary layer thickness at the exit

