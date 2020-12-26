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
%diameter of resonator foot =  0.1369 m
d2 = a2*2;

%% H6 Question Two and Three

% first approach: simplify model because conical angle is very small: 
% cylindrical bore with section S1 and finger hole with area Sh = S1. 

% characteristic impedance of a cylindrical horn
Z0 = rho*c/S1; 
%end correction for this model
Delta = 0.6*a1; 
%new frequency
G4 = 392; 
k2 = G4*2*pi/c;
S = S1;

syms D
assume(D > 0)
Zcy = (1i*rho*G4*2*pi/S)*(L -  D - Delta^2/(D + 2*Delta))== 0;

%imposing Zcy equal to zero in order to find D
 solD = solve(Zcy, D);
 Dvalue = double(solD);
%D is 0.4473 m

%second approach: conical bore with section from finger hole to 
% foot constant and finger hole's section Sh = S1

syms D2
assume (D2 > 0 )
Zin2 = 1i*(rho*c)/(S1)*(sin(k2*(L - D2 - (Delta^2)/(D2 +  ...
2*Delta)))*sin(atan(k2*x1value)))/sin(k2*(L - D2 - (Delta^2)/(D2 +  ...
2*Delta) + (1/k2)*atan(k2*x1value))) == 0;

%imposing Zin2 equal to zero in order to find D2
SolD2 = solve(Zin2, D2);
D2Value = double(SolD2);
%D is 0.4473 m -> very strange!

%third approach: conical shape and small finger hole (to be continued!)







