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
solk = solve(solx1 > 0, param)
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


