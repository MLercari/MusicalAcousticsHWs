close all, clear all, clc;

syms rho c S1 L1 L2 x1 Z0 k;
Zc = (i*rho*c)/S1*(cot(k*L1) + (k*x1)^(-1))^(-1);
Zin = Z0*( (Zc*cos(k*L2) + i*Z0*sin(k*L2)) / (i*Zc*sin(k*L2) + Z0*cos(k*L2)) );

vals = [1.25, 343, pi*(0.1)^2, 0.5305, 0.5, 0.5, 54590.1455];

Zin_val = subs(Zin, [rho, c, S1, L1, L2, x1, Z0], vals);


assume(k>0);
figure(1);
fplot(abs(Zin_val));
hold on;


solk = zeros(1,4);
solk(1) = vpasolve(Zin_val == 0, k, 3.5);
solk(2) = vpasolve(Zin_val == 0, k, 6.1);
solk(3) = vpasolve(Zin_val == 0, k, 9.6);
solk(4) = vpasolve(Zin_val == 0, k, 12.4);

fMins = (343*solk)/(2*pi);

xline(solk(1));

