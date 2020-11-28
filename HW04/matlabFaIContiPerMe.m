close all, clear all, clc;

syms rho c S1 L1 L2 x1 Z0 k f;
Zc = (i*rho*c)/S1*(cot(k*L1) + (k*x1)^(-1))^(-1);
Zin = Z0*( (Zc*cos(k*L2) + i*Z0*sin(k*L2)) / (i*Zc*sin(k*L2) + Z0*cos(k*L2)) );

vals = [1.25, 343, pi*(0.1)^2, 0.5305, 0.5, 0.5, 54590.1455];

Zin_val = subs(Zin, [rho, c, S1, L1, L2, x1, Z0], vals);


assume(k>0);
figure(1);
fplot((k*343)/(2*pi),abs(Zin_val), [0 15]);
ylim([0 5e5]);
hold on;
xlabel('f [Hz]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex');
ylabel('$|Z_{in}|$ [$\frac{Pa \times s}{m^3}$]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex');


solk = zeros(1,4);
solk(1) = vpasolve(Zin_val == 0, k, 3.5);
solk(2) = vpasolve(Zin_val == 0, k, 6.1);
solk(3) = vpasolve(Zin_val == 0, k, 9.6);
solk(4) = vpasolve(Zin_val == 0, k, 12.4);

%fMins = (343*solk)/(2*pi);
fMins = [1.952244289324675e+02,3.329979658736029e+02,5.103338241089317e+02,6.659933466564429e+02];
fMaxs = [1.068046999390690e+02,2.620908149216368e+02,4.212985338627673e+02,5.889605676211254e+02];

str = string(fMins);
for i=1:length(fMins)
    xline(fMins(i), '--', str(i) + ' Hz');
    xline(fMaxs(i));
end



