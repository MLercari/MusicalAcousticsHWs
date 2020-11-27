clear all, close all, clc;
%% Homework 4 - horns and pipes
% Data 

a1 = 0.05; %[m] radius of pipe 
a2 = 0.1; %[m] radius of conical horn 
L1 = 0.5; %[m] length of pipe
L2 = 0.5; % [m] length of conical horn
c = 343; %[m/s] speed of sound in air
rho = 1.25; %[kg/m2] density of air
S1 = pi*(a1)^2; %[m^2] surface of pipe section
S2 = pi*(a2)^2; %[m^2] surface of cone section (ending)
Z01 = rho*c/S1; %characteristic impedance of pipe

f = linspace(0 , 4000, 4000); %[Hz] frequency range 
omega = f.*(2*pi); %[rad/s] frequency range
k = omega./c; %[1/m] wave number
ka = k.*a1; %[-] adimensional ka for pipe only

%% question 1 Determine the frequencies of the maxima impedance
% of the cylindrical pipe only, considering the
%presence of the radiation load.

%first of all, in order to take into account radiation losses, we add a 
% virtual elongation 
%actually, DeltaL1 is a function of ka: see fig. 8.9 of Fletcher(pag. 201)
% let's take a mean value of Delta/a between 0.1 and 0.61

DeltaL1 = 0.61*a1; %[m]

L1corr = L1 + DeltaL1; 

%for an open, not flanged end pipe the input impedance is then
Zin1 = 1i*Z01*tan(k.*L1corr);

%where maxima are at
n = [1:10];

omega_max_1 = (n - 0.5).*(pi*c/L1corr);

%let's take a look at 
figure(1)
plot( f, imag(Zin1),'c', 'lineWidth' , 0.5);
hold on 
plot ( f, imag(1i*Z01*tan(k.*L1)), 'm', 'lineWidth' , 0.5);
xline(162,'-.',{'161.6 Hz'}  , 'lineWidth' , 1.5);
xline(485,'-.',{'484.9 Hz'}  , 'lineWidth' , 1.5);
xline(808,'-.',{'808.2 Hz'}  , 'lineWidth' , 1.5);
xline(1132,'-.',{'1131.5 Hz'}  , 'lineWidth' , 1.5);
legend(["with radiation load", "without radiation load"])
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$Z_{in}$ open pipe [$\frac{Pa \times s}{m^3}$]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylim([ -1e6 1e6]);
xlim([0 1500]);

%% Question 2 Consider now the compound horn. Determine the frequencies 
%of the first four maxima of the input impedance.

%we need to solve the follorwing equation
%the distance from the narrow end and the conic apex is x1

x1 = 0.5; 

syms x
figure(2);
fplot(tan(x*L1corr) - cot(x*L2) - (1/(x*x1)));

S = zeros(1, 4);
S(1) = vpasolve(tan(x*L1corr) - cot(x*L2) - (1/(x*x1)) == 0, x, 1.5); % NB: remember to fplot the lhs of the equation to check for the zeros
S(2) = vpasolve(tan(x*L1corr) - cot(x*L2) - (1/(x*x1)) == 0, x, 4.5);
S(3) = vpasolve(tan(x*L1corr) - cot(x*L2) - (1/(x*x1)) == 0, x, 7.7);
S(4) = vpasolve(tan(x*L1corr) - cot(x*L2) - (1/(x*x1)) == 0, x, 11);

f_max_2 = (S*c)/(2*pi);

figure(2)

plot( f, tan(k.*L1corr) - cot(k.*L2) - (1./(k.*x1)),'b', 'lineWidth' , 0.5);
hold on 
xline(107,'-.',{'106.8 Hz'}  , 'lineWidth' , 1.5);
xline(262,'-.',{'262.1 Hz'}  , 'lineWidth' , 1.5);
xline(422,'-.',{'421.3 Hz'}  , 'lineWidth' , 1.5);
xline(589,'-.',{'589 Hz'}  , 'lineWidth' , 1.5);
yline(0)
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$\tan(kL_2) - \cot(kL_1) - 1/kx_1$" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylim([ -40 40]);
xlim([0 700]);
% syms x 
% eqn = tan(x*L1) - cot(x*L2) - (1/(x*x1)) == 0;
% for n = 1:4 
%   vpasolve(eqn ,x ,'Random',true) 
% end
% %k_max_2 = vpasolve(eqn, k); 

%% Question 3 - Determine the frequencies of the first four minima of the input impedance
%when the numerator of Zin goes to zero (see below)

%% Question 4 - Plot the impedance function in the range [0Hz, 4 kHz]

% compute separately the conical input impedance
%assuming the radiation load with the same correction of the pipe

Zc = (1i*rho*c/S1).*(1./(cot(k.*L1corr) + 1./(k.*x1)));

%then we plug into the inpute impedance of the pipe (with no radiation
%losses)

Zin2 = Z01.*((Zc.*cos(k.*L1) + 1i*Z01*sin(k.*L1))./(1i.*Zc.*sin(k.*L1) + Z01.*cos(k.*L1)));

figure(3)
plot(f, abs(Zin2),'k', 'lineWidth' , 0.7);
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$|Z_{in}|$" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylim([0 1e6]);
xlim([0 4000]);


syms kappa;
ZcSym = (1i*rho*c/S1)*(1/(cot(kappa*L1corr) + 1/(kappa*x1)));
ZSym = Z01*((ZcSym*cos(kappa*L1) + 1i*Z01*sin(kappa*L1))/(1i*ZcSym*sin(kappa*L1) + Z01*cos(kappa*L1)));
min = vpasolve(ZSym == 0, kappa, (2*pi*200)/c);
min = (min*c)/(2*pi);


 %ka = a2pi f / c