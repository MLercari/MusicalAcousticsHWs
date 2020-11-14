%% *Homework 3 - Musical Acoustics*
% *Data*
% Note that $\frac{\nu_{LR}}{E_L} \ne \frac{\nu_{RL}}{E_R}$

%density with clear wood correction
rho = 400; %[Kg/m3]
rho_correction = 0.1;
rho_corr = rho*(1 + rho_correction);

% Longitudinal Young modulus with no corrections
E_L = 10.8e9; %[Pa]

% Correction for clear wood
E_L_young_clear_correction = 0.22;

% Correction for bending stress
E_L_young_bend_correction = 0.1;

% Corrected longitudinal elastic modulus
E_L_corr = E_L*(1 + E_L_young_bend_correction + E_L_young_clear_correction);

% Radial elastic modulus
E_R = E_L_corr*0.078;

%Poisson's ratios
nu_RL = 0.040;
nu_LR = 0.372;
% Question 1 - cutoff frequency
% since 
% 
% $$v_p = c \\\sqrt{1.8 \  f \ h \ c_L} = c \\h = \frac{c^2}{1.8 \ f \ c_L}$$
% 
% and 
% 
% $$c_L = \sqrt{\frac{E_L}{\rho( 1 - \nu_{LR}\nu_{RL})}}$$

% speed of sound in air
c = 343; %[m/s]

%longitudinal speed of sound in wooden thin plate
c_l = sqrt(E_L_corr/(rho_corr*(1 - nu_LR*nu_RL))); %[m/s]

%cutoff frequency 
f_c = 1200; %[Hz]

%thickness of thin plate considering longitudinal waves
h = c^2/(c_l*1.8*f_c);

%% 
% result $h = 9.5  \text{ mm}$
% Question 2- wave propagation direction
% We know that
% 
% $$\tan \theta =   \sqrt { \left ( \frac{v_p^2}{c^2}-1 \bigg ) }\\v_p = \sqrt{f} 
% \sqrt{1.8 \ h \ c_L}$$
% 
% 
f = linspace(1000,10000,10000); %frequency domain

%plate transverse speed as a function of frequency
v_p = sqrt(1.8*h*c_l.*f); %[m/s]

%direction of propagation
tg_theta = sqrt(v_p.^2/c^2 - 1);

%plot

figure(1)
%semilogx(f, real( radtodeg(atan(tg_theta))), 'LineWidth', 2 , "color", "k");
plot(f/1000, real( radtodeg(atan(tg_theta))), 'LineWidth', 2 , "color", "k");
grid on;
xlabel(" f [kHz]")
ylabel(" \theta [deg] ")

