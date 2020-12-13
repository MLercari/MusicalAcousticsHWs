%% QUESTION 1 - BRIDGE IMPEDANCE IN THE TWO-MASS MODEL
clc, clear all, close all;

%% Parameters
k_p = 1.41e5; % [N/m] top plate stiffness
m_p = 0.128*0.385;    % [kg] top plate effective mass
A_p = 0.0375*0.385;   % [m^2] top plate area
R_p = 32; % [N s m^-1] top plate resistance

m_h = 8.04e-4;   % [kg]  air piston mass
A_h = 7.85e-3;   % [m^2] air piston area
R_h = 30;    % [Ns/m^5] sound hole resistance

V = 0.0172;  %[m^3] air cavity volume
rho = 1.2;   % [kg m^-3] air density
c = 343;    % [m/s] sound velocity in air

L_p = m_p/(A_p^2);  % plate inertance
L_h = m_h/(A_h^2);  % sound hole inertance
C_p = A_p^2/k_p;    % plate compliance
C_v = V/(rho*c^2);  % cavity compliance

%plot of the bridge impedance
Fs = 44100;

f = 0:1/Fs:500; %[Hz] range of frequency in Hz
w = 2*pi*f; %[rad/s] range of angular frequency 

Z_p = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1); %[Kg/m^4s] plate impedance
Z_v = 1./(1i.*w*C_v); %[Kg/m^4 s] cavity impedance
Z_h = 1i.*w*L_h + R_h; %[Kg/m^4 s] hole impedance

Z = Z_p + ( Z_v .* Z_h ) ./ (Z_v + Z_h ); %[Kg/m^4 s] bridge impedance

fH = sqrt(1.4*101000*A_h^2/(m_h*V))/(2*pi); %[Hz] Helmholtz resonance frequency

fp = sqrt((k_p + (1/C_v)*A_p^2)/m_p)/(2*pi); %[Hz] closed box bridge resonance frequency 
fa = sqrt(((1/C_v)*A_p^2)/m_p)/(2*pi); %[Hz] no stiffness in the sound board

fres = [1 0 (-fp^2 -fH^2) 0 (fp^(2)*fH^(2) - fa^(2)*fH^(2)) ];
fsol = roots(fres);


figure(1);
subplot 211;
semilogy(f, abs(Z), 'lineWidth' , 1);
grid on
xline(fH, '-.',{'fh = 126.3 Hz'}  , 'lineWidth' , 1);
xline(fp, '-.',{'fp = 270.8 Hz'}  , 'lineWidth' , 1);
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$|Z|$  [$\frac{Kg}{m^4 \times s}$]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')

subplot 212;
plot(f, angle(Z), 'lineWidth' , 1);
grid on
xline(fH, '-.',{'fh = 126.3 Hz'}  , 'lineWidth' , 1);
xline(fp, '-.',{'fp = 270.8 Hz'}  , 'lineWidth' , 1);
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$ \angle Z$  [$rad$]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')






