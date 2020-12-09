%% QUESTION 1 - BRIDGE IMPEDANCE IN THE TWO-MASS MODEL
clc, clear all, close all;

%% Parameters
k_p = 1.41e5 % [N/m] top plate stiffness
m_p = 0.128*0.385    % [kg] top plate effective mass
A_p = 0.0375*0.385   % [m^2] top plate area
R_p = 32 % [N s m^-1] top plate resistance

m_h = 8.04e-4   % [kg]  air piston mass
A_h = 7.85e-3   % [m^2] air piston area
R_h = 30    % [Ns/m^5] sound hole resistance

V = 0.0172  %[m^3] air cavity volume
rho = 1.2;   % [kg m^-3] air density
c = 343;    % [m/s] sound velocity in air

L_p = m_p/(A_p^2);  % plate inertance
L_h = m_h/(A_h^2);  % sound hole inertance
C_p = A_p^2/k_p;    % plate compliance
C_v = V/(rho*c^2);  % cavity compliance

%syms w;
f = 0:0.5:500;
w = 2*pi*f;
%s = 1i*w;
%{
Z = ( s*L_h + R_h ) / ( L_h*C_v*s^2 + R_h*C_v*s );
Z = Z + ( L_p*C_p*s^2 + R_p*C_p*s +1 ) / ( s*C_p );
%}

%Z = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1) + ( 1i.*w*C_v + (1i.*w*L_h + R_h).^(-1) ).^(-1);

Z_p = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1);
Z_v = 1./(1i.*w*C_v);
Z_h = 1i.*w*L_h + R_h;

Z = Z_p + ( Z_v .* Z_h ) ./ (Z_p + Z_h );

figure(1);
subplot 211;
semilogy(f, abs(1./Z));
subplot 212;
plot(f, angle(1./Z));

