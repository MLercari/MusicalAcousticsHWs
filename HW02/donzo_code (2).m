clear all;
close all;
clc
%% Homework 2: Soundboard modeling and string coupling %%
% DATA
Lx = 1; %[m]
Ly = 1.4; %[m]
h = 0.003; %[m]
E = 69*10^9; %[N/m^2]
ni = 0.35;
ro = 2650; %[kg/m^2]

%% Modal frequencies of the lowest five modes

% Supported edges
m = 3;
n = 3;

f_s = zeros(1,5);
c_l = sqrt(E/(ro*(1-ni^2)));

k = 1;
for i = 0:m
    for j = 0:n
    f_s(k) = 0.453 * c_l * h * ((i+1)^2/Lx^2 + (j+1)^2/Ly^2);
    k = k+1;
    end
end
disp(f_s);

%% Free edges
ni1 = 1/6;
D1 = (E*h^3)/(12*(1-ni1^2));
D2 = (E*h^3)/(12*(1-ni^2));

ro_nasa = ro*h; %[kg/m^2]

nasa_coeffs_f = [9.905 9.944 22.245 22.373 27.410];

w_nasa_f = zeros(1,5);

w_nasa_f = nasa_coeffs_f./(Lx^2*sqrt(ro_nasa/D1));

f_f = w_nasa_f/2/pi;

%% Clamped edges

nasa_coeffs_c = [28.0502 45.7196 67.005 75.0586 83.390];

w_nasa_c = zeros(1,5);

w_nasa_c = nasa_coeffs_c./(Lx^2*sqrt(ro_nasa/D2));

f_c = w_nasa_c/2/pi;

%% Q3
t = 0:1/100000:1/1000;
alpha = 1000;

%f_f = [6.96 8.06 16.44 16.7 20.64];
f_s = [11.17 22.5106 33.3883 41.395 44.7191];
%f_c = [21.37 34.85 51.47 57.27 64.71];

fp_f = f_f .* (1 + 0.16 * (0.0004*(exp(-alpha *[t' t' t' t' t']).*exp(1i*2*pi*f_f.*t')/0.004).^2));
fp_c = f_c .* (1 + 0.16 * (0.0004*(exp(-alpha *[t' t' t' t' t']).*exp(1i*2*pi*f_c.*t')/0.004).^2));
fp_s = f_s .* (1 + 0.16 * (0.0004*(exp(-alpha *[t' t' t' t' t']).*exp(1i*2*pi*f_s.*t')/0.004).^2));

figure (1)
for i = 1:5
    subplot(5,1,i)
    plot(t, real(fp_f(:,i)));
    ylabel(['f_',num2str(i),' [Hz]']);
    xlabel('time [ms]');
end

figure (2)
for i = 1:5
    subplot(5,1,i)
    plot(t, real(fp_c(:,i)));
    ylabel(['f_',num2str(i),' [Hz]']);
    xlabel('time [ms]')
end
 
figure (3)
for i = 1:5
    subplot(5,1,i)
    plot(t, real(fp_s(:,i)));
    ylabel(['f_',num2str(i),' [Hz]']);
    xlabel('time [ms]');
end

%% Q4
a = 0.05;
x = 0:0.01:1;
y = 0:0.01:1.4;

Z = a.*sin(4*pi*x).*sin(2*pi*y'./Ly)+ a.*sin(3*pi*x).*sin(4*pi*y'./Ly); %funzione velcità moltilplicata per meshgrid
vel = 1i*2*pi*130.8*Z;

v_a = abs(vel);
M = max(v_a, [], 'all');
[y1,x1]=find(v_a >= M-1e-6);

y1 = y1*0.01;
x1 = x1*0.01;

figure(4)
surf(x,y,Z);  %plot
hold on
for i=1:length(x1)
    plot3(x1(i), y1(i), 0.1,'r*')
end
hold off
xlabel('Asse X');
ylabel('Asse Y');
zlabel('Asse Z');

%% 5) --- spostamento bridge---
C3 = 130.8; %Hz
L = 1.007; %m
mu = 12e-3; %kg/m
eps = -0.01:0.001:0.01; %[0.1 3]

T1=4*L^2*C3^2*mu;
T2=T1*(1+4*eps);
Z0=sqrt(mu*T1);

dyb = a*2*pi*cos(2*pi*x1(1))*sin(4*pi*y1(1)/Ly) + a*4*pi*cos(4*pi*x1(1))*sin(3*pi*y1(1)/Ly);
fb = -T1*dyb - T2*dyb;

Yb = M./fb;

chi=1i*Z0*Yb/pi;
ap=chi+eps+sqrt(eps.^2+chi.^2);
am=chi+eps-sqrt(eps.^2+chi.^2);

figure(5)
plot(2*pi*C3*eps,2*pi*C3*ap,2*pi*C3*eps,2*pi*C3*am)

%% --- spostamento corda ---
eps = linspace(-0.03,0.03,1000); %[0.1 3]

T1=4*L^2*C3^2*mu;
T2=T1*(1+4*eps);
Z0=sqrt(mu*T1);

fb= T1*pi/L + T2*pi/L;

Yb= M./fb;

chi=1i*Z0*Yb/pi;
ap=chi+eps+sqrt(eps.^2+chi.^2);
am=chi+eps-sqrt(eps.^2+chi.^2);

figure(6)
plot(2*pi*C3*eps,2*pi*C3*real(ap),2*pi*C3*eps,2*pi*C3*real(am))
xlabel('f_0 x \epsilon');
ylabel('f_0 x Re[a_\pm]')

figure(7)
plot(2*pi*C3*eps,2*pi*C3*imag(ap),2*pi*C3*eps,2*pi*C3*imag(am)) 
xlabel('f_0 x \epsilon');
ylabel('f_0 x Im[a_\pm]')




