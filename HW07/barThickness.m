%% Script for generating the thickness of a marimba bar as a function of the horizontal position
% Fabio Antonacci, 2020 @ Musical Acoustics, A.A. 2020/2021

close all, clear all, clc;

dx = 0.01; % Step of the horizontal axis
x = 0:dx:10-dx; % Horizontal axis
N = length(x);
z = zeros(size(x)); % Variable for storing the thickness of the bar as a function of the horizontal position
a = 0.1; % To be modified in the range 0.1-->0.9, step 0.1
b = 1;
for n = 1:floor(N/4) % Constant thickness for the first fourth of the length of the bar
z(n) = -2.5;
end

for n = 3*floor(N/4)+1:N % Constant thickness for the last fourth of the length of the bar
z(n) = -2.5;
end


figure
tiledlayout("flow");
for m = 1:9
    a = m*0.1;
    for n = floor(N/4)+1:3*floor(N/4) % Sinusoidal thickness profile in the central part of the bar
    z(n) = -2.5+a*cos((x(n)-5)*2*pi*b/10);
    end
    nexttile
    plot(x,z)
    grid on, grid minor;
    xlim([0 10]);
    ylim([-4 0]);
    xlabel('x [cm]');
    ylabel('y [cm]');
end


%% INHARMONICITY
close all, clc;  % uncomment when running the section multiple times;
M = readmatrix('eigenfrequencies.txt', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
M = M(:, [1 2]);
ind = abs(M(:, 2)) > 100;
M = M(ind, :);

a = unique(M(:,1));
eigfreqs = zeros(29, length(a));

for i = 1:length(a)
    eigfreqs(:, i) = M(M(:,1)==a(i), 2);
end
eigfreqs = eigfreqs(1:5, :);

eigfreqs = normalize(eigfreqs, 'scale', 'first');

figure;
tiledlayout("flow");
hold on;
for i = 1:5
    plot(a*1e3, eigfreqs(i,:), 'Marker', 'd', 'LineStyle', '--', 'LineWidth', 1.0, 'MarkerFaceColor', 'auto');
end
xlabel('a [mm]');
ylabel('f [Hz]');

legend('f_1', 'f_2', 'f_3', 'f_4', 'f_5');


m = 1:1e4;
mn = zeros(5, 9);
for j = 1:9
    for i = 2:5
        [val, ind] = min(abs(eigfreqs(i, j) - m*eigfreqs(i-1, j)));
        mn(i, j) = ind;
    end
end


figure
tiledlayout('flow');
for N = 2:5
    I = zeros(1, 9);
    for j=1:9
        for i = 2:N
            I(1, j) = I(1, j) + abs( eigfreqs(i,j)/eigfreqs(i-1,j) - mn(i, j) );
        end
    end
    nexttile;
    bar(a*1e3, I, 0.4, 'FaceColor', [0.8500 0.3250 0.0980]);
    xlabel('a [mm]');
    ylabel(strcat('I_', num2str(N)));
end

display(I);
%{
bar(a*1e3, I, 0.4, 'FaceColor', [0.8500 0.3250 0.0980]);
xlabel('a [mm]');
ylabel(strcat('I_', num2str(N)));
%}


        
        
        
