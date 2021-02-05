%% Extra of HW7
close all; clear all; clc;

%first and second transversal modes [Hz]
f = [8401.5, 8152.3, 7876, 7580.3, 7265.5, 6931.9, 6585.2, 6228.7, 5855.7];
s = [15489, 15357, 15177, 14963, 14707, 14409, 14076, 13715, 13317];

%first and second transversal modes [Hz] for thickness = 1.2 cm
f_thin = [4934.2, 4555.8, 4146.7, 3714.7, 3275.7, 2821.7, 2375, 1885.9, 1392.6];
s_thin = [11732, 11395, 10969, 10458, 9879.4, 9202.7, 8421.9, 7509.3, 6463.2];

eigen_thin = zeros(5,9);
eigen_thin(1,:) = f_thin;
eigen_thin(2,:) = s_thin;
eigen_thin(3,:) = [19246, 19040, 18758, 18367, 17891, 17224, 16357, 15196, 13791 ];
eigen_thin(4,:) = [26682, 26532, 26357, 26174, 26021, 25805, 25642, 24405, 23659];
eigen_thin(5,:) = [34056, 33852, 33589, 33304, 33077, 32796, 32408, 31843, 31093];


bendFreqs = [
    4.9342    4.5558    4.1467    3.7147    3.2757    2.8217    2.3750    1.8859    1.3926;
   11.7320   11.3950   10.9690   10.4580    9.8794    9.2027    8.4219    7.5093    6.4632;
   19.2460   19.0400   18.7580   18.3670   17.8910   17.2240    16.357   15.1960   13.7910;
   26.6820   26.5320   26.3570   26.1740   26.0210   25.8050   25.6420   24.4050   23.6590;
   34.0560   33.8520   33.5890   33.3040   33.0770   32.7960   32.4080   31.8430   31.0930
   ];

a = 0.001:0.001:0.009;


figure;
tiledlayout("flow");
hold on;
for i = 1:5
    plot(a*1e3, bendFreqs(i,:), 'Marker', 'd', 'LineStyle', '--', 'LineWidth', 1.0, 'MarkerFaceColor', 'auto');
end
xlabel('a [mm]');
ylabel('f [kHz]');

legend('f_1', 'f_2', 'f_3', 'f_4', 'f_5');


m = 1:1e4;
mn = zeros(5, 9);
for j = 1:9
    for i = 2:5
        [val, ind] = min(abs(bendFreqs(i, j) - m*bendFreqs(i-1, j)));
        mn(i, j) = ind;
    end
end

figure
tiledlayout('flow');
for N = 2:5
    I = zeros(1, 9);
    for j=1:9
        for i = 2:N
            I(1, j) = I(1, j) + abs( bendFreqs(i,j)/bendFreqs(i-1,j) - mn(i, j) );
        end
    end
    nexttile;
    bar(a*1e3, I, 0.4, 'FaceColor', [0.8500 0.3250 0.0980]);
    xlabel('a [mm]');
    ylabel(strcat('I_', num2str(N)));
end
   

   
   
   