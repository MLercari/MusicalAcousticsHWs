%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 1
% Recording inspection
% We inspect the recordings for the estimation of the radiance pattern.
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc 

%% Setup

addpath('Functions');
addpath('Recordings');
Fs = 48000; %recording sampling frequency

typeOfSignal = ["sweep", "noise"];       % Noise or sweep

dirRec = ["Recordings/sweep", "Recordings/noise" ];  % Recordings directories

angles = zeros(1,24);

for i=2:length(angles) %compute the vector that contains the angles (in degree)
    angles(i) = 15 + angles(i-1);
end

nMic = 24;                        % Number of microphones
signalEnergy = zeros(2, 24);     % Vector of energy

% Labels will be the angle: store important information of each measument
% into a struct called structMis

for j = 1:length(dirRec) 
    
    if(dirRec(j) == "Recordings/sweep")
        
        Scont = dir(dirRec(j)); %struct with content of the directory
        Scont = Scont(4:end); %remove useless components
        
        %sort the struct directory in ascendent order
        ST = struct2table(Scont); % convert the struct array to a table to perform the sorting
        
        for k = 1:length(Scont)
          ST.name{k} = erase(ST.name{k},'.wav'); %erase the last string part .wav  
          ST.name{k} = str2double(ST.name{k}); %convert the remaining into a number
        end

        sortedST = sortrows(ST, 'name'); %sort the name field of the table
        sortedSCont = table2struct(sortedST); % change it back to struct array
        
        for ii = 1:height(sortedST) %assign name, type and angle and store to a struct
            structMis.fileName{ii} = sortedST.name{ii};
            structMis.angle{ii} = angles(ii);
            structMis.type{ii} = "sweep";
        end
   
    elseif(dirRec(j) == "Recordings/noise") %same for the noise recordings
        Ncont = dir(dirRec(j)); 
        Ncont = Ncont(4:end); 

        NT = struct2table(Ncont); 
        
        for k = 1:length(Ncont)
          NT.name{k} = erase(NT.name{k},'.wav');   
          NT.name{k} = str2double(NT.name{k}); 
        end
        
        sortedNT = sortrows(NT, 'name'); 
        sortedNCont = table2struct(sortedNT);  
        
        for ii = 25:height(sortedNT) + 24 %assign name, type and angle to the struct
            structMis.fileName{ii} = sortedNT.name{ii-24};
            structMis.angle{ii} = angles(ii-24);
            structMis.type{ii} = "noise";
        end

    end
    
end

%%
%in order to find the information for each signal, including the label,
%look into the structure structMis. To access one measurement information
%use:

MesNum = 27; %Sweep measurements are from 1 to 24 and and noise ones from 25 to 48
fn = structMis.fileName(MesNum);
fa = structMis.angle(MesNum);
ft = structMis.type(MesNum);
disp(strcat("the file " , num2str(fn{1}) , ".wav" , " is a " , ft{1} , ...
    " signal with an angle of " , num2str(fa{1}) , " degrees."));

%% Load the signals and compute the energy - compute_energy a function added in "Functions"

for jj = 1:1:height(sortedST)
    structMis.audio{jj} = audioread(strcat(dirRec(1),'/', num2str(structMis.fileName{jj}), '.wav'));
    structMis.energy{jj} = compute_energy( structMis.audio{jj} );
    signalEnergy(1, jj) = structMis.energy{jj};
end

for jj = 25:height(sortedNT) + 24
    structMis.audio{jj} = audioread(strcat(dirRec(2),'/', num2str(structMis.fileName{jj}), '.wav'));
    structMis.energy{jj} = compute_energy( structMis.audio{jj} );
    signalEnergy(2, jj-24) = structMis.energy{jj};
end


%% Plot the results
figure(1)
subplot 121
stem(angles,signalEnergy(1,:), 'filled', 'LineStyle','-.',...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green')
 title('sine sweep')
 xlabel('angle [°]')
 ylabel('Energy Ex')
 
 subplot 122
stem(angles,signalEnergy(2,:), 'filled', 'LineStyle','-.',...
     'MarkerFaceColor','blue',...
     'MarkerEdgeColor','green')
 title('noise')
 xlabel('angle [°]')
 ylabel('Energy Ex')
