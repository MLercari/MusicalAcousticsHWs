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

typeOfSignal = ["noise", "sweep"];       % Noise or sweep

dirRec = ["Recordings/sweep", "Recordings/noise" ];  % Recordings directory

angles = zeros(1,24);
for i=2:length(angles) %compute the vector that contains the angles (in degree)
    angles(i) = 15 + angles(i-1);
end

StructMis = struct(... %struct object that contain any measurement info
    'fileName', {}, ...
    'type', {}, ...
    'angle', {});

nMic = 1;                        % Number of microphones
signalEnergy = zeros(2, 24);     % Vector of energy

for j = 1:length(dirRec)
    
    if(dirRec(j) == "Recordings/sweep")
        Scont = dir(dirRec(j)); %struct with content of the directory
        Scont = Scont(4:end); %remove non useful components
        
        %sort the struct directory in ascendent order
        ST = struct2table(Scont); % convert the struct array to a table
        
        for k = 1:length(Scont)
          ST.name{k} = erase(ST.name{k},'.wav'); %erase the last string part .wav  
          ST.name{k} = str2double(ST.name{k}); %convert the remaining into a number
        end

        sortedST = sortrows(ST, 'name'); %sort the name field of the table
        sortedSCont = table2struct(sortedST); % change it back to struct array
   
    elseif(dirRec(j) == "Recordings/noise")
        Ncont = dir(dirRec(j)); %content of the directory
        Ncont = Ncont(4:end); %remove non useful components
        
        %sort the struct directory in ascendent order
        NT = struct2table(Ncont); % convert the struct array to a table
        
        for k = 1:length(Ncont)
          NT.name{k} = erase(NT.name{k},'.wav'); %erase the last string part .wav  
          NT.name{k} = str2double(NT.name{k}); %convert the remaining into a number
        end
        
        sortedNT = sortrows(NT, 'name'); %sort the name field of the table
        sortedConNt = table2struct(sortedNT); % change it back to struct array   

    end
    
end



% Labels will be the angle 

%% Load the signals and compute the energy

%         for k=1:length(cont) %load each 
%             subDir = strcat(dirRec(j), '/', k, '.wav'); 
%         end
%% Plot the results
