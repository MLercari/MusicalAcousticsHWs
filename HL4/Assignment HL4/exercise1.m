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

for i=2:length(angles)
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
        cont = dir(dirRec(j)); %struct with content of the directory
        cont = cont(4:end); %remove non useful components
        
        %sort the struct directory in ascendent order
        T = struct2table(cont); % convert the struct array to a table
        
        for k = 1:length(cont)
          T.name{k} = erase(T.name{k},'.wav'); %erase the last string part .wav  
        end
        
        
        %T.name = str2double(T.name); %convert the remaining into a number
        sortedT = sortrows(T, 'name'); %sort the name field of the table
        sortedCont = table2struct(sortedT); % change it back to struct array
   
    elseif(dirRec(j) == "Recordings/noise")
        cont = dir(dirRec(j)); %content of the directory
        cont = cont(4:end); %remove non useful components
        
        %sort the struct directory in ascendent order
        T = struct2table(cont); % convert the struct array to a table
        
        for k = 1:length(cont)
          T.name{k} = erase(T.name{k},'.wav'); %erase the last string part .wav  
        end
        
        %T.name = str2double(T.name); %convert the remaining into a number
        sortedT = sortrows(T, 'name'); %sort the name field of the table
        sortedCont = table2struct(sortedT); % change it back to struct array   

    end
    
end



% Labels will be the angle 

%% Load the signals and compute the energy

%         for k=1:length(cont) %load each 
%             subDir = strcat(dirRec(j), '/', k, '.wav'); 
%         end
%% Plot the results
