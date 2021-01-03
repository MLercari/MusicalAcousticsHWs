%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2020                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
fs = 4*44.1e3; %[Hz] sampling frequency
%fs = 4*44.1e3; %[Hz] suggested by Saitis  -> it takes a lot of time
ts = 1/fs;
duration = 8; %[s] signal length
tSamples = duration*fs;
t = linspace(0, duration, tSamples); %temporal vector

% Fundamental note
f1 = 52.8221; %[Hz] frequency of C2

% Boundary          

% String parameters
b1 = 0.5; %air damping coefficient
b2 = 6.25e-9; %internal friction of the string coefficient
T0 = 750;  % [N] tension at rest
L = 1.92;    % [m] length of the string
M_s = 35e-3;    % [kg] mass of the string
mu = M_s/L;   % [kg m^-1] avg linear mass density of the string
eps = 7.5e-6; % [-] string stiffness parameter assumed equal to kappa- string stiff. coefficient
c = sqrt(T0/mu);    % [m/s] string wave velocity
x0 = 0.12*L; %[m] stricking position (see d/L from Fletcher)

% Aliasing condition
if (fs < 2*f1)
    disp("WARNING: SAMPLING FREQUENCY BELOW NYQUIST");
end

% Spatial sampling parameters
X_max = sqrt( 0.5*((c*ts)^2 + 4*b2*ts + sqrt( ((c*ts)^2 + 4*b2*ts)^2 + 16*(eps*ts)^2 )) ); %stability condition (see Saitis)

% Number of maximum spatial steps
gamma = fs/(2*f1);

% Integer values
N_max = ceil(sqrt( (-1 + sqrt(1 + 16*eps*(gamma^2)))/(8*eps) ));
X_min = L/N_max;

% Spatial sampling
N = N_max - 1;
X = L/N;    %as close as possible to Xmin (see Chaigne et al.)
%X = X_max - 0.0005; %according to Saitis
%X = L(521; %suggested value by Saitis 

if ( X > X_max )
    disp("WARNING: WRONG SPATIAL SAMPLING STEP");
end


x = linspace( 1 , L , L/X  ); %spatial coordinate vector 

% FD parameters
lambda = c*ts/X; %Courant parameter
nu = 2*b2*ts/(ts^2); %convenient variable n. 1
mi = eps^2/(c^2*X^2); %convenient variable n.2

% Hammer parameters
Vh0 = 2.5;  % [m/s] initial hammer velocity
af = (ts^2/mu)/(1 + b1*ts); %force coefficient
bh = 1e-4; %[s^-1] fluid damping coefficient
Mh = 4.9e-3; %[Kg] hammer mass
eta = zeros(1 , length(t)); %hammer displacement vector
Fh = zeros(1 , length(t)); %hammer force vector

% Hammer contact window definition
m0 = ceil(x0/X); %number of sampling of hammer stricking coordinate
wh = 0.02; %[m] width of the hammer
w = ceil(wh/X); %sample of the hammer's width
gloc = hann(2*w);
g = zeros(1, length(x));
g = [g(1:(m0 - w)-1) , gloc' , g(((m0 - w)+length(gloc)):end)];
%g(m0-w:m0-w+length(gloc)-1)=gloc';

%PDE Coefficients:

% ai coefficients
a1 = (-lambda^2*mi)/(1 + b1*ts);
a2 = (lambda^2 + 4*lambda^2*mi + nu)/(1 + b1*ts);
a3 = (2 - 2*lambda^2 - 6*lambda^2*mi - 2*nu)/(1 + b1*ts);
a4 = (-1 + b1*ts + 2*nu)/(1 + b1*ts);
a5 = -nu/(1 + b1*ts);

%di coefficients
d1 = 2/(1 + bh*ts/(2*Mh));
d2 = (-1 + bh*ts/(2*Mh))/(1 + bh*ts/(2*Mh));
df = (-ts^2/Mh)/(1+bh*ts/(2*Mh));

% Bridge boundary coefficients
zeta_b = 1000; %normalized impedance at the bridge

%bri coefficients
br1 = (2 - 2*lambda^2*mi - 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br2 = (4*lambda^2*mi + 2*lambda^2)/(1 + b1*ts + zeta_b*lambda);
br3 = (-2*lambda^2*mi)/(1 + b1*ts + zeta_b*lambda);
br4 = (-1 -b1*ts + zeta_b*lambda)/(1 + b1*ts + zeta_b*lambda);
br5 = (ts^2/mu)/(1 + b1*ts + zeta_b*lambda);
brf = (ts^2/mu)/(1 + b1*ts + zeta_b*lambda);

% Left hand (hinged string end) boundary coefficients
zeta_l = 1e20; %normalized impedance

%bli coefficients
bl1 = (2 - 2*lambda^2*mi - 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl2 = (4*lambda^2*mi + 2*lambda^2)/(1 + b1*ts + zeta_l*lambda);
bl3 = (-2*lambda^2*mi)/(1 + b1*ts + zeta_l*lambda);
bl4 = (-1 -b1*ts + zeta_l*lambda)/(1 + b1*ts + zeta_l*lambda);
bl5 = (ts^2/mu)/(1 + b1*ts + zeta_l*lambda);
blf = (ts^2/mu)/(1 + b1*ts + zeta_l*lambda);

% Hammer felt parameters
p = 2.3; 
a = 0.12;
K = 4e8; 

%% Computation of the FD scheme
% Initialization
y = zeros(length(t), length(x));
F = zeros(length(t), length(x));

% Computation loop
for n = 1:length(t)-1
    
  %initial conditions 
  if(n ==1) %at time t= 0 
      %the string is at rest position
      y(n,:) = 0;
      %hammer is set to be in contact with the string
      eta(n) = 0;
      n = n + 1; %in order to prevent error in the space loop
      
  elseif(n ==2) %at time t = dt
      eta(n) = Vh0*ts;
      Fh(n) = K*(abs(eta(n) - y(n, m0)))^p;
      
  else
        % hammer displacement approximation
        eta(n+1) = d1*eta(n) + d2*eta(n-1) + df*Fh(n);
        
        %relative stricking force
        Fh(n) = K*(abs(eta(n) - y(n, m0)))^p; 
  

  end
    
    for m = 1:length(x)
    
      %forcing term definition at any time step
    if(eta(n) < y(n, m0)) %if the hammer is no more in contact with the string
        Fh(n) = 0;
        
    else
        F(n,m) = Fh(n)'*g(m); %windowed with the hanning
    end
    
    %boundary conditions
    if(m==1) %if m = 0
        y(n+1, m) = bl1*y(n, m) + bl2*y(n, m+1) + bl3*y(n, m +2) + ...
            bl4*y(n-1, m) + blf*F(n,m);
    
    elseif (m== length(x)) % if m = M
        y(n+1, m) = br1*y(n, m) + br2*y(n, m-1) + br3*y(n, m -2) + ...
            br4*y(n-1, m) + brf*F(n,m);
        
    elseif( m ==length(x)-1 ) % if m = M -1
        y(n+1,m) = a1*(2*y(n, m+1) - y(n,m) + y(n,m-2)) + a2*(y(n,m+1) + y(n,m-1)) + ...
            a3*y(n,m) + a4*y(n-1,m) + a5*(y(n-1,m+1) + y(n-1,m-1)) + af*F(n,m);
        
    elseif( m == 2) %if m = 1
        y(n+1,m) = a1*(y(n, m +2) - y(n,m) + 2*y(n,m-1)) + a2*(y(n, m+1) + y(n,m-1)) + ...
            a3*y(n,m) + a4*y(n-1,m) + a5*(y(n-1,m+1) + y(n-1,m-1)) + af*F(n,m);
        
    else 
        %backward finite difference scheme "inside the domain"
        y(n+1,m) = a1*(y(n, m+2) + y(n, m-2)) + a2*(y(n, m+1) + y(n, m-1)) + ...
            a3*y(n,m) + a4*y(n-1,m) + a5*(y(n-1, m+1) + y(n-1,m-1)) + af*F(n,m);
    
    end
    
    end
    
end

%% Plot the displacement in time
for i = 1:50:length(t)
    
    %50 step per iteration gives a nice result 
    time = i*ts;
    figure(1);
    plot(x, y(i,:));
    title(['y(x,t) at time : ', num2str(time), ' s']); xlabel('x-axis'); ylabel('y-axis');
    ylim([-6e-11,6e-11])
    
  
end

%% Plot the synthesized signal play it and save it on the disk

% Play the sound
 r0 = ceil((L-x0)/X); %coordinate of the striking specular point
 yspec = [y(:,r0+1)' ; y(:,r0+2)'; y(:,r0+3)'; y(:,r0+4)'; y(:,r0+5)'; y(:,r0+6)'; ...
     y(:,r0-1)'; y(:,r0-2)'; y(:,r0-3)'; y(:,r0-4)'; y(:,r0-5)'; y(:,r0-6)'];
 yaudio = mean(yspec);
 
 plot(t, yaudio); %Plot the estimated sound signal
 
 %sound(yaudio, fs); %Play the sound
 
 filename = 'strike.wav';
 audiowrite(filename, 1e10.*yaudio,fs);


% Save on disk

%EXTRA: plot the time domain of the force in the striking point









