% EFIT_cart_3D_inputfiles.m
% Generates configuration files for EFIT simulations
% Edited version of EFIT_3D_cart_lamb, to work with updated 'EFIT_cart.cpp'
% Input files will not work with previous 'EFIT.cpp' simulations - uses (x,y,z) instead of (y,x,z)
%
% Usage: EFIT_cart_3D_inputfiles(~)
% In: NULL
% Out: NULL
% Dependencies:

% Eric A. Dieckman (W&M)
% 06 September 2011
% Last edited: 08 Sept 2011 EAD


function [ output_args ] = EFIT_cart_3D_inputfiles(~)

%% Material Parameters
% Aluminum
den = 2698; % kg/m^3
cL = 6235; % m/s 6374
cT = 3139; % m/s 2906
mu = den*cT^2; % Pa
lambda = den*cL^2-2*mu; % Pa

% Brass
% den = 8400; % kg/m^3
% cL = 4400; % m/s
% cT = 2200; % m/s
% mu = den*cT^2; % Pa
% lambda = den*cL^2-2*mu; % Pa

% Another way to calculate sound speeds:
% cT=sqrt(mu/den) % results in 3139 m/s for 2024 Al
% cL=sqrt((lambda+2*mu)/den) % results in 6235 m/s for 2024 Al

cmax = cL;
cmin = cT;
clear cL cT

%% Simulation Parameters
fmax = 1.00*10^6; % max frequency (Hz)
wavelength=cmin/fmax;
ds=wavelength/12 % step size (m) -> here we assign number of points per wavelength (>6)

% Size of simulation space (in mm):
x_mm = 360; 
y_mm = 360; 
z_mm = 1.62; 

% Convert size from m to steps:
maxx = round(x_mm/(1000*ds))
maxy = round(y_mm/(1000*ds))
maxz = round(z_mm/(1000*ds))

maxt = 12001; % max simulation time in steps (+1)
outputevery = 100; % output 3D volume how often
dt = 1/(cmax*sqrt(3/(ds^2))) % time step (s)


%% Choose a flaw type (everything should be in number of steps)
%NOTE: y-intercept of line equations can be adjusted in spipe.h code
% Rectangular:
nS     = 1;                      % numref // number of scatterers
rftype = [3];                    % typ    // Reflector type: 3 - Right Rectangular Prism
nsx1   = [round(179.5/(1000*ds))];  % p1     // x-start (x-end will be to side of space)
nsx2   = [round(180.5/(1000*ds))]; % p2     // x-end
nsx3s  = [round(160/(1000*ds))];  % start3 // y-start
nsx3e  = [round(200/(1000*ds))];  % end3   // y-end
rrad   = [round(1.5/(1000*ds))];  % rad    // z-start
rden   = maxz -1; % [round(1.6/(1000*ds))];  % dd     // z-end
rmu    = [0];                    % mu     // null
rlambda= [0];                    % lambda // null

% Spherical:
% nS     = 1;               % numref // numbers scatterers
% rftype = [1];             % typ    // Reflector type: 1 - Sphere
% nsx1   = [round(maxx/2)]; % p1     // x-dim
% nsx2   = [round(maxy/2)]; % p2     // y-dim
% nsx3s  = [round(maxz/2)]; % start3 // z-dim
% nsx3e  = [0];             % end3   // null
% rrad   = [round(maxz/2)-2];% rad    // radius of sphere
% rden   = [-1];            % dd     // density - null (-1 = rigid)
% rmu    = [0];             % mu     // mu - null
% rlambda= [0];             % lambda // lambda - null

% Crack:
% nS=1                    % numref // numbers scatterers
% rftype=[ 1 ]            % typ    // Reflector type: for 2D crack type=0 (in single y-plane) , for  3d crack choose type=1
% nsx1  =[2]              % p1     // 1st slope
% nsx2  =[0.5]            % p2     // 2nd slope for 3d crack
% nsx3s =[round(maxz/2) ] % start3 // start z direction, currently end z is at top surface
% nsx3e =[ numy/2-20]     % end3   // start y direction
% rrad  =[numy/2]         % rad    // end y direction
% rden  =[numx/2-20]      % dd     // start x direction
% rmu    =[numx/2];       % mu     // end x direction
% rlambda= [0];           % lambda // not used


%% Transducer info (everything in number of steps)
ntrans = 1;
transducer_x = [round(72.5/(1000*ds))]; % transducer x position
transducer_y = [round(190/(1000*ds))]; % transducer y position
transducer_z = [maxz-1]; % transducer z position; (maxz - 1) is positioned on top
transducer_rad = [round(.0035/ds)]; % transducer radius

dfpulselen = 5*(1/fmax);  % Pulse length (seconds), 5 cycles
dffreq = fmax;  % Pulse Frequency (Hz)

df(1:maxt) = 0;
dfl = ceil(dfpulselen/dt);
df(1:dfl) = 10^6*sin((0:(dfl-1))*dt*dffreq*2*pi);
drivelen = length(df);

% plot(df)
%dt = 1/(sqrt(3/ds^2))


%% Write files
%==========================================================================================
% Write Inputfiles for simulation - % DO NOT CHANGE THE ORDER OF THIS PART
%==========================================================================================
[fname,pname] = uiputfile('in.file', 'Save Configuration');
fp=fopen('in.file','w');

fprintf(fp, ' %8.0f ' , maxx);      % simparams[0] - num1 (will be +2)
fprintf(fp, ' %8.0f ' , maxy);      % simparams[1] - num2
fprintf(fp, ' %8.0f ' , maxz);      % simparams[2] - num3
fprintf(fp, ' %2.20f ', ds);        % simparams[3] - ds
fprintf(fp, ' %2.20f ', dt);        % simparams[4] - dt

fprintf(fp, ' %15.6f ', den);       % simparams[5] - den
fprintf(fp, ' %15.6f ', lambda);    % simparams[6] - lm
fprintf(fp, ' %15.6f ' , mu);       % simparams[7] - mu
fprintf(fp, ' %8.0f ' , maxt);          % maxt
fprintf(fp, ' %8.0f ' , outputevery);   % outputevery

fprintf(fp, ' %8.0f ' , nS);            % numref
for i = 1:nS                                              % addReflector
    fprintf(fp, ' %8.0f ' ,  rftype(i));    % rpars[0]        typ
    fprintf(fp, ' %15.6f ' , nsx1(i));      % rpars[1]        p1
    fprintf(fp, ' %15.6f ' , nsx2(i));      % rpars[2]        p2
    fprintf(fp, ' %8.0f ' ,  nsx3s(i));     % rpars[3]        start3
    fprintf(fp, ' %8.0f ' ,  nsx3e(i));     % rpars[4]        end3
    fprintf(fp, ' %15.6f ' , rrad(i));      % rpars[5]        rad
    fprintf(fp, ' %15.6f ' , rden(i));      % rpars[6]        dd
    fprintf(fp, ' %15.6f ' , rmu(i));       % rpars[7]        mu
    fprintf(fp, ' %15.6f ' , rlambda(i));   % rpars[8]        lambda
end

%fprintf(fp, ' %s ', [ pname ]);        % working directory
fclose(fp);

   
[fname,pname] = uiputfile('trans.file', 'Save Configuration');
fp=fopen('trans.file','w');
fprintf(fp, ' %8.0f ' , ntrans);                    % numtrans
for i=1:ntrans
    fprintf(fp, ' %8.0f ' , transducer_x(i) );      % tparams[0] // tposx; // transducer x location 
    fprintf(fp, ' %8.0f ' , transducer_y(i) );      % tparams[1] // tposy; // transducer y location
    fprintf(fp, ' %2.20f ', maxz-2 );               % tparams[2] // tposz; // transducer z location  --> because always on top
    fprintf(fp, ' %2.20f ', transducer_rad(i));     % tparams[3] // trad;  // transducer radius
    fprintf(fp, ' %2.20f ', drivelen);           % tparams[4] // drivelen; // length of drive function
    fprintf(fp, ' %15.6f ', df(1:maxt));            % drive[i] where i = tparams[4]  % MUST BE ALTERED TO INCLUDE MULTIPLE Drive functions
end
   
fclose(fp);
