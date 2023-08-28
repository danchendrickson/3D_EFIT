function [ output_args ] = run3DEFITsa_rectangle( input_args )

% Material and Pipe Parameters
% % =====================================================================

%******************************
%2024 Aluminum
den=2780;  % kg/m^3
cL = 6235; % m/s
cT = 3139; % m/s
mu = den*cT^2; % Pa
lambda = den*cL^2-2*mu; %Pa  ---also, note that cmin = 3139;   %m/s, longitudinal
% cT=sqrt(mu/den)  %results in 3139 m/s for 2024 Al
% cL=sqrt((lambda+2*mu)/den) %results in  6235 m/s for 2024 Al

% ----- These need to be changed with possible mode velocities!!!
%       Otherwise, simulation space step sizes are not large/small 
%       enough to accomidate the Lamb wave modes we generate!
cmax=cL;
cmin=cT;
% ----------------------------------------------------------------

% Above are bulk material properties
%***********************
% BRASS:
% den2=8400;   % kg/m^3
% cL = 4400;  % m/s
% cT = 2200;  % m/s
% mu2 = den*cT^2;   % Pa
% lambda2 = den*cL^2-2*mu; % Pa
%***********************
clear cL cT
% ======================================================================
% Simulation Parameters
% ======================================================================

%----Box---------
fmax = 1.25*10^6;    %Hz
%----Plate-------
% fmax = 2.15*10^6;  %Hz
%----------------
wavelength=cmin/fmax;
ds=wavelength/15; %step size in meters per step -- ds units is meters per step --
% here we assign certain number of points per wavelength

% Size of simulation space
%----Box---------
% x_mm = 100;
% y_mm = 100;
% z_mm = 100;
%----Plate-------
x_mm = 500;
y_mm = 300;
z_mm = 2.0;
%----------------

%Al plate - size must be number of steps
%space size: 305 mm (12 in) X 305 mm (12 in) X 3.154 mm
maxx = round(x_mm/1000/ds); %turned into # steps
maxy = round(y_mm/1000/ds);
maxz = round(z_mm/1000/ds);

maxt = 5001;        % max simulation time in number of time steps (+1)
outputevery = 100;    % output 3D volume every 2 time steps
                    % dt is in seconds per step
dt = 1/(cmax*sqrt(3/(ds^2))); % time step (sec)

% ---------------------Rectangular Reflector-------------------------------
%BELOW: converting flaw size in meters to #steps,
%so nsx1, etc. is in #steps
nS     = 0;                      % numref // numbers scatterers
rftype = [3];                    % typ    // Reflector type: 3 - Right Rectangular Prism
nsx1   = [round(50  /1000/ds)];  % p1     // x-start (x-end will be to side of space)
nsx2   = [round(100 /1000/ds)];  % p2     // x-end
nsx3s  = [round(25  /1000/ds)];  % start3 // y-start
nsx3e  = [round(75  /1000/ds)];  % end3   // y-end
rrad   = [round(40  /1000/ds)];  % rad    // z-start
rden   = [round(60  /1000/ds)];  % dd     // z-end
rmu    = [0];                    % mu     // null
rlambda= [0];                    % lambda // null
% -------------------------------------------------------------------------

% ---------------------Sphereical Reflector--------------------------------
% %BELOW: converting flaw size in meters to #steps,
% %so nsx1, etc. is in #steps
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
% -------------------------------------------------------------------------

% ------------------------Crack Reflector----------------------------------
% %BELOW: converting flaw size in meters to #steps,
% %so nsx1, etc. is in #steps
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
% -------------------------------------------------------------------------

%NOTE: y-intercept of line equations can be adjusted in spipe.h code
%itself


%==========================================================================================
% Write Inputfile for 3Dcefit simulation - do not change the order of this!
%==========================================================================================

[fname,pname] = uiputfile('in.file', 'Save Configuration');
fp=fopen('in.file','w');
%DO NOT CHANGE THE ORDER OF THIS PART

fprintf(fp, ' %8.0f ' , maxy);      % simparams[0] - num2
fprintf(fp, ' %8.0f ' , maxx);      % simparams[1] - num1 (will be +2)
fprintf(fp, ' %8.0f ' , maxz);      % simparams[2] - num3

fprintf(fp, ' %2.20f ', ds);        % simparams[3] - ds
fprintf(fp, ' %2.20f ', ds);        % simparams[4] - dp
fprintf(fp, ' %2.20f ', dt);        % simparams[5] - dt
fprintf(fp, ' %15.6f ', den);       % simparams[6] - den
fprintf(fp, ' %15.6f ', lambda);    % simparams[7] - lm
fprintf(fp, ' %15.6f ' , mu);       % simparams[8] - mu
fprintf(fp, ' %8.0f ' , 0);         % simparams[10] - rbeg
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


% Pitch Transducer info and Drive Function
% =====================================================================
ntrans=1;
%transmitter and receiver were stepped in 2mm increments
xpost=150;   %in mm
ypost=150;   %in mm
%***********************************
%EXPT
transducer_x   = round((xpost/1000)/ds);       % transducer x posintion (ds units)
transducer_y   = round((ypost/1000)/ds);       % transducer y posintion (ds units)
%****************************************
transducer_z   = maxz-1;   % transducer z position (ds units), positioned on top

transducer_rad = round(.0035/ds);   % transducer radius (ds units)  

dfpulselen = 5*(1/fmax);  % Pulse length (seconds), 5 cycles
dffreq     = fmax;  % Pulse Frequency (Hz)

df(1:maxt) = 0;
dfl = ceil(dfpulselen/dt);
df(1:dfl) = 10^6*sin((0:(dfl-1))*dt*dffreq*2*pi);

% plot(df)

drivelen=length(df);

%dt = 1/(sqrt(3/ds^2))

%==========================================================================================
% Write Inputfile for 3Dcefit simulation - do not change the order of this!
%==========================================================================================
   
[fname,pname] = uiputfile('trans.file', 'Save Configuration');
fp=fopen('trans.file','w');
fprintf(fp, ' %8.0f ' , ntrans);                    % numtrans
for i=1:ntrans
    fprintf(fp, ' %8.0f ' , transducer_x(i) );      % tparams[0] // tposx; // transducer x location 
    fprintf(fp, ' %8.0f ' , transducer_y(i) );      % tparams[1] // tposy; // transducer y location
    fprintf(fp, ' %2.20f ', maxz-2 );               % tparams[2] // tposz; // transducer z location  --> because always on top
    fprintf(fp, ' %2.20f ', transducer_rad(i));     % tparams[3] // trad;  // transducer radius
    fprintf(fp, ' %2.20f ', drivelen(i));           % tparams[4] // drivelen; // length of drive function
    fprintf(fp, ' %15.6f ', df(1:maxt));            % drive[i] where i = tparams[4]  % MUST BE ALTERED TO INCLUDE MULTIPLE Drive functions
end
   
fclose(fp);

%==========================================================================================
% Run the KZK simulation
%==========================================================================================
%dos(['java -Xmx1000m afit3D_sa_m user.cfg']);