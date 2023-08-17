% script to make workin with Day 1 T1 waveforms easier.


% Look at 1m separation direct waves
num = [1 3 4 5 7 6];
pts = [20240;
       20240;
       20240;
       20240;
       20240;
       20240];
mag1 = [200;
        200;
        200;
        400;
        400;
        400];
mag2 = [5;
        5;
        5;
        20;
        20;
        20];
label = [' Transverse 1 m       ';
         '    Italian           ';
         ' Large Russian        ';
         'Longitudinal 1 m      ';
         '    Italian           ';
         ' Large Russian        '];
    

samplingtime= .001; 
% Total duration of sampled signal, in sec, millisec, or microsec.
samplerate=20000000;  
% Sample rate in Hz, KHz, or MHz, respectively.
n=samplingtime*samplerate; 
% Number of points in signal
% % frequency=20;
x=[0:(1/samplerate):samplingtime];
x = x(1:n);

% Initial values of filter parameters
centerfrequency=500000;%364700;  
% Center frequency of filter, in Hz.
center=centerfrequency*samplingtime; 
%  center harmonic (fourier component)
frequencywidth=500000;  
% Frequency width of filter, in Hz.
width=frequencywidth*samplingtime;
%  width of filter (in harmonics)
shape=20; 
% filter shape (higher numbers = sharper cutoff)
mode=0;  
% mode=0 for band-pass filter, 
mode=1 for band-reject (notch) filter
   
   
% opens the files, only need to run once per session
for i = 1:6
    wave = openfile_be(['070830_4_1_' num2str(num(i))],...
        pts(i),100,'first_try');

    dist = ones(7,1)*1000;

    v_s2 = 4.7;
    v_s3 = 3.65;
    v_a3 = 3.1;
    v_s1 = 2.549;
    v_a4 = 2.55;
    % v_s4 = 1.141;

    sr = 20;
    rs = 2;

    S2 = (dist(i)/v_s2);
    S3 = (dist(i)/v_s3);
    A3 = (dist(i)/v_a3);
    S1 = (dist(i)/v_s1);
    A4 = (dist(i)/v_a4);
%     S4 = (dist(i)/v_s4);
    
    rawdata = wave-mean(wave);
    
    % Fourier filter 
    fy=fft(rawdata');
    % Compute Fourier transform of signal wave
    % Compute filter shape
    lft1=[1:(length(fy)/2)+1];
    lft2=[(length(fy)/2):length(fy)];
    ffilter1=ngaussian(lft1,center+1,width,shape);
    ffilter2=ngaussian(lft2,length(fy)-center+1,width,shape);
    ffilter=[ffilter1,ffilter2];
    modestring='Band-pass mode:';
    if mode==1, 
        ffilter=1-ffilter; 
        modestring='Band-reject (notch) mode:';
    end
    if length(fy)>length(ffilter)
        ffilter=[ffilter ffilter(1)];
    end
    ffy=fy.*ffilter(1:length(fy));  
    % Multiply filter by Fourier transform of signal
    ry=real(ifft(ffy)); 
    % Inverse transform to recover filtered signal 'ry'

    rawdata = resample(rawdata,1,rs);
    ffdata = resample(ry,1,rs);
    datatofilt = ffdata;%rawdata;

    
    % -- filt data
    wvtpf = 'coif3';        % Wavelet Name
    swatoremove = [];       % aproxamation levels to remove
    swdtoremove = [1:4];      % detail levels to remove
    numlvls = 5;            % Number of Levels to Use

    datatofilt = datatofilt(1:(length(datatofilt)...
        -rem(length(datatofilt),2^numlvls)));      
    % clip raw data to appropriate size for wavelet transform
    [swa,swd] = swt(datatofilt, numlvls, wvtpf); 
    % stationary wavelet transform
    swa(swatoremove,:)=0;                   
    % remove some approxamations 
    swd(swdtoremove,:)=0;               
    % remove some details 
    filtdata = iswt(swa, swd, wvtpf)';     
    % inverse stationary wavelet transform 

    % -- thumbprint
    wtpwidth = 1500;    % Window width
    wtpleft = floor(S2*sr/rs-700);
    if wtpleft <= 0
        wtpleft = 1;
    end
    wtpright = wtpwidth+wtpleft-1;
    wvttp = 'mexh';     % Wavelet Name
    ns = 50;            % Number of Levels to Use
    nr = 10;             % Number of Ridges
    rw = .05;            % Ridge Width
    gv = .6;            % Valley multiplier
    op = 0;             % Show: both(0),valleys(1),peaks(2)

    if (length(filtdata) >= wtpright)
        datatothumbprint = [filtdata(wtpleft:wtpright)];
    else
        datatothumbprint = [filtdata(wtpleft:length(filtdata))];
    end
    datatothumbprint = datatothumbprint.*...
        tukeywin(length(datatothumbprint),.25);%';

    thumbprintpeaks   = getThumbprint( ...
        abs(datatothumbprint), wvttp, ns, (1), nr, rw, 2 ); 
    % get thumbprint for peaks
    thumbprintvalleys = getThumbprint( ...
        abs(datatothumbprint), wvttp, ns, (1), nr, rw, 3 );
    % get thumbprint for valleys

    if  op == 1                             
        % Show only Valleys
        thumbprint = thumbprintvalleys.*gv;                             
    elseif  op == 2                         
        % Show only Peaks
        thumbprint = thumbprintpeaks;
    else
        % Show Both
        thumbprint = thumbprintpeaks + thumbprintvalleys.*gv;
    end


%     plot all rawdata in figure 1
    figure(1), subplot(6,1,i), hold off
    plot([1:length(rawdata)]/sr*rs,rawdata,'g')
    axis([0 500 -mag1(i) mag1(i)])
    hold on
    plot([1:length(ffdata)]/sr*rs,ffdata,'k')
%     ylabel(label(i,:),'Rotation',0.0)
    line(ones(801,1)*S2,(-400:400),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(801,1)*S3,(-400:400),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(801,1)*A3,(-400:400),'Color','b',...
        'LineWidth',2,'LineStyle','--')
    line(ones(801,1)*S1,(-400:400),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(801,1)*A4,(-400:400),'Color','b',...
        'LineWidth',2,'LineStyle','--')
%     line(ones(201,1)*S4,(-100:100),'Color','r',...
% 'LineWidth',2,'LineStyle','--')

%     plot all filtdata with window marks
    figure(2), subplot(6,1,i), hold off
    plot([1:length(filtdata)]/sr*rs,filtdata,'k')
%     ylabel(label(i,:),'Rotation',0.0)
    axis([0 500 -mag2(i) mag2(i)])
    hold on
    line(ones(101,1)*wtpleft*rs/sr,(-50:50),...
        'Color','k','LineWidth',2)
    line(ones(101,1)*wtpright*rs/sr,(-50:50),...
        'Color','k','LineWidth',2)
    line(ones(101,1)*S2,(-50:50),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(101,1)*S3,(-50:50),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(101,1)*A3,(-50:50),'Color','b',...
        'LineWidth',2,'LineStyle','--')
    line(ones(101,1)*S1,(-50:50),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(101,1)*A4,(-50:50),'Color','b',...
        'LineWidth',2,'LineStyle','--')
%     line(ones(21,1)*S4,(-10:10),'Color','r',...
% 'LineWidth',2,'LineStyle','--')
    
%     show all thumbprints together
    figure(3), subplot(6,1,i), hold off
    imshow(thumbprint)
%     ylabel(label(i,:),'Rotation',0.0)
    hold on


end
