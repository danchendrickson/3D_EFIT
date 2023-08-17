% June 6, test 1
% open and plot
% add lines for mode arrivals
% fourier filter
% filter data
% thumbprint filtered data
% make where circular feature is centered


names = ['nom1_40db_1w_0s';
         'mag1_40db_1w_0s';
         'nom2_60db_1w_1s';
         'mag1_60db_1w_1s';
         'mag2_60db_1w_1s';
         'mag3_60db_1w_1s';
         'nom3_70db_1w_2s';
         'mag1_70db_1w_2s';
         'mag2_70db_1w_2s';
         'mag3_70db_1w_2s';
         'mag4_70db_1w_2s';
         'nom4_70db_2w_3s';
         'mag1_70db_2w_3s';
         'mag2_70db_2w_3s';
         'mag3_70db_2w_3s';
         'mag4_70db_2w_3s'];

 dist = [5.5+10;
         5.5+10;
         5.5+31;
         5.5+31;
         5.5+31;
         5.5+31;
         5.5+59;
         5.5+59;
         5.5+59;
         5.5+59;
         5.5+59;
         5.5+83+3;
         5.5+83+3;
         5.5+83+3;
         5.5+83+3;
         5.5+83+3]*25.4; % in mm
     
 mag = [25;
        25;
        100';
        100;
        100;
        100;
        100;
        100;
        100;
        100;
        100;
        50;
        50;
        50;
        50;
        50];

 arr = [488;
        482;
        540';
        534;
        532;
        524;
        493;
        475;
        472;
        472;
        482;
        527;
        437;
        450;
        428;
        436];
        
 v_s1 = 4.8;
 v_a0 = 3.1;
 v_a1 = 2.85;
 v_s0 = 2.65;
 sr = 25;
 rs = 5;

samplingtime=.002; 
% Total duration of sampled signal, in sec, millisec,
% or microsec.
samplerate=25000000;  
% Sample rate in Hz, KHz, or MHz, respectively.
n=samplingtime*samplerate; % Number of points in signal
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
% mode=1 for band-reject (notch)
% filter

 t = 1;
 for i = [1,3,7,12]
% for i = 1:16
     S1 = ((dist(i)/v_s1)*sr)/rs;
     A0 = ((dist(i)/v_a0)*sr)/rs;
     A1 = ((dist(i)/v_a1)*sr)/rs;
     S0 = ((dist(i)/v_s0)*sr)/rs;
     signal = load(['D:\Jill\NDE\Projects\REMORA\Aberdeen\'...
         'June 6 2006\test1\' names(i,:)]);
     wave = zeros(1,50048);
     for k = 1:10
         wave = wave + signal(k,:);
     end
     wave = wave/10;
     wave = wave - mean(wave);     
     wave = wave(1:50000);
     
        % Fourier filter 
        fy=fft(wave); % Compute Fourier transform of signal wave
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
        if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];
        end
        ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier 
%         transform of signal
        ry=real(ifft(ffy)); % Inverse transform to recover filtered 
%         signal 'ry'

        rawdata = resample(wave,1,rs);
        fourierdata = resample(ry,1,rs);
        rawdata = rawdata(1:6000);
        fourierdata = fourierdata(1:6000);
        datatofilt = rawdata;
             
        % -- filt data
            wvtpf = 'coif3';        % Wavelet Name
            swatoremove = [];       % aproxamation levels to remove
            swdtoremove = [1:4];      % detail levels to remove
            numlvls = 5;            % Number of Levels to Use

        datatofilt = datatofilt(1:(length(datatofilt)-rem(length ...
            (datatofilt),2^numlvls)));        
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
            wtpleft = floor(S1-500);
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

        datatothumbprint = [abs(filtdata(wtpleft:wtpright))];
        datatothumbprint = datatothumbprint.*tukeywin(length ... 
            (datatothumbprint),.25);%';

        thumbprintpeaks   = getThumbprint( ...
            datatothumbprint, wvttp, ns, (1), nr, rw, 2 );
        % get thumbprint for peaks
        thumbprintvalleys = getThumbprint( ...
            datatothumbprint, wvttp, ns, (1), nr, rw, 3 );
        % get thumbprint for valleys

        if  op == 1                             % Show only Valleys
            thumbprint = thumbprintvalleys.*gv;                             
        elseif  op == 2                         % Show only Peaks
            thumbprint = thumbprintpeaks;
        else                                    % Show Both
            thumbprint = thumbprintpeaks + thumbprintvalleys.*gv;
        end
        
        %search for begining of feature
        arrival_s1 = 0;
        for k = 1:length(thumbprint)
            for j = ns:-1:130
                if (thumbprint(j,k) ~= 0)
                    arrival_s1 = k;
                    break;
                end
            end
            if (arrival_s1 > 0)
                break;
            end
        end
        
%     switch units to microseconds
    S1 = (S1*rs)/sr;
    A0 = (A0*rs)/sr;
    A1 = (A1*rs)/sr;
    S0 = (S0*rs)/sr;
    x = ([1:6000]*rs)/sr;

    figure(1),subplot(4,1,t)
    plot(x(1:length(rawdata)),rawdata,'k')
    hold on
    line(ones(4001,1)*S1,(-2000:2000),'Color','r','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*A0,(-2000:2000),'Color','b','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*A1,(-2000:2000),'Color','b','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*S0,(-2000:2000),'Color','r','LineWidth',2, ...
        'LineStyle','--')
%     line(ones(4001,1)*wtpleft*rs/sr,(-2000:2000),'Color','k', ...
% 'LineWidth',2)
%     line(ones(4001,1)*wtpright*rs/sr,(-2000:2000),'Color','k', ...
% 'LineWidth',2)
    hold off
    figure(2),subplot(4,1,t)
    plot(x(1:length(fourierdata)),fourierdata,'k')
    hold on
    line(ones(4001,1)*S1,(-2000:2000),'Color','r','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*A0,(-2000:2000),'Color','b','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*A1,(-2000:2000),'Color','b','LineWidth',2, ...
        'LineStyle','--')
    line(ones(4001,1)*S0,(-2000:2000),'Color','r','LineWidth',2, ...
        'LineStyle','--')
    hold off
    figure(3),subplot(4,1,t)
    plot(x(1:length(filtdata)),filtdata,'k')
    hold on
    line(ones(mag(i)*2+1,1)*S1,(-mag(i):mag(i)),'Color','r', ...
        'LineWidth',2,'LineStyle','--')
    line(ones(mag(i)*2+1,1)*A0,(-mag(i):mag(i)),'Color','b',...
        'LineWidth',2,'LineStyle','--')
    line(ones(mag(i)*2+1,1)*A1,(-mag(i):mag(i)),'Color','b',...
        'LineWidth',2,'LineStyle','--')
    line(ones(mag(i)*2+1,1)*S0,(-mag(i):mag(i)),'Color','r',...
        'LineWidth',2,'LineStyle','--')
    line(ones(mag(i)*2+1,1)*wtpleft*rs/sr,(-mag(i):mag(i)),...
        'Color','k','LineWidth',2)
    line(ones(mag(i)*2+1,1)*wtpright*rs/sr,(-mag(i):mag(i)),...
        'Color','k','LineWidth',2)
    hold off
    figure(4),subplot(4,1,t)
    plot((wtpleft*rs/sr:1/sr*rs:wtpright*rs/sr),filtdata...
        (wtpleft:wtpright),'k')
    axis([wtpleft*rs/sr wtpright*rs/sr -mag(i) mag(i)]);
    
    figure(5),subplot(4,1,t)
    imshow(thumbprint)
    hold on
    line(ones(200,1)*(arr(i)),(1:200),'Color','r',...
        'LineWidth',2,'LineStyle','-')
    hold off
    
    pts(t) = (arr(i)+wtpleft)*rs/sr;
    expt(t) = S1;
    
    t = t+1;

        
 end
figure(6)
plot([1:4],expt,'ko','MarkerSize',10,'MarkerFaceColor','k')
hold on
plot([1:4],pts,'ro','MarkerSize',10,'MarkerFaceColor','r')
    