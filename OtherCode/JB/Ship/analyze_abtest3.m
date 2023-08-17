% June 7, test 3
% open and plot
% add lines for mode arrivals
% fourier filter
% filter data
% thumbprint filtered data
% make where circular feature is centered


names = ['a-b_30db_1=6_11';
         'a-c_50db_1=6_11';
         'a-d_50db_1=6_11';
         'a-e_70db_1=6_11';
         'a-f_70db_1=6_11';
         'a-g_80db_1=6_11';
         'a-h_80db_1=6_11';
         'h-a_80db_1=6_11';
         'h-z_90db_800_11';
         'h-y_90db_800_11';
         'h-x_90db_1=6_11';
         'h-w_90db_800_28';
         'h-v_90db_800_30'];

 dist = [31-4;    %propagation distance
         58-4;
         85-4;
         109-4;
         112-4+4;
         112-4+19;
         112-4+24;
         24+112-4;
         24+112+13+4;
         24+112+13+31;
         24+112+13+58;
         24+112+13+85;
         24+112+13+109]*25.4; % in mm
     
 mag = [10;%magnitude of the windowed filtered signal
        20;
        20;
        50;
        50;
        50;
        60;
        150;
        100;
        150;
        100;
        120;
        150];

 arr = [745;%arrival of S1 in points taken from previouse runs
        937;
        966;
        1088;
        853;
        1274;
        824;
        296;%800;
        410;%323;
        358;%853;
        361;%1191;
        482;%1208;
        612];%1161];
    
fig = [1,1,7;%helps separate out different distances and figures
       1,2,7;
       1,3,7;
       1,4,7;
       1,5,7;
       1,6,7;
       1,7,7;
       2,1,7;%helps separate out different distances and figures
       2,2,7;
       2,3,7;
       2,4,7;
       2,5,7;
       2,6,7];
        
 v_s1 = 4.8;
 v_a0 = 3.1;
 v_a1 = 2.85;
 v_s0 = 2.65;
 sr = 25;
 rs = 5;

samplingtime=.002; 
% Total duration of sampled signal, in sec, millisec, or microsec.
samplerate=25000000;  
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
% mode=1 for band-reject (notch) filter

 t = 1;
%  for i = [1,3,7,12]
for i = 8:13
     S1 = ((dist(i)/v_s1)*sr)/rs;
     A0 = ((dist(i)/v_a0)*sr)/rs;
     A1 = ((dist(i)/v_a1)*sr)/rs;
     S0 = ((dist(i)/v_s0)*sr)/rs;
     signal = load(names(i,:));
     [a,b] = size(signal)
     wave = zeros(1,b);
     for k = 1:a
         wave = wave + signal(k,:);
     end
     wave = wave/10;
     wave = wave - mean(wave);     
%      wave = wave(1:50000);
     
        % Fourier filter 
        fy=fft(wave);
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

        rawdata = resample(wave,1,rs);
        fourierdata = resample(ry,1,rs);
%         rawdata = rawdata(1:6000);
%         fourierdata = fourierdata(1:6000);
        datatofilt = rawdata;
             
        % -- filt data
            wvtpf = 'coif3';% Wavelet Name
            swatoremove = [];% aproxamation levels to remove
            swdtoremove = [1:4];% detail levels to remove
            numlvls = 5;   % Number of Levels to Use

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
            wtpleft = floor(S1-500);
            if wtpleft <= 0
                wtpleft = 1;
            end
            wtpright = wtpwidth+wtpleft-1;
            wvttp = 'mexh';% Wavelet Name
            ns = 50;       % Number of Levels to Use
            nr = 10;       % Number of Ridges
            rw = .05;     % Ridge Width
            gv = .6;       % Valley multiplier
            op = 0;  % Show: both(0),valleys(1),peaks(2)

        if (length(filtdata) >= wtpright)
            datatothumbprint = [filtdata(wtpleft:wtpright)];
        else
            datatothumbprint = [filtdata(...
                wtpleft:length(filtdata))];
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
    x = ([1:length(rawdata)]*rs)/sr;


    figure(18+fig(i,1)-1),subplot(7,1,fig(i,2))
    imshow(thumbprint)
    xLabel(['S1 arrival ' num2str(...
        ((arr(i))+wtpleft)*rs/sr) ' us'])
    hold on
    line(ones(200,1)*(arr(i)),(1:200),'Color',...
        'r','LineWidth',2,'LineStyle','-')
    hold off
 
    pts(t) = (arr(i)+wtpleft)*rs/sr;
    expt(t) = S1;

    t = t+1;

        
     figure(17+fig(i,1)-1),subplot(6,1,fig(i,2))
     plot(x(1:length(rawdata)),rawdata,'k')
     hold on
     line(ones(4001,1)*S1,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','--')
     line(ones(4001,1)*A0,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','--')
     line(ones(4001,1)*A1,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','--')
     line(ones(4001,1)*S0,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','--')
     line(ones(4001,1)*wtpleft*rs/sr,(-2000:2000),...
     'Color','k','LineWidth',2)
     line(ones(4001,1)*wtpright*rs/sr,(-2000:2000),...
         'Color','k','LineWidth',2)
     hold off
     figure(17),hold on
     plot(dist(i),S1,'ko','MarkerSize',10,'MarkerFaceColor','k')
     plot(dist(i),(arr(i)+wtpleft)*rs/sr,'ro','MarkerSize',...
         10,'MarkerFaceColor','r')
     hold off
 end
