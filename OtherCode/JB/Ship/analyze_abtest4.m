% June 8, test 2
% open and plot
% add lines for mode arrivals
% fourier filter
% filter data
% thumbprint filtered data
% make where circular feature is centered


names = ['L1-L2_32_60db';
         'L1-L3_32_70db';
         'L4-L3_32_80db';
         'L5-L3_32_90db'];

 dist = [3967;    %propagation distance
         6435;
         8767;
         10556]; % in mm
     
 mag = [50;%magnitude of the windowed filtered signal
        30;
        50;
        100];

 arr = [783;%arrival of S1 in points taken from previouse runs
        763;
        1043;
        821];
    
fig = [1,1,4;%helps separate out different distances and figures
       1,2,4;
       1,3,4;
       1,4,4];
        
 v_s1 = 4.8;
 v_a0 = 3.1;
 v_a1 = 2.85;
 v_s0 = 2.65;
 sr = 25;
 rs = 5;

 eS1 = [846;
      1328;
      1758;
      2175]*sr/rs;
  
samplingtime= .004; 
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
for i = 1:4
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
            wtpleft = floor(eS1(i)-700)
            if wtpleft <= 0
                wtpleft = 1;
            end
            wtpright = wtpwidth+wtpleft-1;
            wvttp = 'mexh';% Wavelet Name
            ns = 50;     % Number of Levels to Use
            nr = 10;     % Number of Ridges
            rw = .05;     % Ridge Width
            gv = .6;     % Valley multiplier
            op = 0;    % Show: both(0),valleys(1),peaks(2)

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

   %     switch units to microseconds
    eS1(i) = (eS1(i)*rs)/sr;
    S1 = (S1*rs)/sr;
    A0 = (A0*rs)/sr;
    A1 = (A1*rs)/sr;
    S0 = (S0*rs)/sr;
    x = ([1:length(rawdata)]*rs)/sr;

     figure(10+fig(i,1)-1),subplot(4,1,fig(i,2))
     plot([1:length(rawdata)]*1/sr*rs,rawdata,'k')
     axis([0 4000 -2000 2000]);
     hold on
     line(ones(4001,1)*eS1(i),(-2000:2000),'Color',...
         'g','LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*wtpleft*rs/sr,(-2000:2000),...
         'Color','k','LineWidth',2)
     line(ones(4001,1)*wtpright*rs/sr,(-2000:2000),...
         'Color','k','LineWidth',2)
     hold off

    figure(11+fig(i,1)-1),subplot(4,1,fig(i,2))
    imshow(thumbprint)
    xLabel(['S1 arrival ' num2str(((arr(i))+wtpleft)*rs/sr) ' us'])
    hold on
    line(ones(200,1)*(arr(i)),(1:200),'Color','r',...
        'LineWidth',2,'LineStyle','-')
    hold off

    t = t+1;

        
%      figure(fig(i,1)),subplot(4,1,fig(i,2))
     figure(i),subplot(3,1,3)
     plot([1:length(rawdata)]*1/sr*rs,rawdata,'k')
%      plot((wtpleft:wtpright),rawdata(wtpleft:wtpright))
     axis([0 4000 -2000 2000]);
     hold on
     line(ones(4001,1)*S1,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*A0,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*A1,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*S0,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','-')
     hold off
     subplot(3,1,2)
     plot(x(1:length(filtdata)),filtdata,'k')
     axis([0 4000 -mag(i) mag(i)]);
%      axis([wtpleft*rs/sr wtpright*rs/sr -mag(i) mag(i)]);
     hold on
     line(ones(4001,1)*S1,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*A0,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*A1,(-2000:2000),'Color','b',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*S0,(-2000:2000),'Color','r',...
         'LineWidth',2,'LineStyle','-')
     line(ones(4001,1)*wtpleft*rs/sr,(-2000:2000),...
         'Color','k','LineWidth',2)
     line(ones(4001,1)*wtpright*rs/sr,(-2000:2000),...
         'Color','k','LineWidth',2)
     hold off
     subplot(3,1,1)
     imshow(thumbprint)
     hold on
     line(ones(200,1)*(arr(i)),(1:200),'Color','r',...
         'LineWidth',2,'LineStyle','-')
    hold off
 end
