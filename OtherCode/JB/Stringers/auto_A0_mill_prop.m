pathname = 'D:\Jill\NDE\Projects\Oceana\Corey Tests\Milling\averaged\';

dist = .8;
delay = 0;
numStep = 10;
% numStep = 14;
rs = 4;

steps = [1.5982,1.539,1.4544,1.2402,1.0048,0.786,...
    0.7546,0.683,0.6206,0.5382,0.488];

for i = 1:numStep+1
    if i > 10 
        filename{i} = ['Step_',num2str(i-1),'\T80NP3OS',num2str(i-1)];
    else
        filename{i} = ['Step_0',num2str(i-1),'\T80NP3OS0',num2str(i-1)];
    end
end
    % -- Compute expected mode arrivals
sr = 80;
v_s0 = 5.057;
v_a0 = 3.186;
SO = ((10^3*(dist/v_s0)*sr)-delay*sr)/rs;
AO = ((10^3*(dist/v_a0)*sr)-delay*sr)/rs;

num = length(filename{i});
first = zeros(1,num);
for i = 1:numStep+1    
    % -- Import Data
        sig_length = 60000;
        name = [pathname,filename{i}]
        fid = fopen(name);
    rawdata = fread(fid,sig_length,'int16');
    fclose(fid);
    rawdata = rawdata(1:40000);
    rawdata = rawdata-mean(rawdata);
    rawdata = resample(rawdata,1,rs);
    % -- filt data
        wvtpf = 'coif3';        % Wavelet Name
        swatoremove = [];       % aproxamation levels to remove
        swdtoremove = [1:4];      % detail levels to remove
        numlvls = 3;            % Number of Levels to Use

    rawdata = rawdata(1:(length(rawdata)-rem(...
        length(rawdata),2^numlvls)));       
    % clip raw data to appropriate size for wavelet transform
    [swa,swd] = swt(rawdata, numlvls, wvtpf); 
    % stationary wavelet transform
    swa(swatoremove,:)=0;                 
    % remove some approxamations 
    swd(swdtoremove,:)=0;              
    % remove some details 
    filtdata = iswt(swa, swd, wvtpf)';     
    % inverse stationary wavelet transform 

%         values = gettimes2(filtdata,1,80,i,numStep);
        wavelet_number = 12;        % 2^wavelet_number = 8192
        thres = 50;
        single_length = length(filtdata);
        gstart = 00;
        begin = 80;

        values.wave = filtdata;
        ppp = 2*2^wavelet_number
        filtdata = filtdata(begin:(begin+ppp-1));
        x1 = [begin:(begin+ppp-1)];
        [swa,swd] = swt(abs(filtdata),wavelet_number,'coif3');
        swd(1:(wavelet_number-3),:)=0;
        env1 = iswt(swa,swd,'coif3');
        values.line = [x1(:) env1(:)];
        clear x1 xx1 peaks1 valley1 env1

    % -- thumbprint
%             wtpwidth = 1500;    % Window width
%             wtpleft = 1;%floor(AO-100);
%             wtpright = wtpwidth+wtpleft-1;
        wvttp = 'haar';     % Wavelet Name
        ns = 200;            % Number of Levels to Use
        nr = 15;             % Number of Ridges
        rw = .05;            % Ridge Width
        gv = .6;            % Valley multiplier
        op = 0;             % Show: both(0),valleys(1),peaks(2)

    datatothumbprint = [zeros(2000,1); values.line(2001:8192,2)];
    datatothumbprint = datatothumbprint.*tukeywin(...
        length(datatothumbprint),.25);%';

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

%         Search thumbprint for double feature
    arrival1 = 0;
    arrival2 = 0;
    for k = 2100:length(thumbprint)
        v = 0;
        p = 0;
        tp = 0;
        tv = 0;
        for j = 1:ns
            if ((thumbprint(j,k)==gv) && (thumbprint(j,k) ~= p))
                tp = tp+1;
                p = gv;
            elseif (thumbprint(j,k) == 0)
                p = 0;
            end
            if ((thumbprint(j,k)==gv) && (thumbprint(j,k) ~= v))
                tv = tv+1;
                v = gv;
            elseif (thumbprint(j,k) == 0)
                v = 0;
            end
        end
        if ((tp >= 2) && (arrival1 == 0))
            arrival1 = k;
        end
        if ((tv >= 5) && (arrival2 == 0))
            arrival2 = k;
        end
    end 



%         figure out arrival in micro seconds
    arr1(i) = (arrival1)/(sr/rs);
    arr2(i) = (arrival2)/(sr/rs);

    if (i > 6)
        figure(4),hold on,subplot(6,1,i-6)
    else
        figure(1),hold on ,subplot(6,1,i)
    end
    plot([1:10000]/(sr/rs),rawdata)
    axis([0 500 -2000 2000])
    if (i > 6)
        figure(5),hold on,subplot(6,1,i-6)
    else
        figure(2),hold on,subplot(6,1,i)
    end
    imshow(thumbprint)
    hold on
    title(['Thumb prints Incremental Milling ' num2str(steps(i)) 'mm'])
    YLabel([num2str(arr2(i)) 'us'],'Rotation',0)

    figure(3),hold on,subplot(numStep+1,1,i),
    plot(datatothumbprint)

end

v_a0 = [3.177945
        3.174709
        3.168145
        3.142365
        3.086967
        2.995234
        2.970463
        2.923188
        2.879336
        2.787392
        2.724847];
v_s0 = [5.092499
        5.130234
        5.164895
        5.252934
        5.32018
        5.365323
        5.374716
        5.38747
        5.395083
        5.405279
        5.414055];
    
s0 = (400./v_s0)+(400/v_s0(1));
a0 = (400./v_a0)+(400/v_a0(1));

