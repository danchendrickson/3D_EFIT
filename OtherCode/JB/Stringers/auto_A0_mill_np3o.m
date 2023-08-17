pathname = 'D:\Jill\NDE\Projects\Oceana\Corey Tests\Milling\averaged\';

dist = .8;
delay = 0;
numStep = 10;
% numStep = 14;
rs = 3;

steps = [1.5982,1.539,1.4544,1.2402,1.0048,...
    0.786,0.7546,0.683,0.6206,0.5382,0.488];

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
        rawdata = rawdata-mean(rawdata);
        rawdata = resample(rawdata,1,rs);
        % -- filt data
            wvtpf = 'coif5';        % Wavelet Name
            swatoremove = [];       % aproxamation levels to remove
            swdtoremove = [1:4];      % detail levels to remove
            numlvls = 5;            % Number of Levels to Use

        rawdata = rawdata(1:(length(rawdata)-rem(...
            length(rawdata),2^numlvls)));
        % clip raw data to appropriate size for wavelet transform
        [swa,swd] = swt(rawdata, numlvls, wvtpf);
        % stationary wavelet transform
        swa(swatoremove,:)=0;       
        % remove some approxamations 
        swd(swdtoremove,:)=0;  
        % remove some details 
        filtdata = iswt(swa, swd, wvtpf); 
        % inverse stationary wavelet transform 

        % -- thumbprint
            wtpwidth = 1500;    % Window width
            wtpleft = floor(AO-100);
            wtpright = wtpwidth+wtpleft-1;
            wvttp = 'morl';     % Wavelet Name
            ns = 75;            % Number of Levels to Use
            nr = 5;             % Number of Ridges
            rw = .1;            % Ridge Width
            gv = .6;            % Valley multiplier
            op = 0;             % Show: both(0),valleys(1),peaks(2)

        datatothumbprint = filtdata(wtpleft:wtpright);
        datatothumbprint = datatothumbprint.*tukeywin(...
            length(datatothumbprint),.25)';

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
        for k = 1:wtpwidth
            v = 0;
            p = 0;
            tp = 0;
            tv = 0;
            for j = 1:ns
                if ((thumbprint(j,k)==1) && (thumbprint(j,k) ~= p))
                    tp = tp+1;
                    p = 1;
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
            if ((tp >= 3))% || (tv >= 3))
                arrival = k;
                break;
            end
        end 
        
        arr(i) = (arrival+wtpleft)/(sr/rs);
        
        figure(1),hold on,subplot(numStep+1,1,i),imshow(thumbprint)
        hold on
        line(ones(75,1)*arrival,(1:75),'Color','b','LineWidth',2)
        title(['A0 arrival for Mill thickness ' num2str(steps(i)) 'mm'])
        YLabel([num2str(arr(i)-3) 'us'],'Rotation',0)
        
        
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

a0 = (400./v_a0)+(400/v_a0(1));
figure(3), plot(steps,a0,'ko','MarkerSize',5,'MarkerFaceColor','k')
hold on, plot(steps,arr-3,'bo','MarkerSize',5,'MarkerFaceColor','b')