
pathname = ['D:\Jill\NDE\Projects\Oceana\Corey Tests\' ...
    'Corrosion\T Stringers\averaged\'];

dist = .8;
delay = 0;
% numStep = 10;
numStep = 14;
rs = 3;

steps = [1.5906;
         1.5998;
         1.5968;
         1.6106;
         1.6238;
         1.67;
         1.6928;
         1.723;
         1.8446;
         1.8768;
         1.8986;
         1.9892;
         2.1358;
         2.2002];
arrt = [5:13];
arr2 = [343 322 352 345 431 472 494 499 555];
t = 1;

for i = 1:numStep
    if i >= 10 
        filename{i} = [num2str(i),'-R\T80NP3O',num2str(i),'R'];
    else
        filename{i} = ['0',num2str(i),'-R\T80NP3O0',num2str(i),'R'];
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
    for i = 1:numStep;
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

        rawdata = rawdata(1:(length(rawdata)-rem( ...
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
        datatothumbprint = datatothumbprint.*tukeywin( ...
            length(datatothumbprint),.25)';

        thumbprintpeaks   = getThumbprint( datatothumbprint,...
            wvttp, ns, (1), nr, rw, 2 );  % get thumbprint for peaks
        thumbprintvalleys = getThumbprint( datatothumbprint, ...
            wvttp, ns, (1), nr, rw, 3 );  % get thumbprint for valleys

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
        
%         figure out arrival in micro seconds
        arr(i) = (arrival+wtpleft)/(sr/rs);
        
        if (i > 7)
            figure(2),subplot(7,1,i-7),hold off
        else
            figure(1),subplot(7,1,i),hold off
        end
        plot([1:length(rawdata)]/(sr/rs),rawdata,'k')
        axis([0 500 -1000 1000])
        
        if (i > 7)
            figure(4),subplot(7,1,i-7),hold off
        else
            figure(3),subplot(7,1,i),hold off
        end
        plot([1:length(filtdata)]/(sr/rs),filtdata,'k')
        hold on,plot([wtpleft:wtpright]/(sr/rs),...
            datatothumbprint,'r'),hold off
        axis([0 500 -500 500])
        
        if (i > 7)
            figure(6),subplot(7,1,i-7),hold off
        else
            figure(5),subplot(7,1,i),hold off
        end
        imshow(thumbprint(:,(1:1500)))
        hold on
        line(ones(75,1)*arrival,(1:75),'Color','b','LineWidth',2)
        YLabel([num2str(arr(i)-3) 'us'],'Rotation',0)
        if (i >= 5 && i~=14)
            line(ones(75,1)*arr2(t),(1:75),'Color','r','LineWidth',2)
            t = t+1;
        end

        
        
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
figure(7)
hold on, plot([1:14],arr-3,'bo','MarkerSize',5,'MarkerFaceColor','b')
plot(arrt,(arr2+wtpleft)/(sr/rs),'ro','MarkerSize',5,'MarkerFaceColor','r')



% times = [252.1 255.5 255.6 254.8 260.2 259.4 260.5 260.3 263.5 265 ...
%     265.8 266 268.1 269.4];
% vel = 400./(times - (400/3.177945));
% figure, plot([1:14],vel,'bo','MarkerSize',5,'MarkerFaceColor','b')