function make_fingers(filtdata,i)%file)

%fid = fopen(file, 'r');
%wave = fread(fid,'int16');
%fclose(fid);

    % -- plot thumbprint
    wvttp = 'haar';	% Wavelet Name
    ns = 200;           % Number of Levels to Use
    nr = 15;            % Number of Ridges
    rw = .05;           % Ridge Width
    gv = .6;            % Valley multiplier

    datatothumbprint = filtdata;
    datatothumbprint = datatothumbprint.*tukeywin(length(...
        datatothumbprint),.25);%';
    
    % get thumbprint for peaks
    thumbprintpeaks   = getthumbprint( datatothumbprint, ...
        wvttp, ns, (1), nr, rw, 2 );
    % get thumbprint for valleys
    thumbprintvalleys = getThumbprint( datatothumbprint, ...
        wvttp, ns, (1), nr, rw, 3 );  
    thumbprint = thumbprintpeaks + thumbprintvalleys.*gv;
    figure(4), hold on
    subplot(5,1,i)
    imshow(thumbprint)
    if i == 1
        title('DFWT Thumbprints for Increasing Propagation Lengths');
    end
    
    arrival1 = 0;
    arrival2 = 0;
     for k = 2100:length(filtdata)
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
         
    line(ones(200,1)*arrival1,(1:200),'Color','r','LineWidth',3)
    line(ones(200,1)*arrival2,(1:200),'Color','b','LineWidth',3)
