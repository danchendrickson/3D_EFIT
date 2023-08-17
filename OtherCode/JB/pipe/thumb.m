function thumbprint = thumb(filtdata,i,j)
% -- thumbprint
    wtpwidth = 700;    % Window width
    wtpleft = 1;%floor(SO+1500);
    wtpright = wtpwidth+wtpleft-1;
    wvttp = 'gaus2';     % Wavelet Name
    ns = 50;            % Number of Levels to Use
    nr = 15;             % Number of Ridges
    rw = .03;            % Ridge Width
    gv = .6;            % Valley multiplier
    op = 2;             % Show: both(0),valleys(1),peaks(2)

datatothumbprint = filtdata(wtpleft:wtpright);
datatothumbprint=datatothumbprint.*tukeywin(length(datatothumbprint),.25)';

thumbprintpeaks   = getThumbprint( datatothumbprint, wvttp, ns,  ...
    (1), nr, rw, 2 );  % get thumbprint for peaks
thumbprintvalleys = getThumbprint( datatothumbprint, wvttp, ns,  ...
    (1), nr, rw, 3 );  % get thumbprint for valleys

if  op == 1                             % Show only Valleys
thumbprint = thumbprintvalleys.*gv;                             
elseif  op == 2                         % Show only Peaks
thumbprint = thumbprintpeaks;
else                                    % Show Both
thumbprint = thumbprintpeaks + thumbprintvalleys.*gv;
end
figure(i), subplot(9,1,j)
% plot(datatothumbprint)
imshow(thumbprint)
title(num2str(9*(i-1)+j))