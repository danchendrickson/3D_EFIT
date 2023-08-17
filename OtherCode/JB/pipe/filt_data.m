function filtdata = filt_data(rawdata)

resam_rawdata = resample(rawdata,1,20);
% -- filt data
    wvtpf = 'coif3';        % Wavelet Name
    swatoremove = [];       % aproxamation levels to remove
    swdtoremove = 1:3;      % detail levels to remove
    numlvls = 5;            % Number of Levels to Use

% clip raw data to appropriate size for wavelet transform
resam_rawdata = resam_rawdata(1:length(resam_rawdata)-rem(length  ...
    (resam_rawdata),2^numlvls));        
% stationary wavelet transform
[swa,swd] = swt(resam_rawdata, numlvls, wvtpf);	
% remove some approxamations 
swa(swatoremove,:)=0;    
% remove some details
swd(swdtoremove,:)=0;     
% inverse stationary wavelet transform 
filtdata = iswt(swa, swd, wvtpf);     

plot(filtdata,'r')