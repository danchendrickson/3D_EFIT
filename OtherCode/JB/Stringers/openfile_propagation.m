% function openfile

%[file,path] = uigetfile('*.*');
file = '1400kHz_w';
path = 'D:\Jill\NDE\Projects\Oceana\tstiff\';
distance = [50 60 70 80 90];

for i = 1:5
    name = strcat(path,file,num2str(distance(i)))
    load(name);
    
    data = [zeros(2000,1); total];
        
    dist = distance(i)*.01;
    plot_filter
    figure(2)
    subplot(5,1,i), hold on
    plot([zeros(2000,1);d.line(:,2)],'g')
%    make_thumbnails(abs(filtdata),i)
%     make_fingers(abs(filtdata),i)
    make_fingers([zeros(2000,1);d.line(:,2)],i)
    
    
end


