function wave = open_plot(file,kind,num,c)%,wave)

dir = 'D:\Jill\NDE\Projects\DOT pipes\Coatings data\';

fid = fopen([dir file],'r');
switch kind
    case 1  % 1 transmitter all receivers
        for j = 1:num
            w = fread(fid,30001,'int16');
            plot([1:30000]*20e-6,w,c)
            axis([0 .0015 -600 600])
            [j]
            pause(.1)
            wave{j} = w;
        end
    case 2 % parallel waves
        figure
        for j = 1:num
            for i = 1:num
                w = fread(fid,30001,'int16');
%                 size(w)
                if rem(j,2) && i==j
                    plot([1:30001]/20e6,w,c)
                    axis([0 .0015 -600 600])
                    [j,i]
                    pause(.01)
                    wave{j} = w;
                elseif ~rem(j,2) && i==num-j+1
                    plot([1:30001]/20e6,w,c)
                    axis([0 .0015 -600 600])
                    [j,i]
                    pause(.01)
                    wave{j} = w;
                end
            end
        end
    case 3  % single wave
        figure(num)
        hold on
        wave = fread(fid,'int16');
        plot(wave,c)
        axis([0 30000 -500 500])
%         pause(.01)
    case 4
        for i = 1:num
            wave{i} = fread(fid,30001,'int16');
            plot([1:30001]/20e6,wave{i},c)
            axis([0 .0015 -900 900])
            pause(.01)
        end        
end
fclose(fid);
