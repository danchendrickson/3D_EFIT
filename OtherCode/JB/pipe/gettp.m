% figure

% c3cp = open_plot('tomo data\new setup\c4\proj2',2,72,'k'); %HBE-95
% c3cp = open_plot('tomo data\new setup\c3\p2w',2,72,'k'); %Bitumastic 50
c3cp = open_plot('tomo data\new setup\c2\p2_1',2,72,'k'); %Protal+Bit

figure
for i = 1:72
    k3cp.resam{i} = resample(c3cp{i},1,20);
    k3cp.filtdata{i} = filt_data(c3cp{i});
    k3cp.data{i} = c3cp{i};
    pause(.01)
end


% for i = 1:72
% [k42.xx1{i} k42.peaks{i}] = find_peaks(k42.filtdata{i},42,'b');
% end
% for i = 1:72
% figure(221),hold on
% plot(i*ones(length(k22.xx1(i))),k22.xx1{i},'ko',  ...
%           'MarkerFaceColor','c','MarkerSize',5)
% end

figure
for i = 1:8
    for j = 1:9
        k3cp.tp{9*(i-1)+j} = thumb(k3cp.filtdata{9*(i-1)+j},i,j);
        pause(.01)
%         F(9*(i-1)+j) = getframe;
    end
end
% movie2avi(F,'m32t','compression','None')