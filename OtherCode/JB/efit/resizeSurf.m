for i = 1:114
    rs(:,i) = resample(A.thickness(1:880,i),329,880,0);
end
% figure, contourf(rs)
for i = 1:329
    rs2(i,:) = resample(rs(i,:),229,114,0);
end
% figure, contourf(rs2(1:3280,:))

clear rs
ds = 1.22e-4;

rs3 = rs2(1:328,:);
rs3 = rs3(:,1:228);
for i = 1:328
    for j = 1:228
        rs4(i,j) = ceil(rs3(i,j)/ds)+1;
    end
end

rs5 = zeros(328,228);
for x = 1:328
    for y = 1:228   
        rs5(x,y) = rs4(x,228-y+1);
    end
end

clear rs2
figure(1),contourf(rs5')

H = fspecial('disk',5);
rs5 = imfilter(rs5,H,'replicate');
figure(6), contourf(rs5')

rs5 = rs5 +1;
rs5(1,:) = 15;
rs5(328,:) = 15;
rs5(:,228) = 15;

s2 = rs5(1:166,1:228);%*0+8;
size(s2)
s2 = [ones(1,228)*15;s2];
s2 = [ones(167,1)*0,s2];
size(s2)
s2f = ceil(s2);
figure(2),contourf(s2f',[6:15])
csvwrite('surf_2',s2f')

s4 = rs5(165:328,1:228);%*0+13;
size(s4)
s4 = [s4;ones(1,228)*15];
s4 = [ones(165,1)*0,s4];
s4 = [s2f((166:167),:);s4];
size(s4)
s4f = ceil(s4);
figure(4),contourf(s4f',[6:15])
csvwrite('surf_4',s4f')