for i = 1:114
    rs(:,i) = resample(A.thickness(1:880,i),3290,880,0);
end
% figure, contourf(rs)
for i = 1:3290
    rs2(i,:) = resample(rs(i,:),229,114,0);
end
% figure, contourf(rs2(1:3280,:))

clear rs
ds = 1.22e-4;

rs3 = rs2(4:3284,:);
rs3 = rs3(:,1:228);
for i = 1:3281
    for j = 1:228
        rs4(i,j) = ceil(rs3(i,j)/ds)+1;
    end
end

rs5 = zeros(3281,228);
for x = 1:3281
    for y = 1:228   
        rs5(x,y) = rs4(x,228-y+1);
    end
end


clear rs2 rs3
figure(1),contourf(rs5')

H = fspecial('disk',10);
rs6 = imfilter(rs5,H,'replicate');
figure(2), contourf(rs6')

rs6 = rs6 +1;
rs6(1,:) = 15;
rs6(3281,:) = 15;
rs6(:,228) = 15;


% % % s = rs4(1:167,1:120);
    s2 = rs6(1:165,1:121);
    s2 = [ones(1,121)*15;s2];
    s2 = [ones(166,1)*0,s2];
    s28 = size(s2)
    s2f = ceil(s2);
    figure(3),contourf(s2f',[6:15])
    csvwrite('surf_28',s2f')
    clear s2

% % % s = rs4(163:494,1:120);
    s4 = rs6(164:492,1:121);
    s4 = [ones(329,1)*0,s4];
    s4f = [s2f((163:164),:);s4];
    s32 = size(s4f)
    s4f = ceil(s4f);
    figure(4),contourf(s4f',[6:15])
    csvwrite('surf_32',s4f')
    clear s4
    
% % % s = rs4(493:821,1:120);
    s2 = rs6(493:821,1:121);
    s2 = [ones(329,1)*0,s2];
    s2f = [s4f((330:331),:);s2];
    s36 = size(s2f)
    s2f = ceil(s2f);
    figure(5),contourf(s2f',[6:15])
    csvwrite('surf_36',s2f')
    clear s2
    
% % % s = rs4(819:1149,1:120);
    s4 = rs6(822:1150,1:121);
    s4 = [ones(329,1)*0,s4];
    s4f = [s2f((330:331),:);s4];
    s40 = size(s4f)
    s4f = ceil(s4f);
    figure(6),contourf(s4f',[6:15])
    csvwrite('surf_40',s4f')
    clear s4
    
% % % s = rs4(1148:1476,1:120);
    s2 = rs6(1151:1479,1:121);
    s2 = [ones(329,1)*0,s2];
    s2f = [s4f((330:331),:);s2];
    s44 = size(s2f)
    s2f = ceil(s2f);
    figure(7),contourf(s2f',[6:15])
    csvwrite('surf_44',s2f')
    clear s2
    
% % % s = rs4(1476:1804,1:120);
    s4 = rs6(1480:1808,1:121);
    s4 = [ones(329,1)*0,s4];
    s4f = [s2f((330:331),:);s4];
    s48 = size(s4f)
    s4f = ceil(s4f);
    figure(8),contourf(s4f',[6:15])
    csvwrite('surf_48',s4f')
    clear s4
    
% % % s = rs4(1804:2132,1:120);
    s2 = rs6(1809:2137,1:121);
    s2 = [ones(329,1)*0,s2];
    s2f = [s4f((330:331),:);s2];
    s52 = size(s2f)
    s2f = ceil(s2f);
    figure(9),contourf(s2f',[6:15])
    csvwrite('surf_52',s2f') 
    clear s2
    
% % % s = rs4(2132:2460,1:120);
    s4 = rs6(2138:2466,1:121);
    s4 = [ones(329,1)*0,s4];
    s4f = [s2f((330:331),:);s4];
    s56 = size(s4f)
    s4f = ceil(s4f);
    figure(10),contourf(s4f',[6:15])
    csvwrite('surf_56',s4f')
    clear s4
    
% % % s = rs4(2460:2788,1:120);
    s2 = rs6(2467:2795,1:121);
    s2 = [ones(329,1)*0,s2];
    s2f = [s4f((330:331),:);s2];
    s60 = size(s2f)
    s2f = ceil(s2f);
    figure(11),contourf(s2f',[6:15])
    csvwrite('surf_60',s2f')
    clear s2
    
% % % s = rs4(2788:3116,1:120);
    s4 = rs6(2796:3124,1:121);
    s4 = [ones(329,1)*0,s4];
    s4f = [s2f((330:331),:);s4];
    s64 = size(s4f)
    s4f = ceil(s4f);
    figure(12),contourf(s4f',[6:15])
    csvwrite('surf_64',s4f')
    clear s4 s2f
    
% % % s = rs4(3116:3279,1:120);
    s2 = rs6(3117:3281,1:121);
    s2 = [s2;ones(1,121)*14];
    s2f = [ones(166,1)*0,s2];
    s68 = size(s2f)
    s2f = ceil(s2f);
    figure(13),contourf(s2f',[6:15])
    csvwrite('surf_68',s2f')
    clear s2
    
% % % % % % % % % % % % % % % % % % % 
% % % s = rs4(1:167,1:120);
    s2 = rs6(1:165,120:228);
    s2 = [ones(1,109)*14;s2];
    s2 = [s2,ones(166,1)*15];
    s29 = size(s2)
    s2f = ceil(s2);
    figure(14),contourf(s2f',[6:15])
    csvwrite('surf_29',s2f')
    clear s2

% % % s = rs4(163:494,1:120);
    s4 = rs6(164:492,120:228);
    s4 = [s4,ones(329,1)*15];
    s4f = [s2f((163:164),:);s4];
    s33 = size(s4f)
    s4f = ceil(s4f);
    figure(15),contourf(s4f',[6:15])
    csvwrite('surf_33',s4f')
    clear s4
    
% % % s = rs4(493:821,1:120);
    s2 = rs6(493:821,120:228);
    s2 = [s2,ones(329,1)*15];
    s2f = [s4f((330:331),:);s2];
    s37 = size(s2f)
    s2f = ceil(s2f);
    figure(16),contourf(s2f',[6:15])
    csvwrite('surf_37',s2f')
    clear s2
    
% % % s = rs4(819:1149,1:120);
    s4 = rs6(822:1150,120:228);
    s4 = [s4,ones(329,1)*15];
    s4f = [s2f((330:331),:);s4];
    s41 = size(s4f)
    s4f = ceil(s4f);
    figure(17),contourf(s4f',[6:15])
    csvwrite('surf_41',s4f')
    clear s4
    
% % % s = rs4(1148:1476,1:120);
    s2 = rs6(1151:1479,120:228);
    s2 = [s2,ones(329,1)*15];
    s2f = [s4f((330:331),:);s2];
    s45 = size(s2f)
    s2f = ceil(s2f);
    figure(18),contourf(s2f',[6:15])
    csvwrite('surf_45',s2f')
    clear s2
    
% % % s = rs4(1476:1804,1:120);
    s4 = rs6(1480:1808,120:228);
    s4 = [s4,ones(329,1)*15];
    s4f = [s2f((330:331),:);s4];
    s49 = size(s4f)
    s4f = ceil(s4f);
    figure(19),contourf(s4f',[6:15])
    csvwrite('surf_49',s4f')
    clear s4
    
% % % s = rs4(1804:2132,1:120);
    s2 = rs6(1809:2137,120:228);
    s2 = [s2,ones(329,1)*15];
    s2f = [s4f((330:331),:);s2];
    s53 = size(s2f)
    s2f = ceil(s2f);
    figure(20),contourf(s2f',[6:15])
    csvwrite('surf_53',s2f') 
    clear s2
    
% % % s = rs4(2132:2460,1:120);
    s4 = rs6(2138:2466,120:228);
    s4 = [s4,ones(329,1)*15];
    s4f = [s2f((330:331),:);s4];
    s57 = size(s4f)
    s4f = ceil(s4f);
    figure(21),contourf(s4f',[6:15])
    csvwrite('surf_57',s4f')
    clear s4
    
% % % s = rs4(2460:2788,1:120);
    s2 = rs6(2467:2795,120:228);
    s2 = [s2,ones(329,1)*15];
    s2f = [s4f((330:331),:);s2];
    s61 = size(s2f)
    s2f = ceil(s2f);
    figure(22),contourf(s2f',[6:15])
    csvwrite('surf_61',s2f')
    clear s2
    
% % % s = rs4(2788:3116,1:120);
    s4 = rs6(2796:3124,120:228);
    s4 = [s4,ones(329,1)*15];
    s4f = [s2f((330:331),:);s4];
    s65 = size(s4f)
    s4f = ceil(s4f);
    figure(23),contourf(s4f',[6:15])
    csvwrite('surf_65',s4f')
    clear s4 s2f
    
% % % s = rs4(3116:3279,1:120);
    s2 = rs6(3117:3281,120:228);
    s2 = [s2;ones(1,109)*15];
    s2f = [s2,ones(166,1)*15];
    s69 = size(s2f)
    s2f = ceil(s2f);
    figure(24),contourf(s2f',[6:15])
    csvwrite('surf_69',s2f')
    clear s2
