function [wave_out] = openfile_be(fname,size,num,label)
% This function solves the little endian vs big endian problem
%   Problem:
%       LabView uses the Big Endian: (00000000 00000001)
%       Matlab & C uses the Little Endian: (00000001 00000000)
%         In Big Endian format, the most significant byte (MSB) of a
%         multi-byte number is written first, then the second MSB, and
%         so on down to the least significant byte (LSB). However in 
%         Little Endian, it is reversed such that the LSB is written 
%         first, then the second LSB on up to the MSB.
%   Function:
%       wave_out = openfile_be(fname,size,num)
%         Opens fid for fname
%         Reads in 'int16' 'ieee-be' data
%         If num = 0, reads whole file and removes average
%         Else, reads in data of size and averages over num


fid = fopen(fname,'r');
if num == 0
    wave = fread(fid,'int16','ieee-be');
    wave_out = wave-mean(wave);
else
%     fseek(fid,15024000,'bof');
    wave = zeros(size,1);
    for i = 1:num
        w = fread(fid,size,'int16','ieee-be');
%         w = w - mean(w);
%         plot(w),pause(.001)
        wave = wave + w;
    end
%     wave = zeros(80240,1);
%     for i = 1:num
%         w = fread(fid,80240,'int16','ieee-be');
% %         w = w - mean(w);
%         plot(w),pause(.001)
%         wave = wave + w;
%     end
    wave_out = wave/num;
    clear w
end
fclose(fid);
clear wave
% figure,plot(wave_out)
% title(label);

% paramFourierTest

% filt_wave = ry;