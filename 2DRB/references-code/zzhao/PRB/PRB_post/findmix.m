clear 
clc
data= importdata('T.dat');
first_data = data(:, 1)';
second_data = data(:, 2);
[~, locs] = findpeaks(-second_data);  % 找到局部最小值的位置
disp(locs);  % 输出局部最小值的位置
itc(size(locs))=0;
for i = 1: size(locs)
    itc(i) = first_data(locs(i));
end

fileID = fopen('itc.dat', 'w');  % 打开文件，'w' 表示写入模式

% 将数据写入文件
for i = 1:size(locs)
    fprintf(fileID, '%d \n', itc(i));  % 每行写入3个浮点数，换行
end

fclose(fileID);  
