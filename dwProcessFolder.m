clear, clc

% -------------------------
% parameter

fpath = '/home/mc457/files/CellBiology/IDAC/Marcelo/Seidman/Toepfer/DemoVideoCrops';

% -------------------------
% batch processing

mpaths = listfiles(fpath,'.avi');
for i = 1:length(mpaths)
    fprintf('\n')
    disp(mpaths{i});
    dwProcessVideo(mpaths{i});
end