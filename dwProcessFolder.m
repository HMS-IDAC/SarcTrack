clear, clc

% -------------------------
% parameter

fpath = '~/Downloads/SarcTrackSampleVideos';

% -------------------------
% batch processing

mpaths = listfiles(fpath,'.avi');
for i = 1:length(mpaths)
    fprintf('\n')
    disp(mpaths{i});
    dwProcessVideo(mpaths{i});
end