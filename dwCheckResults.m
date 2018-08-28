clear, clc

% -------------------------
% parameter

path = '~/Downloads/SarcTrackSampleVideos/Sample2.avi';

% -------------------------
% average distance / frequency estimation

[fpath,fname] = fileparts(path);
load([fpath filesep fname '_fpd.mat']);
figureQSS
plot(fpd.x,fpd.y2Fit,'r'), hold on
plot(fpd.x,fpd.ySin,'g'),
plot(fpd.x,fpd.ySaw,'b'),
plot(fpd.x,fpd.ySawFit,'k'), hold off
legend('y2Fit','ySin','ySaw', 'ySawFit')

% -------------------------
% read/plot stats

tpath = [fpath filesep fname '_DWStats.csv'];
T = readtable(tpath);
A = table2array(T);

if ~isempty(A)
    % -------------------------
    % bar plots

    prms = A';
    labels = cell(1,7*size(prms,2));
    for i = 1:size(prms,2)
        labels{i} = 'contraction time';
        labels{size(prms,2)+i} = 'relaxation time';
        labels{2*size(prms,2)+i} = 'offset from average';
        labels{3*size(prms,2)+i} = 'min ds';
        labels{4*size(prms,2)+i} = 'max ds';
        labels{5*size(prms,2)+i} = 'min ds fit';
        labels{6*size(prms,2)+i} = 'max ds fit';
    end
    figureQSS
    boxplot([prms(1,:) prms(2,:) prms(3,:)],labels(1:3*size(prms,2)))
    figureQSS
    boxplot([prms(4,:) prms(5,:) prms(6,:) prms(7,:)],labels(3*size(prms,2)+1:end))

    % -------------------------
    % histograms

    figureQSS
    titles = {'contraction time','relaxation time','offset from average','min ds','max ds','min ds fit','max ds fit'};
    for i = 1:7
        subplot(1,7,i)
        histogram(prms(i,:))
        title(titles{i})
        fprintf('median %s: %f\n', titles{i}, median(prms(i,:)));
    end
    fprintf('median min/max: %f\n', median(prms(4,:)./prms(5,:)));
    fprintf('median min/max fit: %f\n', median(prms(6,:)./prms(7,:)));
else
    msgbox('Stats table is empty.');
end

% -------------------------
% read/plot dists

tpath = [fpath filesep fname '_DWDists.csv'];
T = readtable(tpath);
A = table2array(T);

if ~isempty(A)
    figure
    histogram(A(:),20), title('dists')
else
    msgbox('Dists table is empty.');
end