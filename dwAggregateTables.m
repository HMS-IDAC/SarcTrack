clear, clc

% aggregates only stats tables

% -------------------------
% parameter

fpath = '~/Downloads/SarcTrackSampleVideos';


% -------------------------
% read tables

mpaths = listfiles(fpath,'.avi');

As = [];
PFs = cell(length(mpaths),3);
for i = 1:length(mpaths)
    [~,n] = fileparts(mpaths{i});
    disp(n)
    T = readtable([fpath filesep n '_DWStats.csv']);
    A = table2array(T);
    As = [As; A];

    PFs{i,1} = n;
    aPF = table2array(readtable([fpath filesep n '_DWPrdFrq.csv']));
    if ~isempty(aPF)
        PFs{i,2} = aPF(1);
        PFs{i,3} = aPF(2);
    end
end

if ~isempty(As)
    % -------------------------
    % write aggregated tables

    T = array2table(As,'VariableNames',T.Properties.VariableNames);
    writetable(T,[fpath filesep '_AggDWStats.csv']);
    writetable(cell2table(PFs,'VariableNames',{'Video','Period','Frequency'}),[fpath filesep '_AggDWPrdFrq.xls']);
    
    % -------------------------
    % bar plots

    prms = As';
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