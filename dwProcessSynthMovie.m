%% parameters

ds = 9:0.2:11;
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8;
hopsize = 7;
halfwindowsize = 2;
magthreshold = [];
convtype = 'real+';

%% read frames (from memory)
disp('reading frames')

nFrames = size(S,3);
Frames = cell(1,nFrames);
for i = 1:nFrames
    Frames{i} = S(:,:,i);
end
disp('done')

%% read frames (from disk)
% disp('reading frames')
% 
% v = VideoReader('~/Downloads/SarcTrackSampleVideos/Synth/Sample_10_22.avi');
% nFrames = round(v.Duration*v.FrameRate);
% Frames = cell(1,nFrames);
% count = 0;
% while hasFrame(v)
%     if mod(count+1,round(nFrames/10)) == 1
%         fprintf('.')
%     end
%     count = count+1;
%     frame = readFrame(v);
%     I = double(rgb2gray(frame))/255;
%     Frames{count} = I;
% end
% fprintf('\n')

%% compute double-wavelet outputs
disp('fitting double-wavelet kernels')

cM = cell(1,nFrames);
cA = cell(1,nFrames);
cX = cell(1,nFrames);
cY = cell(1,nFrames);
pfpb = pfpbStart(nFrames);
parfor fIndex = 1:nFrames
    pfpbUpdate(pfpb);
    [M3,A3,X3,Y3] = coefficientsmatrix23(adapthisteq(Frames{fIndex}),ds,...
        'NumAngles',nangs,...
        'WavStretch',stretch,...
        'WavScale',scale,...
        'MagThreshold',magthreshold,...
        'HopSize',hopsize,...
        'HalfWindowSize',halfwindowsize,...
        'ConvType',convtype);

    cM{fIndex} = M3;
    cA{fIndex} = A3;
    cX{fIndex} = X3;
    cY{fIndex} = Y3;
end

%% display fit

% for i = 1:nFrames
%     disp(i)
%     D = dwDrawFit(cM{i},cA{i},cX{i},cY{i},ds,adapthisteq(Frames{i}));
%     imshow(D)
%     pause(0.1)
% end

%% aggregate outputs
disp('aggregating outputs')

V = zeros(size(cM{1},1),size(cM{1},2),nFrames);
A = zeros(size(V));
X = zeros(size(V));
Y = zeros(size(V));
mx = zeros(size(V,1),size(V,2));
for i = 1:nFrames
    [m,im] = max(cM{i},[],3);
    V(:,:,i) = im;
    mx = max(m,mx);
    
    A3 = cA{i};
    X3 = cX{i};
    Y3 = cY{i};
    for row = 1:size(im,1)
        for col = 1:size(im,2)
            A(row,col,i) = A3(row,col,im(row,col));
            X(row,col,i) = X3(row,col,im(row,col));
            Y(row,col,i) = Y3(row,col,im(row,col));     
        end
    end
end

M = imbinarize(mx);
s = std(ds(V).*repmat(M,[1 1 size(V,3)]),0,3);
bw = s > 0.1;
% imshow(bw)
disp('done')

%% inspect single 'track'

[r,c] = find(bw);
index = 56; % should be between 1 and length(r)
l = squeeze(V(r(index),c(index),:));
dsl = ds(l)';

for iFrame = 1:nFrames
    I = repmat(S(:,:,iFrame),[1 1 3]);
    row0 = X(r(index),c(index),iFrame);
    col0 = Y(r(index),c(index),iFrame);
    angle = A(r(index),c(index),iFrame);
    d = ds(V(r(index),c(index),iFrame));
    for k = -d/2:d/2
        row = round(row0+k*cos(angle+pi/2));
        col = round(col0+k*sin(angle+pi/2));
        if row >= 1 && row <= size(I,1) && col >= 1 && col <= size(I,2)
            I(row,col,:) = 0;
            I(row,col,2) = 1;
        end
    end
    subplot(1,2,1)
    imshow(I)
    subplot(1,2,2)
    plot(1:nFrames,dsl,'.'), hold on
    plot(iFrame,dsl(iFrame),'o'), hold off
    pause(0.1)
end

%% estimate frequency
disp('estimating frequency')

[r,c] = find(bw);
avgDsts = zeros(nFrames,1);
for index = 1:length(r)
    l = squeeze(V(r(index),c(index),:));
    dsl = ds(l)';
    avgDsts = avgDsts+dsl;
end
avgDsts = avgDsts/length(r);

x = (0:nFrames-1)';
mx = prctile(avgDsts,90);
mn = prctile(avgDsts,10);
yAvg = (avgDsts-mean(avgDsts))/std(avgDsts);
y2Fit = 2*((avgDsts-mn)/(mx-mn)-0.5);

syAvg = smooth(yAvg,15);
f = fit(x,syAvg,'sin1');
ySin = f.a1*sin(f.b1*x+f.c1);

c = pi/2;
r = pi/2;
o = 0;
ySaw = f.a1*sawtooth(f.b1*x+pi/2+f.c1,c,r);

f2m = @(cro) -corr( f.a1*sawtooth(f.b1*x+cro(3)+pi/2+f.c1,cro(1),cro(2)) , y2Fit );

cro0 = [c; r; o];
lb = [0 0 -pi/4];
ub = [pi pi pi/4];
cro = fmincon(f2m,cro0,[],[],[],[],lb,ub,[],optimoptions('fmincon','Display','off'));
ySawFit = f.a1*sawtooth(f.b1*x+cro(3)+pi/2+f.c1,cro(1),cro(2));

% f plot data, used by dwCheckResults
fpd.x = x;
fpd.y2Fit = y2Fit;
fpd.ySin = ySin;
fpd.ySaw = ySaw;
fpd.ySawFit = ySawFit;

figureQSS
plot(fpd.x,fpd.y2Fit,'r'), hold on
plot(fpd.x,fpd.ySin,'g'),
plot(fpd.x,fpd.ySaw,'b'),
plot(fpd.x,fpd.ySawFit,'k'), hold off
legend('y2Fit','ySin','ySaw', 'ySawFit')


%% fit sawtooth curves
disp('fitting sawtooth curves')

[r,c] = find(bw);
prms = [];
idcs = [];
dsls = [];
figureQSS
for index = 1:length(r)
    if mod(index,round(length(r)/10)) == 1
        fprintf('.')
    end
    l = squeeze(V(r(index),c(index),:));

    x = (0:nFrames-1)';
    dsl = ds(l)';
    
%     plot(dsl)
%     pause(0.1)
%     continue
    
    mx = prctile(dsl,90);
    mn = prctile(dsl,10);
    if mx > mn
        y = 2*((dsl-mn)/(mx-mn)-0.5);
        f2m = @(prm) -corr( f.a1*sawtooth(f.b1*x+prm(3)+pi/2+f.c1,prm(1),prm(2)) , y );
        options = optimoptions('fmincon','Display','off');
        prm = fmincon(f2m,cro,[],[],[],[],lb,ub,[],options);
        ySawFit = f.a1*sawtooth(f.b1*x+prm(3)+pi/2+f.c1,prm(1),prm(2));
        z = ((ySawFit/2)+0.5)*(mx-mn)+mn; % original range

        rmse = sqrt(sum((y-ySawFit).^2)/nFrames);

        if rmse < 1 && max(dsl) < max(ds) && min(dsl) > min(ds)
            prms = [prms [prm; min(dsl); max(dsl); min(z); max(z)]];
            idcs = [idcs index];
            dsls = [dsls dsl];
            
            plot(fpd.x,y,'r'), hold on
            plot(fpd.x,ySawFit,'g'), hold off
            legend('y2Fit','ySawFit')
            pause(0.1)
        end
    end
end
close all
fprintf('\n')

%% draw outputs
disp('writing images')

outFit = '~/Downloads/Outputs';
if ~exist(outFit,'dir')
    mkdir(outFit);
end

l = zeros(nFrames,length(idcs));
dsl = zeros(nFrames,length(idcs));
a = zeros(nFrames,length(idcs));
x = zeros(nFrames,length(idcs));
y = zeros(nFrames,length(idcs));
for i = 1:length(idcs)
    index = idcs(i);
    l(:,i) = squeeze(V(r(index),c(index),:));
    dsl(:,i) = ds(l(:,i));
    a(:,i) = squeeze(A(r(index),c(index),:));
    x(:,i) = squeeze(X(r(index),c(index),:));
    y(:,i) = squeeze(Y(r(index),c(index),:));
end

I = Frames{1};
C = cool(size(I,1));
CB = zeros(size(I,1),20,3);
for i = 1:size(I,1)
    c = C(size(I,1)-i+1,:);
    CB(i,:,:) = repmat(reshape(c,[1 1 3]),[1 size(CB,2)]);
end
CB = [0.5*ones(size(I,1),5,3) CB];
C = cool(length(ds));
for fIndex = 1:nFrames
    if mod(fIndex,round(nFrames/10)) == 1
        fprintf('.');
    end
    F = Frames{fIndex};
    F = repmat(F,[1 1 3]);
    for j = 1:length(idcs)
        d = dsl(fIndex,j);
        c = C(l(fIndex,j),:);
        row0 = x(fIndex,j);
        col0 = y(fIndex,j);
        angle = a(fIndex,j);
        for k = -d/2:d/2
            row = round(row0+k*cos(angle+pi/2));
            col = round(col0+k*sin(angle+pi/2));
            if row >= 1 && row <= size(F,1) && col >= 1 && col <= size(F,2)
                F(row,col,:) = 0;
                F(row,col,:) = 0.5*reshape(c,[1 1 3]);
            end
        end
        row1 = row0+d/2*cos(angle+pi/2);
        col1 = col0+d/2*sin(angle+pi/2);
        row2 = row0-d/2*cos(angle+pi/2);
        col2 = col0-d/2*sin(angle+pi/2);
        for k = -d/2:d/2
            row = round(row1+k*cos(angle));
            col = round(col1+k*sin(angle));
            if row >= 1 && row <= size(F,1) && col >= 1 && col <= size(F,2)
                F(row,col,:) = 0;
                F(row,col,:) = reshape(c,[1 1 3]);
            end
            row = round(row2+k*cos(angle));
            col = round(col2+k*sin(angle));
            if row >= 1 && row <= size(F,1) && col >= 1 && col <= size(F,2)
                F(row,col,:) = 0;
                F(row,col,:) = reshape(c,[1 1 3]);
            end
        end
    end
    imwrite([F CB],[outFit filesep sprintf('Frame%03d.png',fIndex)]);
end
fprintf('\n')

%% plot stats

prms0 = prms;
if ~isempty(prms0)
%     prd = f.b1*nFrames; % incorrect period
    prd = 2*pi/f.b1; % correct period
    c = prms0(1,:)/(2*pi)*prd;
    r = prms0(2,:)/(2*pi)*prd;
    o = prms0(3,:)/(2*pi)*prd;
    
    T = array2table([c' r' o' prms0(4:7,:)'],'VariableNames',{'contraction_time','relaxation_time','offset_from_average','min_ds','max_ds','min_ds_fit','max_ds_fit'});
    A = table2array(T);
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
    
    
    
    vn = cell(1,size(dsls,2));
    for i = 1:size(dsls,2)
        vn{i} = sprintf('track%05d',idcs(i));
    end
    T = array2table(dsls,'VariableNames',vn);
    A = table2array(T);
    figure
    histogram(A(:),20), title('dists')

    T = array2table([prd 1/prd],'VariableNames',{'period','frequency'});
    disp(T)

end

%% plot average fit curve, compare to ground truth

hPad = 10;
p = ones(hPad,1);
t = [p*A(1,1); A(:,1); p*A(end,1)];
at = t;
for i = 2:size(A,2)
    cs = zeros(1,2*hPad);
    for j = 1:2*hPad
%         plot(t,'r'), hold on
        u = [A(1,i)*ones(j,1); A(:,i); A(end,i)*ones(2*hPad-j,1)];
%         plot(u,'g'), hold off
%         pause
        cs(j) = corr(t,u);
    end
    [~,j] = max(cs);
    u = [A(1,i)*ones(j,1); A(:,i); A(end,i)*ones(2*hPad-j,1)];
%     plot(t,'r'), hold on
%     plot(u,'g'), hold off
%     pause
    at = at+u;
end
at = at/size(A,2);
% plot(at)


i0 = round((length(at)-length(dsts))/2);
at = at(i0:i0+length(dsts)-1);

t = [p*dsts(1); dsts'; p*dsts(end)];
cs = zeros(1,2*hPad);
for j = 1:2*hPad
    u = [at(1)*ones(j,1); at; at(end)*ones(2*hPad-j,1)];
    cs(j) = corr(t,u);
end
[~,j] = max(cs);
u = [at(1)*ones(j,1); at; at(end)*ones(2*hPad-j,1)];

plot(t,'g'), hold on
plot(u,'b'), hold off
axis([-1 length(t)+1 9.25 10.75])
legend(sprintf('gt (%d, %d)', f1, f2),sprintf('avg fit (%d)',size(A,2)))