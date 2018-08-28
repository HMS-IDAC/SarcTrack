function dwProcessVideo(path)

% clear, clc
tic

%% parameters

ds = 10:0.2:15;
stretch = 0.5; % stretch of morlet wavelet
scale = 2; % scale of morlet wavelet
nangs = 8;
hopsize = 5;
halfwindowsize = 2;
magthreshold = [];
convtype = 'real+';

%% read frames
disp('reading frames')

v = VideoReader(path);
nFrames = round(v.Duration*v.FrameRate);
Frames = cell(1,nFrames);
count = 0;
while hasFrame(v)
    if mod(count+1,round(nFrames/10)) == 1
        fprintf('.')
    end
    count = count+1;
    frame = readFrame(v);
    I = double(rgb2gray(frame))/255;
    Frames{count} = I;
end
fprintf('\n')

%% compute double-wavelet outputs
disp('fitting double-wavelet kernels')
% disp('---')
% disp('parallel processing / difficult to estimate time')
% disp('---')

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
bw = s > 0.25;

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
[fpdPath,fpdName] = fileparts(path);
save([fpdPath filesep fpdName '_fpd.mat'], 'fpd');


%% fit sawtooth curves
disp('fitting sawtooth curves')

[r,c] = find(bw);
prms = [];
idcs = [];
dsls = [];
for index = 1:length(r)
    if mod(index,round(length(r)/10)) == 1
        fprintf('.')
    end
    l = squeeze(V(r(index),c(index),:));

    x = (0:nFrames-1)';
    dsl = ds(l)';
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
        end
    end
end
fprintf('\n')

%% draw outputs
disp('writing images')

[rpath,fname] = fileparts(path);
outFit = [rpath filesep fname '_DWFit'];
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

%% write tables
disp('writing tables')

[rpath,fname] = fileparts(path);
outPathS = [rpath filesep fname '_DWStats.csv'];
outPathD = [rpath filesep fname '_DWDists.csv'];
outPathPF = [rpath filesep fname '_DWPrdFrq.csv'];


if ~isempty(prms)
    prd = f.b1*nFrames;
    c = prms(1,:)/(2*pi)*prd;
    r = prms(2,:)/(2*pi)*prd;
    o = prms(3,:)/(2*pi)*prd;

    T = array2table([c' r' o' prms(4:7,:)'],'VariableNames',{'contraction_time','relaxation_time','offset_from_average','min_ds','max_ds','min_ds_fit','max_ds_fit'});
    writetable(T,outPathS);
    
    vn = cell(1,size(dsls,2));
    for i = 1:size(dsls,2)
        vn{i} = sprintf('track%05d',idcs(i));
    end
    T = array2table(dsls,'VariableNames',vn);
    writetable(T,outPathD);
    
    T = array2table([prd 1/prd],'VariableNames',{'period','frequency'});
    writetable(T,outPathPF);
else
    writetable(array2table([]),outPathS);
    writetable(array2table([]),outPathD);
end

%%
toc

end