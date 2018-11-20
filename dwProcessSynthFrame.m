%% parameters

frameIndex = 1;

ds = 9:0.2:11; % range of distances (in pixels)
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8;
hopsize = 7;
halfwindowsize = 2;
magthreshold = [];
convtype = 'real+';

%% read frame, show wavelets

I = S(:,:,frameIndex);

mr1 = normalize(smorlet2(ds(1),stretch,scale,90));
mr2 = normalize(smorlet2(ds(end),stretch,scale,90));
J = normalize(I);
J(1:size(mr1,1),1:size(mr1,2)) = mr1;
J(size(mr1,1)+1:size(mr1,1)+size(mr2,1),1:size(mr2,2)) = mr2;
imtool(J)

%% process

tic
[M3,A3,X3,Y3] = coefficientsmatrix23(adapthisteq(I),ds,...
    'NumAngles',nangs,...
    'WavStretch',stretch,...
    'WavScale',scale,...
    'MagThreshold',magthreshold,...
    'HopSize',hopsize,...
    'HalfWindowSize',halfwindowsize,...
    'ConvType',convtype);
toc

%% display

V = zeros(size(M3,1),size(M3,2));
A = zeros(size(V));
X = zeros(size(V));
Y = zeros(size(V));
[M,V] = max(M3,[],3);

for row = 1:size(V,1)
    for col = 1:size(V,2)
        A(row,col) = A3(row,col,V(row,col));
        X(row,col) = X3(row,col,V(row,col));
        Y(row,col) = Y3(row,col,V(row,col));     
    end
end

m = reshape(M,1,[]);
l = reshape(V,1,[]);
a = reshape(A,1,[]);
x = reshape(X,1,[]);
y = reshape(Y,1,[]);
dsl = ds(l);

C = cool(size(I,1));
CB = zeros(size(I,1),20,3);
for i = 1:size(I,1)
    c = C(size(I,1)-i+1,:);
    CB(i,:,:) = repmat(reshape(c,[1 1 3]),[1 size(CB,2)]);
end
CB = [0.5*ones(size(I,1),5,3) CB];
C = cool(length(ds));
F = repmat(I,[1 1 3]);

for j = 1:length(l)
    if m(j) > 0.5
        d = dsl(j);
        c = C(l(j),:);
        row0 = x(j);
        col0 = y(j);
        angle = a(j);
        for k = -d/2:d/2% [-d/2:-d/4 d/4:d/2]
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
        for k = -d/4:d/4
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
end
figure
imshow([F CB])