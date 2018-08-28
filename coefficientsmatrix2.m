function [M,A,X,Y,kernels,iKernels] = coefficientsmatrix2(I,d,nangs,stretch,scale,hopsize,halfwindowsize,magthreshold,convtype)
% I should be double, and in the range [0,1]

orientations = (0:nangs-1)*180/nangs;
norientations = length(orientations);

windowsize = 2*halfwindowsize;

nr = min([floor((size(I,1)-windowsize)/hopsize)+1 size(I,1)]);
nc = min([floor((size(I,2)-windowsize)/hopsize)+1 size(I,2)]);

C = zeros(nr,nc);
A = zeros(nr,nc);
X = zeros(nr,nc);
Y = zeros(nr,nc);

RI3 = zeros(size(I,1),size(I,2),norientations);

mr = smorlet2(d,stretch,scale,0);
PI = padarray(I,(size(mr)-1)/2,'replicate');
% PI = I;
% PI = I;
% PI = gpuArray(I);
kernels = cell(1,norientations);
for j = 1:norientations
    orientation = orientations(j);

    mr = smorlet2(d,stretch,scale,orientation);
    
    if strcmp(convtype,'real')
        RI3(:,:,j) = conv2(PI,mr,'valid');
%         RI3(:,:,j) = conv2(PI,mr,'same');
%         RI3(:,:,j) = imfilter(PI,mr,'same');
%         RI3(:,:,j) = gather(imfilter(PI,mr,'same'));
        kernel = mr;
    elseif strcmp(convtype,'real+')
        R = conv2(PI,mr,'valid');
%         R = conv2(PI,mr,'same');
%         R = imfilter(PI,mr,'same');
%         R = gather(imfilter(PI,mr,'same'));
        RI3(:,:,j) = R.*(R > 0);
        kernel = mr.*(mr > 0);
    end
    
    kernels{j} = kernel;
end
rows = round(linspace(halfwindowsize+1,size(I,1)-halfwindowsize,nr));
cols = round(linspace(halfwindowsize+1,size(I,2)-halfwindowsize,nc));

iKernels = zeros(nr,nc);
for k = 1:nr
    for l = 1:nc
        row = rows(k);
        col = cols(l);
        M3 = RI3(row-halfwindowsize:row+halfwindowsize,col-halfwindowsize:col+halfwindowsize,:);
        [M,IM] = max(M3,[],3);
        [mC,imC] = max(M);
        [mR,imR] = max(mC);
        C(k,l) = mR;
        rm = imC(imR);
        cm = imR;
        A(k,l) = (IM(rm,cm)-1)*pi/nangs+pi/2;
        X(k,l) = row+rm-(halfwindowsize+1);
        Y(k,l) = col+cm-(halfwindowsize+1);
        iKernels(k,l) = IM(rm,cm);
    end
end

M = C;
% bw = imregionalmax(C);
% if isempty(magthreshold)
%     M = C.*bw.*imbinarize(C);
% else
%     M = C.*bw.*(C > magthreshold*max(C(:)));
% end

end