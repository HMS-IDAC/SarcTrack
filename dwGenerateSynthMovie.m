%% parameters

dMax = 0.4; % maximum distance from contracted to relaxed, in pixels
f1 = 1; % contraction time factor
f2 = 2; % relaxation time factor (f2 = 2 and f1 = 1 means relaxation takes twice as long as contraction)

%% generate frames

% controls width of gaussian window (default: scale)
sigma = 2;

% orientation (in radians) of component
wOrientation = 0;
theta = -(wOrientation-90)/360*2*pi;

% controls elongation in direction perpendicular to wave
stretch = 1;
gamma = 1/(1+stretch);

% width and height of kernel
support = 75;

% wavelength (default: 4*sigma)
npeaks = 1;
lambda = 1/npeaks*4*sigma;

% phase offset (in radians)
psi = 0;

xmin = -support;
xmax = -xmin;
ymin = xmin;
ymax = xmax;

xdomain = xmin:xmax;
ydomain = ymin:ymax;

[x,y] = meshgrid(xdomain,ydomain);

xx = linspace(-3*pi,3*pi,300); % 300
dsts = 10+dMax/2*sawtooth(xx,f1*pi/4,f2*pi/4);
% contraction portion f1*pi/4 (of pi), relaxation portion f2*pi/4 (of pi)

S = zeros(size(x,1),size(x,2),length(dsts)-19);

for i = 10:length(dsts)-10 % frame
    disp(i)
    mr = zeros(size(x));
    
    for j = -40:20:40

        dst = dsts(i);
        for d = linspace(-4*dst,-dst,4)
            xprime = cos(theta)*(x-j)+sin(theta)*(y-d);
            yprime = -sin(theta)*(x-j)+cos(theta)*(y-d);
            expf = exp(-0.5/sigma^2*(xprime.^2+gamma^2*yprime.^2));
            mr = mr+expf.*cos(2*pi/lambda*xprime+psi);
        end
        
        dst = dsts(i-9);
        for d = linspace(0,4*dst,5)
            xprime = cos(theta)*(x-j)+sin(theta)*(y-d);
            yprime = -sin(theta)*(x-j)+cos(theta)*(y-d);
            expf = exp(-0.5/sigma^2*(xprime.^2+gamma^2*yprime.^2));
            mr = mr+expf.*cos(2*pi/lambda*xprime+psi);
        end

        % mr = imrotate(mr,orientation,'bicubic','crop');
        
    end
    mr = imrotate(mr,90,'bicubic','crop');
    
    mr = 0.1+0.8*normalize(mr.*(mr > 0));
    rnd = randn(size(mr));
    rnd = rnd.*(rnd > -1);
    rnd = rnd.*(rnd < 1);
    S(:,:,i-9) = mr+0.05*rnd;

%     S(:,:,i-9) = normalize(mr.*(mr > 0));
end

timeLapseViewTool(S);
figure
plot(dsts), title('distances')
dsts = dsts(10:length(dsts)-10);

% save distances
save(sprintf('~/Downloads/SarcTrackSampleVideos/Synth/Sample_%02d_%d%d.mat',dMax*10,f1,f2),'dsts');

%%

% nr = size(S,1);
% nc = size(S,2);
% S4 = zeros(2*nr,2*nc,size(S,3));
% 
% for i = 1:size(S,3)
%     I = S(:,:,i);
%     S4(1:nr,1:nc,i) = I;
%     J = imrotate(I,45,'bicubic','crop');
%     S4(nr+1:end,1:nc,i) = J;
%     J = imrotate(I,90,'bicubic','crop');
%     S4(1:nr,nc+1:end,i) = J;
%     J = imrotate(I,135,'bicubic','crop');
%     S4(nr+1:end,nc+1:end,i) = J;
%     
%     I = 0.1+0.8*S4(:,:,i);
%     rnd = randn(size(I));
%     rnd = rnd.*(rnd > -1);
%     rnd = rnd.*(rnd < 1);
%     S4(:,:,i) = I+0.05*rnd;
% end
% 
% timeLapseViewTool(S4)

%% save movie

vidObj = VideoWriter(sprintf('~/Downloads/SarcTrackSampleVideos/Synth/Sample_%02d_%d%d.avi',dMax*10,f1,f2));
open(vidObj);

for k = 1:size(S,3)
    disp(k)
    writeVideo(vidObj,uint8(255*S(:,:,k)));
end

close(vidObj);
disp('done')

%% check saved movie

v = VideoReader(sprintf('~/Downloads/SarcTrackSampleVideos/Synth/Sample_%02d_%d%d.avi',dMax*10,f1,f2));
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

T = list2stack(Frames);
timeLapseViewTool(T);