function mr = smorlet2(d,stretch,scale,orientation)

% controls width of gaussian window (default: scale)
sigma = scale;

% orientation (in radians) of component
wOrientation = 0;
theta = -(wOrientation-90)/360*2*pi;

% controls elongation in direction perpendicular to wave
gamma = 1/(1+stretch);

% width and height of kernel
support = round((2.5*sigma/gamma+d)/2)*2+1;

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

xprime = cos(theta)*x+sin(theta)*(y-d/2);
yprime = -sin(theta)*x+cos(theta)*(y-d/2);
expf = exp(-0.5/sigma^2*(xprime.^2+gamma^2*yprime.^2));
mr1 = expf.*cos(2*pi/lambda*xprime+psi);
% mi1 = expf.*sin(2*pi/lambda*xprime+psi);

xprime = cos(theta)*x+sin(theta)*(y+d/2);
yprime = -sin(theta)*x+cos(theta)*(y+d/2);
expf = exp(-0.5/sigma^2*(xprime.^2+gamma^2*yprime.^2));
mr2 = expf.*cos(2*pi/lambda*xprime+psi);
% mi2 = expf.*sin(2*pi/lambda*xprime+psi);

mr = mr1+mr2;
% mi = mi1+mi2;

mr = imrotate(mr,orientation,'bicubic','crop');
% mi = imrotate(mi,orientation,'bicubic','crop');

% mean = 0
mr = mr-sum(sum(mr))/numel(mr);
% mi = mi-sum(sum(mi))/numel(mi);

% norm = 1
mr = mr./sqrt(sum(sum(mr.*mr)));
% mi = mi./sqrt(sum(sum(mi.*mi)));

end