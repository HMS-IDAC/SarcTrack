function [M3,A3,X3,Y3,cKernels,cIKernels] = coefficientsmatrix23(I,wavDists,varargin)

% list of wavelet coefficients (magnitudes, angles, locations)
% ----------
% INPUTS
% expects double, range [0,1], grayscale image input I
% wavDist: distance between parallel wavelets
% varargin:
%   'NumAngles', default 32: number of angles considered in [0,pi) or
%       [0,2*pi) (depending on 'ConvType')
%   'WavStretch', default 0: stretch of wavelets
%   'WavScale', default 1: scale of wavelets
%   'HopSize', default 3: one coefficient is output per 'HopSize' pixels
%   'HalfWindowSize', default 1: half size of window around which best
%       coefficient is looked for; should be no larger than (HopSize-1)/2
%   'MagThreshold', default 0.01: threshold coefficient should pass to be
%       considered
%   'ConvType', default 'complex': type of wavelet used in convolution
%       'real' -- good for ridge-type edges
%       'real+' -- better than 'real' in some cases
% ----------
% OUTPUTS
% m: magnitudes, in range [0,1]
% a: angles, in radians
% x: row location of coefficient (Yeah I know, this is different from Matlab's
%   coordinate system. Sorry.)
% y: column location of coefficients

ip = inputParser;
ip.addParameter('NumAngles',32);
ip.addParameter('WavStretch',0);
ip.addParameter('WavScale',1);
ip.addParameter('HopSize',3);
ip.addParameter('HalfWindowSize',1);
ip.addParameter('MagThreshold',0.01);
ip.addParameter('ConvType','real');

ip.parse(varargin{:});
p = ip.Results;

cM = cell(1,length(wavDists));
cA = cell(1,length(wavDists));
cX = cell(1,length(wavDists));
cY = cell(1,length(wavDists));
cKernels = cell(1,length(wavDists));
cIKernels = cell(1,length(wavDists));
% tic
for i = 1:length(wavDists)
    wavDist = wavDists(i);
    [M,A,X,Y,kernels,iKernels] = coefficientsmatrix2(I,wavDist,p.NumAngles,p.WavStretch,p.WavScale,p.HopSize,p.HalfWindowSize,p.MagThreshold,p.ConvType);
    cM{i} = M;
    cA{i} = A;
    cX{i} = X;
    cY{i} = Y;
    cKernels{i} = kernels;
    cIKernels{i} = iKernels;
end
% toc

M3 = zeros(size(M,1),size(M,2),length(wavDists));
A3 = zeros(size(M,1),size(M,2),length(wavDists));
X3 = zeros(size(M,1),size(M,2),length(wavDists));
Y3 = zeros(size(M,1),size(M,2),length(wavDists));
for i = 1:length(wavDists)
    M3(:,:,i) = cM{i};
    A3(:,:,i) = cA{i};
    X3(:,:,i) = cX{i};
    Y3(:,:,i) = cY{i};
end

end
