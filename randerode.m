function dst = randerode(src,kernel,R)
%% randerode carries out pseudo-morphological erosions on random pixels within an image
%
% B = randerode(A,se,R)
%
% Details   This script runs a morphological erosion-like process on a
%           random select of pixels in im using the structural element se.
% Inputs    A - 2D image (grayscale or BW)
%           se - structural element; must be flat
%           R - ratio of pixels to use, e.g. 0.5 means use 50% of pixels
% Outputs   B - eroded 2D image (grayscale or BW)
%
% Examples:
% A = rand(20);%Starting image - random
% se = strel('disk',4,0);%Structural element - disk (N.B. flat)
% B = opterode(A,se,0.5), returns grayscale image B same size as A
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%% Pad with Zeros (c.f. MATLAB's imerode default of Inf)
[sy,sx] = size(src);
[sky,skx] = size(kernel);

B = zeros(sy+(2*sky),sx+(2*skx));
B(sky+(1:sy),skx+(1:sx)) = src;

dst = B;

%% Select Random Pixels
id = randperm(numel(src),numel(src)*R);
[idy,idx] = ind2sub([sy,sx],id);

idx = idx+skx;
idy = idy+sky;

%% Strel Knowledge
[idc,idr] = meshgrid(1:skx,1:sky);
idr = idr(:);
idc = idc(:);

%% 'Erode' at Chosen Pixels
% Maxes out memory!
for kx = 1:length(idr)
	shift = [floor(idr(kx)-sky/2),floor(idc(kx)-skx/2)];
    currentMin = dst(idy,idx);
    kernelValue = kernel(idr(kx),idc(kx));
    dst(idy,idx) = min(B(idy+shift(1),idx+shift(2))-kernelValue,currentMin);
end

%% Remove Padding
dst = dst(sky+(1:sy),skx+(1:sx));
end
