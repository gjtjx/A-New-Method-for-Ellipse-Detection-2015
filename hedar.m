function data = hedar (im,maxLength,resLength,resAngular,th)
%% hedar, Parametrisable Object Detection in 2D with Hilbert-powered step detection (developed for ellipses)
%
% data = hedar (im,maxLength,resLength,resAngular,th)
%
% Details   This function contains the parametrisable object detection
%           algorithm, developed for ellipse detection, as implemented
%           for 'A New Method for Ellipse Detection' by Carl J.
%           Nelson, Philip T. G. Jackson and Boguslaw Obara in 2015. This
%           is an advanced method that utilised the Hilbert transform
%           to find steps in the erosion signals.
%
% Inputs    im - a 2D binary or greyscale image
%           maxLength - maximum major axis length to search for (in
%           pixels; default is the length of the smallest image dimension;
%           if this argument is a vector then minLength = maxLength(1),
%           maxLength=maxLength(2), thus a minimum axis length can be set)
%           resLength - the resolution to search for axis lengths between
%           minLength and maxLength (default is 1)
%           resAngular - the resolution to search for axis orientation
%           between 0 and 180 degrees (default is 45 degrees)
%           th - the threshold (within 0 and 1) over which any step in
%           signal must be to contribute to the ellipse detection stage
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Examples:
% data = hedar (im,20,5), runs hedar on im looking for ellipses of maximum
% size 20 pixels and resolution of +/- 5 pixels and with angular
% resolution of 45 degrees (the default), considers all signals with a step
% of >=20% (default) the maximum to be contributing
% data = hedar (im,[10,20],[],5,0), runs hedar on im looking for ellipses of
% axis size between 10 and 20 pixels and angular resolution of 5 degrees,
% considers all signals with a step of >=0% the maximum to be contributing,
% i.e. all steps
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%
% See also POD2, ELLIPTICALHOUGH, PODEXPERIMENTS

%% Inputs
if nargin<5; th = 0.2; end
if nargin<4; resAngular = 1; end
if nargin<3; resLength = 1; end
if nargin<2; maxLength = min(size(im)); end
if length(maxLength)>1
    minLength = maxLength(1);
    maxLength = maxLength(2);
else
    minLength = 1;
end
%% Set Up
[m,n] = size(im);
lStepNumber = ceil((maxLength-minLength+1)/resLength);
lengthSteps = linspace(minLength,maxLength,lStepNumber);
aStepNumber = ceil(180/resAngular);
angularSteps = linspace(0,180,aStepNumber+1);
angularSteps(end) = [];
clear maxLength minLength resLength
%% Granulometric Signals
granulometricSignals = zeros(m, n, lStepNumber,aStepNumber);
for lStep = 1:length(lengthSteps)
    for aStep = 1:length(angularSteps)
        %Structuring Element
        r = floor(lengthSteps(lStep)/2); a = angularSteps(aStep);
        p1 = [-r*sind(a),r*cosd(a)];
        p2 = [r*sind(a),-r*cosd(a)];
        pmin = min([p1;p2]);
        p1 = round(p1-pmin)+1; p2 = round(p2-pmin)+1;
        pmax = max([p1;p2]);
        se = false(pmax); se(p1(1),p1(2)) = true; se(p2(1),p2(2)) = true;
        clear r a p1 p2 pmin pmax
        %Erosion
        granulometricSignals(:,:,lStep,aStep) = opterode(im,se);
    end
end
clear lStep aStep lStepNumber aStepNumber se
%% Find Drop in Signals using Hilbert
signals = reshape(granulometricSignals,[],size(granulometricSignals,3),size(granulometricSignals,4));
hilbdrop=zeros(m*n,size(granulometricSignals,4));
gap=zeros(m*n,size(granulometricSignals,4));
for ti = 1:size(signals,3)
    signal = squeeze(signals(:,:,ti))';
    h = hilbert([signal(1,:);signal;signal(end,:)]);
    h = h(2:end-1,:);
    hi = imag(h);
    [~,hmi] = max(hi);
    mask = repmat((1:size(granulometricSignals,3))',1,size(signal,2));
    mask2 = repmat(hmi,size(signal,1),1);
    mask = (mask<mask2);
    if hmi~=0
        high = sum(signal.*mask)./hmi;
    else
        high = 0;
    end
    if hmi~=size(signals,3)
        low = sum(signal.*(~mask))./-(hmi-size(signals,3));
    else
        low = 0;
    end
    gap(:,ti) = high - low;
    hilbdrop(:,ti) = squeeze(hmi)';
    clear signal ti h hi hmi mask mask2 high low
end
hilbdrop = reshape(hilbdrop,size(granulometricSignals,1),size(granulometricSignals,2),size(granulometricSignals,4));
clear signals
%% Major & Minor Axes
[majorAxis,majorOrientation] = max(hilbdrop,[],3);
X = 1:size(gap,1);
idx = sub2ind(size(gap),X',majorOrientation(:));
gap = gap(idx);
gap = reshape(gap,size(granulometricSignals,1),size(granulometricSignals,2));
minorOrientation = mod(majorOrientation-((90+resAngular)/resAngular),180/resAngular)+1;
[y,x] = ndgrid(1:m,1:n);
idr = sub2ind(size(hilbdrop),y(:),x(:),minorOrientation(:));
minorAxis = hilbdrop(idr);
minorAxis = reshape(minorAxis,size(majorAxis));
clear minorOrientation hilbdrop resAngular x y idr idx X granulometricSignals
%% Background Check
backgroundM = (majorAxis==1); backgroundm = (minorAxis==1);
backgroundG = (gap<=th);
background = backgroundM | backgroundm | backgroundG;
majorAxis = majorAxis-1; minorAxis = minorAxis-1;
majorOrientation = majorOrientation-1;
majorAxis(background) = 1; majorAxis(majorAxis==0) = 1;
majorOrientation(background) = 1; majorOrientation(majorOrientation==0) = 1;
minorAxis(background) = 1; minorAxis(minorAxis==0) = 1;
clear backgroundM backgroundm backgroundG background
%% Centroids
minorAxisMasked =  minorAxis.* (minorAxis>1);
majorAxisMasked =  majorAxis.* (majorAxis>1);
regionalPeaks = majorAxisMasked.*minorAxisMasked;
% regionalPeaks = medfilt2(regionalPeaks, [3,3]);
regionalMaxima = logical(regionalPeaks - imreconstruct(regionalPeaks-1,regionalPeaks));
centroids = regionprops(regionalMaxima,'Centroid');
centroids = reshape(cell2mat(struct2cell(centroids))',2,[])';
clear minorAxisMasked majorAxisMasked regionalPeaks regionalMaxima
%% Compile Data
data = zeros(size(centroids,1),5);
data(:,1:2) = round(centroids);
idc = sub2ind([m,n],data(:,1),data(:,2));
data(:,3) = lengthSteps(majorAxis(idc));
data(:,4) = lengthSteps(minorAxis(idc));
data(:,5) = angularSteps(majorOrientation(idc));
clear idc centroids m n lengthSteps majorAxis minorAxis angularSteps majorOrientation idc
end
