function data = pod2 (bw,maxLength,resLength,resAngular)
%% pod2, Parametrisable Object Detection in 2D (developed for ellipses)
%
% data = pod2 (bw,maxLength,resLength,resAngular)
%
% Details   This function contains the parametrisable object detection
%           algorithm, developed for ellipse detection, as implemented 
%           for 'A New Method for Ellipse Detection' by Carl J.
%           Nelson, Philip T. G. Jackson and Boguslaw Obara in 2015. This
%           is the simple, finite difference method described and is
%           designed for low noise, binary images.
% Inputs    bw - a 2D, low noise, binary image
%           maxLength - maximum major axis length to search for (in
%           pixels; default is the length of the smallest image dimension;
%           if this argument is a vector then minLength = maxLength(1), 
%           maxLength=maxLength(2), thus a minimum axis length can be set)
%           resLength - the resolution to search for axis lengths between
%           minLength and maxLength (default is 1)
%           resAngular - the resolution to search for axis orientation
%           between 0 and 180 degrees (default is 45 degrees)
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Examples:
% data = pod2 (bw,20,5), runs pod2 on bw looking for ellipses of maximum
% size 20 pixels and resolution of +/- 5 pixels and with angular
% resolution of 45 degrees (the default)
% data = pod2 (bw,[10,20],[],5), runs pod2 on bw looking for ellipses of
% axis size between 10 and 20 pixels and angular resolution of 5 degrees
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
% 
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%
% See also PODH, ELLIPTICALHOUGH, PODEXPERIMENTS

%% Inputs
if nargin<4; resAngular = 15; end
if nargin<3; resLength = 1; end
if nargin<2; maxLength = min(size(bw)); end
if length(maxLength)>1
    minLength = maxLength(1);
    maxLength = maxLength(2);
else
    minLength = 1;
end

%% Set Up
[m,n] = size(bw);
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
        granulometricSignals(:,:,lStep,aStep) = opterode(bw,se);
    end
end
clear lStep aStep lStepNumber aStepNumber se
%% Find Drops in Signals
granulometricPrime = diff(granulometricSignals,1,3);
[~,primeMinima] = min(granulometricPrime,[],3);
primeMinima = squeeze(primeMinima);
clear granulometricPrime granulometricSignals
%% Major & Minor Axes
[majorAxis,majorOrientation] = max(primeMinima,[],3);
minorOrientation = mod(majorOrientation-((90+resAngular)/resAngular),180/resAngular)+1;
[y,x] = ndgrid(1:size(bw,1),1:size(bw,2));
idr = sub2ind(size(primeMinima),y(:),x(:),minorOrientation(:));
minorAxis = primeMinima(idr);
minorAxis = reshape(minorAxis,size(majorAxis));
clear minorOrientation primeMinima resAngular x y idr
%% Background Check
backgroundM = (majorAxis==1); backgroundm = (minorAxis==1);
background = backgroundM | backgroundm;
majorAxis = majorAxis-1; minorAxis = minorAxis-1;
majorAxis(background) = 1; majorAxis(majorAxis==0) = 1;
majorOrientation(background) = 1; majorOrientation(majorOrientation==0) = 1;
minorAxis(background) = 1; minorAxis(minorAxis==0) = 1;
clear backgroundM backgroundm background
%% Centroids
minorAxisMasked =  minorAxis.* (minorAxis>1);
majorAxisMasked =  majorAxis.* (majorAxis>1);
regionalPeaks = majorAxisMasked.*minorAxisMasked;
regionalPeaks = medfilt2(regionalPeaks, [3,3]);
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
