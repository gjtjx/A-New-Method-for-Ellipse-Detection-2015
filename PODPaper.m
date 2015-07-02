% Input: BW, binary image, and optional settings

%% Default Settings
angularResolution = 1; lengthResolution = 1;
maxLength = min(size(bw)); minLength = 1;

%% Set Up
numLengths = ceil((maxLength-minLength+1)/lengthResolution);
lengthSteps = linspace(minLength,maxLength,numLengths);
numAngles = ceil(180/angularResolution);
angularSteps = linspace(0,180,numAngles+1); angularSteps(end) = [];

%% Granulometric Signals
for lStep = 1:length(lengthSteps)
    for aStep = 1:length(angularSteps)
        %% Structuring Element
        r = floor(lengthSteps(lStep)/2);
        a = angularSteps(aStep);
        p1 = [-r*sind(a),r*cosd(a)];
        p2 = [r*sind(a),-r*cosd(a)];
        pmin = min([p1;p2]);
        p1 = round(p1-pmin)+1;
        p2 = round(p2-pmin)+1;
        pmax = max([p1;p2]);
        se = false(pmax); se(p1(1),p1(2)) = true; se(p2(1),p2(2)) = true;
        
        %% Erosion
        granulometricSignals(:,:,lStep,aStep) = opterode(bw,se);
    end
end

%% Find Drops in Signals
granulometricPrime = diff(granulometricSignals,1,3);
[~,primeMinima] = squeeze(min(granulometricPrime,[],3));

%% Major & Minor Axes
[majorAxis,majorOrientation] = max(primeMinima,[],3);
minorOrientation = mod(majorOrientation-((90+angularResolution)/angularResolution),180/angularResolution)+1;
[y,x] = ndgrid(1:size(bw,1),1:size(bw,2));
idr = sub2ind(size(primeMinima),y(:),x(:),minorOrientation(:));
minorAxis = primeMinima(idr); minorAxis = reshape(minorAxis,size(majorAxis));

%% Background Check
majorAxis = majorAxis-1;
minorAxis = minorAxis-1;
majorOrientation = majorOrientation-1;
majorAxis(majorAxis==0|minorAxis==0) = 1;
majorOrientation(majorAxis==0|minorAxis==0) = 1;
minorAxis(majorAxis==0|minorAxis==0) = 1;

%% Centroids
minorAxisMasked =  minorAxis.*(minorAxis>1);
majorAxisMasked =  majorAxis.*(majorAxis>1);
regionalPeaks = majorAxisMasked.*minorAxisMasked;
regionalPeaks = medfilt2(regionalPeaks, [3,3]);
regionalMaxima = logical(regionalPeaks - imreconstruct(regionalPeaks-1,regionalPeaks));
centroids = regionprops(regionalMaxima,'Centroid');

% Output: data, array of ellipse metrics
