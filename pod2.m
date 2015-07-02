function granulometricSignals = pod2 (bw,maxLength,resLength,resAngular)
%% Parametrised Object Detection in 2D - Developed for Ellipses
%% Inputs
if nargin<4; resAngular = 1; end
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
return
clear lStep aStep lStepNumber aStepNumber se
%% Find Drops in Signals
granulometricPrime = diff(granulometricSignals,1,3);
[~,primeMinima] = min(granulometricPrime,[],3);
primeMinima = squeeze(primeMinima);
clear granulometricPrime granulometricSignals
%% Major & Minor Axes
[majorAxis,majorOrientation] = max(primeMinima,[],3);
minorOrientation = mod(majorOrientation-((90+resAngular)/resAngular),180/resAngular)+1+1;
minorOrientation(minorOrientation==181) = 1;%TODO: fix in line above
[y,x] = ndgrid(1:size(bw,1),1:size(bw,2));
idr = sub2ind(size(primeMinima),y(:),x(:),minorOrientation(:));
minorAxis = primeMinima(idr);
minorAxis = reshape(minorAxis,size(majorAxis));
clear minorOrientation primeMinima resAngular x y idr
%% Background Check
backgroundM = (majorAxis==1); backgroundm = (minorAxis==1);
background = backgroundM | backgroundm;
majorAxis = majorAxis-1; minorAxis = minorAxis-1;
majorOrientation = majorOrientation -1;%Can this be justified?
majorAxis(background) = 1; majorAxis(majorAxis==0) = 1;
majorOrientation(background) = 1; majorOrientation(majorOrientation==0) = 1;
minorAxis(background) = 1; minorAxis(minorAxis==0) = 1;
clear backgroundM backgroundm background
%% Centroids
minorAxisMasked =  minorAxis.* (minorAxis>1);
majorAxisMasked =  majorAxis.* (majorAxis>1);
regionalPeaks = majorAxisMasked.*minorAxisMasked;
%regionalPeaks = medfilt2(regionalPeaks, [3,3]);
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
% %% Check No Detected Ellipse is Within Another
% con = [];
% if size(data,1)>1
%     for i=1:size(data,1)
%         bwi = ellipse2(size(bw),data(i,1:2),data(i,3),data(i,4),data(i,5));
% %         % Remove Ellipses that are connected to the image boundary
% %         if max(imclearborder(bwi,8))==0
% %             con = [con;i];
% %         end
%         for j=i+1:size(data,1)
%             distSq = (data(i,1)-data(j,1))^2 + (data(i,2)-data(j,2))^2;
%             majSq = ((data(i,3)+data(i,4))^2)/4;
%             if distSq<majSq
%                 bwj = ellipse2(size(bw),data(j,1:2),data(j,3),data(j,4),data(j,5));
%                 score = sum(bwi(:) & bwj(:))/ min(sum(bwi(:)>0),sum(bwj(:)>0));
%                 if score>0.8
%                     con = [con;i];
%                     break
%                 end
%             end
%         end
%     end
%     data(con,:) = [];
%     clear i j bwi bwj distSq majSq con
% end
end