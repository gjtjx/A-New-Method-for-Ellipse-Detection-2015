%% Scripts for producing the data maps for 'A New Method for Ellipse Detection' by Nelson et al. (2015)
%
% Please run each section seperately.
%
% Details   This script contains two sections, each one can be run
% idependently and will produce the data maps as shown in Figures 7 and 8
% for 'A New Method for Ellipse Detection' by Nelson et al. (2015).
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%
% See also HEDAREXPERIMENTS

%% Create Uncertainty Map from Experiment 5:- Accuracy over m, M and theta
% Load Data
fileID = fopen('experiment5.dat','r');%N.B. assumes file is in working directory
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
Major = dataArray{:, 3}; Minor = dataArray{:, 4}; Rotation = dataArray{:, 5};
PODJaqqard = dataArray{:, 7}; PODHJaqqard = dataArray{:, 9}; HoughJaqqard = dataArray{:, 11};
% Create Maps
podmap = zeros(64); podhmap = zeros(64); houghmap = zeros(64); nummap = zeros(64);
for i = 1:length(Major)
    if Major(i)==max(Major)%Creates the maps for the largest Major axis length in the data
        r = Minor(i); a = Rotation(i)+90;
        %Identify pixels to assign to
        p1 = [-r*sind(a),r*cosd(a)]; p2 = [r*sind(a),-r*cosd(a)];
        p1 = round(p1)+33; p2 = round(p2)+33;
        se = false(64); se(p1(1),p1(2)) = true; se(p2(1),p2(2)) = true;
        [r,c] = find(se);
        %Assign data to all maps
        for j=1:2
            podmap(r(j),c(j)) = podmap(r(j),c(j)) + PODJaqqard(i);
            podhmap(r(j),c(j)) = podhmap(r(j),c(j)) + PODHJaqqard(i);
            houghmap(r(j),c(j)) = houghmap(r(j),c(j)) + HoughJaqqard(i);
            nummap(r(j),c(j)) = nummap(r(j),c(j)) + 1;
        end
    end
end
% Normalise Pixels
triv = (nummap==0);
podmap = 1-(podmap./nummap); podmap(triv) = 0;
podhmap = 1-(podhmap./nummap); podhmap(triv) = 0;
houghmap = 1-(houghmap./nummap); houghmap(triv) = 0;
% Find Global Min and Max
all = [podmap(podmap>0);podhmap(podhmap>0);houghmap(houghmap>0)];
mx = max(all);
mn = min(all);
%Invert for 'Uncertainty' Map
podmap = 1-podmap;
podhmap = 1-podhmap;
houghmap = 1-houghmap;
tmp = mn;
mn = 1-mx;
mx = 1-tmp;
% Plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(131), imagesc(podmap), caxis([mn,mx]), axis image off;
subplot(132), imagesc(podhmap), caxis([mn,mx]), axis image off;
subplot(133), imagesc(houghmap), caxis([mn,mx]), axis image off;

%% Create Detection Number Map from Experiment 4:- Clustering and Overlapping Ellipses
% Load Data
fileID = fopen('experiment4.dat','r');%N.B. assumes file is in working directory
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
distx = dataArray{:, 1}; disty = dataArray{:, 2};
PODn = dataArray{:, 3}; PODHn = dataArray{:, 5}; Houghn = dataArray{:, 7};
% Create Maps
sz = 211;%max(distx(:)) - min(distx(:));
pos = 106;
podmap = zeros(sz); podhmap = zeros(sz); houghmap = zeros(sz); nummap = zeros(sz);
distances = []; podmean = []; podhmean = []; houghmean = []; nummean = [];
for i = 1:size(distx,1)
    podmap(pos+distx(i),pos+disty(i)) = podmap(pos+distx(i),pos+disty(i)) + PODn(i);
    podhmap(pos+distx(i),pos+disty(i)) = podhmap(pos+distx(i),pos+disty(i)) + PODHn(i);
    houghmap(pos+distx(i),pos+disty(i)) = houghmap(pos+distx(i),pos+disty(i)) + Houghn(i);
    nummap(pos+distx(i),pos+disty(i)) = nummap(pos+distx(i),pos+disty(i)) + 1;
    dst = round(sqrt((distx(i))^2 + (disty(i))^2));
    idm = find(distances==dst);
    if isempty(idm)
        distances(end+1) = dst;
        podmean(end+1) = PODn(i);
        podhmean(end+1) = PODHn(i);
        houghmean(end+1) = Houghn(i);
        nummean(end+1) = 1;
    else
        podmean(idm) = podmean(idm)+PODn(i);
        podhmean(idm) = podhmean(idm)+PODHn(i);
        houghmean(idm) = houghmean(idm)+Houghn(i);
        nummean(idm) = nummean(idm) + 1;
    end
end
% Normalise Pixels
triv = (nummap==0);
podmap = (podmap./nummap); podmap(triv) = 0;
podhmap = (podhmap./nummap); podhmap(triv) = 0;
houghmap = (houghmap./nummap); houghmap(triv) = 0;
% Normalise and sort lists
triv = (nummean==0); [distances,ids] = sort(distances);
podmean = (podmean./nummean); podmean(triv) = 0; podmean = podmean(ids);
podhmean = (podhmean./nummean); podhmean(triv) = 0; podhmean = podhmean(ids);
houghmean = (houghmean./nummean); houghmean(triv) = 0; houghmean = houghmean(ids);

% Plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(231), imagesc(abs(podmap-2)); colormap parula; axis image off; caxis([0,10]);
subplot(232), imagesc(abs(podhmap-2)); colormap parula; axis image off; caxis([0,10]);
subplot(233), imagesc(abs(houghmap-2)); colormap parula; axis image off; caxis([0,10]);
subplot(234), plot(distances,podmean-2); xlim([0,sqrt(2*105^2)]); ylim([-2,10]);
subplot(235), plot(distances,podhmean-2); xlim([0,sqrt(2*105^2)]); ylim([-2,10])
subplot(236), plot(distances,houghmean-2); xlim([0,sqrt(2*105^2)]); ylim([-2,10])
