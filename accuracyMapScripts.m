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
% See also PODEXPERIMENTS

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
sz = max(distx(:)) - min(distx(:));
pos = floor(sz/2);
podmap = zeros(sz); podhmap = zeros(sz); houghmap = zeros(sz); nummap = zeros(sz);
for i = 1:length(distx)
    podmap(pos+distx(i),pos+disty(i)) = podmap(pos+distx(i),pos+disty(i)) + PODn(i);
    podhmap(pos+distx(i),pos+disty(i)) = podhmap(pos+distx(i),pos+disty(i)) + PODHn(i);
    houghmap(pos+distx(i),pos+disty(i)) = houghmap(pos+distx(i),pos+disty(i)) + Houghn(i);
    nummap(pos+distx(i),pos+disty(i)) = nummap(pos+distx(i),pos+disty(i)) + 1;
end
% Normalise Pixels
triv = (nummap==0);
podmap = (podmap./nummap); podmap(triv) = 0;
podhmap = (podhmap./nummap); podhmap(triv) = 0;
houghmap = (houghmap./nummap); houghmap(triv) = 0;
% Lavel 'Correct', 'Close', 'Not Computed' and 'Incorrect' (respectively)
podmap2 = podmap; podhmap2 = podhmap; houghmap2 = houghmap;
podmap2(podmap==2) = 2; podhmap2(podhmap==2) = 2; houghmap2(houghmap==2) = 2;%Correct, i.e. identified two
podmap2(podmap==1) = 1; podhmap2(podhmap==1) = 1; houghmap2(houghmap==1) = 1;%Identified one
podmap2(nummap==0) = 2; podhmap2(nummap==0) = 2;  houghmap2(nummap==0) = 2;%Not computed
podmap2(podmap>=3) = 0; podhmap2(podhmap>=3) = 0; houghmap2(houghmap>=3) = 0;%Incorrect

% Plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(131), imagesc(podmap2); colormap gray; caxis([0,2]); axis image off;
subplot(132), imagesc(podhmap2); colormap gray; caxis([0,2]); axis image off;
subplot(133), imagesc(houghmap2); colormap gray; caxis([0,2]); axis image off;