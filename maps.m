% Open Experiment 5
clear all, clear all, close all, clc
fileID = fopen('experiment5.dat','r');
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
X = dataArray{:, 1};
Y = dataArray{:, 2};
Major = dataArray{:, 3};
Minor = dataArray{:, 4};
Rotation = dataArray{:, 5};
PODn = dataArray{:, 6};
PODJaqqard = dataArray{:, 7};
PODHn = dataArray{:, 8};
PODHJaqqard = dataArray{:, 9};
Houghn = dataArray{:, 10};
HoughJaqqard = dataArray{:, 11};

podmap = zeros(64);
podhmap = zeros(64);
houghmap = zeros(64);
nummap = zeros(64);
for i = 1:length(X)
    if Major(i)==max(Major)
        r = Minor(i); a = Rotation(i)+90;
        p1 = [-r*sind(a),r*cosd(a)];
        p2 = [r*sind(a),-r*cosd(a)];
        p1 = round(p1)+33; p2 = round(p2)+33;
        se = false(64); se(p1(1),p1(2)) = true; se(p2(1),p2(2)) = true;
        [r,c] = find(se);
        for j=1:2
            podmap(r(j),c(j)) = podmap(r(j),c(j)) + PODJaqqard(i);
            podhmap(r(j),c(j)) = podhmap(r(j),c(j)) + PODHJaqqard(i);
            houghmap(r(j),c(j)) = houghmap(r(j),c(j)) + HoughJaqqard(i);
            nummap(r(j),c(j)) = nummap(r(j),c(j)) + 1;
        end
    end
end
triv = (nummap==0);
podmap = 1-(podmap./nummap); podmap(triv) = 0;
podhmap = 1-(podhmap./nummap); podhmap(triv) = 0;
houghmap = 1-(houghmap./nummap); houghmap(triv) = 0;
all = [podmap(podmap>0);podhmap(podhmap>0);houghmap(houghmap>0)];
mx = max(all);
mn = min(all);

% podmap = 1-podmap;
% podhmap = 1-podhmap;
% houghmap = 1-houghmap;
% tmp = mn;
% mn = 1-mx;
% mx = 1-tmp;
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(podmap);
colormap jet, caxis([mn,mx]), axis image off;
print('podacc','-dpng','-r0'); close
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(podhmap);
colormap jet, caxis([mn,mx]), axis image off;
print('podhacc','-dpng','-r0'); close
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(houghmap');
colormap jet, caxis([mn,mx]), axis image off;
print('houghacc','-dpng','-r0'); close
% figure, subplot(131), imagesc(podmap), caxis([mn,mx]), axis image off;
% subplot(132), imagesc(podhmap), caxis([mn,mx]), axis image off;
% subplot(133), imagesc(houghmap), caxis([mn,mx]), axis image off;
cells = sum(nummap(:)>0);
disp(sum(podmap(:))/cells)
disp(sum(podhmap(:))/cells)
disp(sum(houghmap(:))/cells)

%% Experiment 4, Map
clear all, clear all, close all, clc
fileID = fopen('experiment4_complete.dat','r');
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
distx = dataArray{:, 1};
disty = dataArray{:, 2};
PODn = dataArray{:, 3};
PODHn = dataArray{:, 5};
Houghn = dataArray{:, 7};

podmap = zeros(401);
podhmap = zeros(401);
houghmap = zeros(401);
nummap = zeros(401);
for i = 1:length(distx)
    podmap(201+distx(i),201+disty(i)) = podmap(201+distx(i),201+disty(i)) + PODn(i);
    podhmap(201+distx(i),201+disty(i)) = podhmap(201+distx(i),201+disty(i)) + PODHn(i);
    houghmap(201+distx(i),201+disty(i)) = houghmap(201+distx(i),201+disty(i)) + Houghn(i);
    nummap(201+distx(i),201+disty(i)) = nummap(201+distx(i),201+disty(i)) + 1;
end
triv = (nummap==0);
podmap = (podmap./nummap); podmap(triv) = 0;
podhmap = (podhmap./nummap); podhmap(triv) = 0;
houghmap = (houghmap./nummap); houghmap(triv) = 0;
podmap2 = podmap; podmap2(podmap==0) = 0; podmap2(podmap==1) = 1; podmap2(podmap==2) = 2; podmap2(podmap>=3) = 0; podmap2(triv) = 2;
podhmap2 = podhmap; podhmap2(podhmap==0) = 0; podmap2(podhmap==1) = 1; podhmap2(podhmap==2) = 2; podhmap2(podhmap>=3) = 0; podmap2(triv) = 2;
houghmap2 = houghmap; houghmap2(houghmap==0) = 0; houghmap2(houghmap==1) = 1; houghmap2(houghmap==2) = 1; houghmap2(houghmap>=3) = 0; podmap2(triv) = 2;
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(podmap2);
colormap gray; caxis([0,2]); axis image off;
pause(0.1)
print('podprox','-dpng','-r0'); close
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(podhmap2);
colormap gray; caxis([0,2]); axis image off;
pause(0.1)
print('podhprox','-dpng','-r0'); close
figure('units','normalized','outerposition',[0 0 1 1]), imagesc(houghmap2);
colormap gray; caxis([0,2]); axis image off;
pause(0.1)
print('houghprox','-dpng','-r0'); close