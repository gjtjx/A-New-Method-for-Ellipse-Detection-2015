fileID = fopen('experiment5c.dat','r');%N.B. assumes file is in working directory
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
Major = dataArray{:, 3}; Minor = dataArray{:, 4}; Rotation = dataArray{:, 5};
N = dataArray{:, 6}; Jaqqard = dataArray{:, 7};
% Create Maps
jmap = zeros(64); nummap = zeros(64); zmap = zeros(64);
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
            if N(i)>0
                jmap(r(j),c(j)) = jmap(r(j),c(j)) + Jaqqard(i);
                nummap(r(j),c(j)) = nummap(r(j),c(j)) + 1;
            else
                zmap(r(j),c(j)) = 1;
            end
        end
    end
end
% Normalise Pixels
triv = (nummap==0);
jmap = 1-(jmap./nummap); jmap(triv) = 0;
% Find Global Min and Max
all = [jmap(jmap>0)]
mx = max(all);
mn = min(all);
%Invert for 'Uncertainty' Map
jmap = 1-jmap;
jmap(triv | zmap) = 0;
tmp = mn;
mn = 1-mx;
mx = 1-tmp;
% Plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(121), imagesc(jmap), caxis([mn,mx]), axis image off;
subplot(122), imagesc(zmap), axis image off;
