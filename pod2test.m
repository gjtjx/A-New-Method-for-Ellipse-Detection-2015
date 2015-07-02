clear, clc, close all
N=1;
bw = zeros(201,201,N);
for i=1:N
    bws = zeros(201,201);
    major=randi(10,1)+10;
    minor=randi(round(major-5),1)+5;
    rot=randi(180,1);
    idx = BOEllipse2D(201,201,100,100,major,minor,degtorad(rot-90));
    bws(idx) = 1;
    bw(:,:,i) = imfill(bws, 'holes');
end
bw = max(bw,[],3);
% Hough
tic;
params = struct;
params.uniformWeights = true;
params.maxMajorAxis = round((2*20)+1);
edges = edge(bw,'canny');
dataED = ellipseDetection(edges,params);
toc
imapp = zeros(201,201,size(dataED,1));
for i=1:size(dataED,1)
    imapps = zeros(201,201);
    idx = BOEllipse2D(size(bw,1),size(bw,2),dataED(i,1),dataED(i,2),dataED(i,3),dataED(i,4),degtorad(dataED(i,5)));
    imapps(idx) = 1;
    imapp(:,:,i) = imfill(imapps,'holes');
end
imapp = sum(imapp,3);
imapp = (imapp-min(imapp(:)))/(max(imapp(:))-min(imapp(:)));
imout = cat(3,bw,imapp,zeros(201,201));
figure, subplot(131), imshow(bw), subplot(132), imshow(imapp), subplot(133), imshow(imout)
% POD
tic
data = pod2 (bw,params.maxMajorAxis,1,1);
toc
imapp = zeros(201,201,size(data,1));
for i=1:size(data,1)
    imapps = zeros(201,201);
    idx = BOEllipse2D(size(bw,1),size(bw,2),data(i,1),data(i,2),data(i,3)/2,data(i,4)/2,degtorad(data(i,5)-90));     
    imapps(idx) = 1;
    imapp(:,:,i) = imfill(imapps,'holes');
end
imapp = sum(imapp,3);
imapp = (imapp-min(imapp(:)))/(max(imapp(:))-min(imapp(:)));
imout = cat(3,bw,imapp,zeros(201,201));
figure, subplot(131), imshow(bw), subplot(132), imshow(imapp), subplot(133), imshow(imout)