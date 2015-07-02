clear, clc, close all
bw = zeros(201,201);
major=randi(90,1)+10;
minor=randi(major-10,1)+10;
rot=randi(180,1);

idx = BOEllipse2D(201,201,100,100,major,minor,degtorad(rot-90));
bw(idx)=1;
bw = imfill(bw, 'holes');
[data1,granulometricSignals] = podh (bw,round(0.75*(2*major+1)),1,1);
[data2,granulometricSignals] = podh (bw,round(0.8*(2*major+1)),1,1);
[data3,granulometricSignals] = podh (bw,round(0.85*(2*major+1)),1,1);
return
bw = awgn(bw,0.5,'measured');
bw = imadjust(bw);
% [dataN,granulometricSignalsN] = pod2 (bw,100,1,1);
bw2 = imfilter(bw,fspecial('gaussian',[5 5], 1),'symmetric');
bw2 = imadjust(bw2);
tic
[dataHB,granulometricSignalsHB] = podh (bw2,100,1,1,0.2);
toc
% figure,imshow(bw),axis off
% print 'imnoise' -dpng -r300
tic
% dataED = ellipseDetection(bw2);
toc
return
%%%%
imapp = zeros(201,201);
for i=1:size(dataN,1)
    slice = zeros(201,201);
    idx = BOEllipse2D(size(bw,1),size(bw,2),dataN(i,1),dataN(i,2),dataN(i,3)/2,dataN(i,4)/2,degtorad(dataN(i,5)-90));
    slice(idx) = 1;
    slice = imfill(slice,'holes');
    imapp(logical(slice)) = i;
end
figure, imagesc(imapp)
% print 'podnoise' -dpng -r300,axis off
%%%%
return
imapp = zeros(201,201);
for i=1:size(dataHB,1)
    slice = zeros(201,201);
    idx = BOEllipse2D(size(bw,1),size(bw,2),dataHB(i,1),dataHB(i,2),dataHB(i,3)/2,dataHB(i,4)/2,degtorad(dataHB(i,5)-90));
    slice(idx) = 1;
    slice = imfill(slice,'holes');
    imapp(logical(slice)) = i;
end
figure, imagesc(imapp)
% print 'podhnoise' -dpng -r300,axis off
%%%%
imapp = zeros(201,201);
for i=1:size(dataED,1)
    slice = zeros(201,201);
    idx = BOEllipse2D(size(bw,1),size(bw,2),dataED(i,1),dataED(i,2),dataED(i,3)/2,dataED(i,4)/2,degtorad(dataED(i,5)-90));
    slice(idx) = 1;
    slice = imfill(slice,'holes');
    imapp(logical(slice)) = i;
end
figure, imagesc(imapp)
% print 'ellipseDetection' -dpng -r300,axis off