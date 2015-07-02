function [data]=mainEllipseDetector(image)
%% Ellipse Detection by Diferential Evolution
%This code was based on the content of the research article 
% "An Improved Computer Vision Method for White Blood Cells Detection"
%Erik Cuevas, Margarita D�az, Miguel Manzanares, Daniel Zaldivar and
%Marco Perez- Cisneros, Computational and Mathematical Methods in Medicine,
%vol. 2013, Hindawi Publishing Corporation

%image >>> is a greyscale image
%outellipses >>> is a 3-dimensional matrix that stores the five points of the
%best ellipses founded

%% Set-Up
warning('off','all');
image = (image-min(image(:)))/(max(image(:))-min(image(:)));
thresmin = 0; thresmax = 0.5*graythresh(image);
win_size=5;
%win_size is a value that define how thick the edge-map will be, this
%number is typically 5
[IB]=Imaborders(image,win_size,thresmin,thresmax);  % the edge-map is obtained
VTRe=-0.75; % Value to Reach, minimum fitness value to accept an ellipse candidate
NP=30; %Total population
Itera=5000; %Total number of generations
bestfit=-1;
fitthresh = -0.5;
finito=0;
ellipfound=0; %Total number of ellipses saved and found in the image

%The method consists basically into divided the original image on sub-images
%containing smaller structures, so that the search is faster.

IB1=bwlabel(IB); %The edge-map its labeled to do the search fo every structure in the image
[S1,S2]=size(IB1);
emptyIB=zeros(S1,S2);
total=max(max(IB1));
figure, subplot(122), imshow(image);

for ii=1:total
    emptyIB=zeros(S1,S2);
    [a(:,2),a(:,1)]=find(IB1==ii);
    [taux,~]=size(a);
    for jj=1:taux
    	emptyIB(a(jj,2),a(jj,1))=1; %each label in the edge-map must be store in a single image
    end
    clear a
    [vecbordXY(:,1),vecbordXY(:,2)]=find(emptyIB);
    [searchspace,~]=size(vecbordXY); %The search space is defined by the points on the edge map for every structure
    if searchspace>=S1*S2*0.0001;
        %devec3 corresponds to the code of Differential Evolution Algorithm (Storn &
        %Price 1995) in this case we use the one provided by the Authors on
        %his web-site
        %You can exchange your own Optimizacion Algorithm in this place
        %In this case we use a version of Differential Evolution algorithm
        %By Storn % Price (1996), with the 8 strategy described in the code
        %of the respective algorithm
        [points,bestfit,~]=devec3('funcDE',VTRe,5,[1 1 1 1 1],[searchspace searchspace searchspace searchspace searchspace],emptyIB,NP,Itera,0.8,0.8,8,1);
%The five points founded to be the best so far for Differential Evolution,
%must be identified in de edge map
        points=round(points);
        p1bw=vecbordXY(points(1),:);
        p2bw=vecbordXY(points(2),:);
        p3bw=vecbordXY(points(3),:);
        p4bw=vecbordXY(points(4),:);
        p5bw=vecbordXY(points(5),:);
%This points must be switched, in the case of images the x and y axis are
%opposite
        p1=[p1bw(1,2),p1bw(1,1)];
        p2=[p2bw(1,2),p2bw(1,1)];
        p3=[p3bw(1,2),p3bw(1,1)];
        p4=[p4bw(1,2),p4bw(1,1)];
        p5=[p5bw(1,2),p5bw(1,1)];
        [SS1,SS2]=size(emptyIB);
%With the 5 points an ellipse trajectory is generated
        [ellipse,rmin,rmax,X0,Y0,theta]=DEllipse(p1,p2,p3,p4,p5);
%Only if this ellipse candidate has a fitness value less than -0.65 
%(that is 65% of coincidences) is drawn
        if bestfit>=fitthresh
            ima2=zeros(SS1,SS2,3);
            ima2(:,:,1)=emptyIB;
            [times,~]=size(ellipse);
            for iii=1:times
                    if(ellipse(iii,1)<=0)||(ellipse(iii,2)<=0)||(ellipse(iii,1)>SS2)||(ellipse(iii,2)>=SS1)
                    else
                        ima2(ellipse(iii,2),ellipse(iii,1),2)=1; %the ellipse candidate is drawn and depicted in pink color 
                        %point by point of its trajectory
                    end
            end
            ima2(:,:,3)=IB;
            subplot(121), imshow(ima2);
            use='';

%After the ellipse candidate is shown, it allows to the user decide if it keeps it or not            
            while(not(strcmp(use,'y'))&&not(strcmp(use,'Y'))&&not(strcmp(use,'n'))&&not(strcmp(use,'N')))             
                use = input('Would you like to save this ellipse? [Y]/n: ', 's');
                if(isempty(use))
                    use='Y';
                    
                end
            end
            if ((strcmp(use,'y'))||(strcmp(use,'Y')))
                use=1;
            else
                use=0;
            end
           
%If the user decides to save th ellipse candidate as part of final results, this is store in "allellipses"            
            if (use==1)
                ellipfound=ellipfound+1;
                allellipses(:,:,ellipfound)=[p1;p2;p3;p4;p5];
                data(ellipfound,:) = [X0,Y0,rmax,rmin,theta];
                for i2i=1:S1
                    for j2j=1:S2
                        if(IB1(i2i,j2j)==ii)
                            IB(i2i,j2j)=0;
                        end
                    end
                end 
            end
        else %if the ellipse candidate has a fitness value better than -0.65 (more than 65% of coincidences) then this
            %ellipse is automatically saved             
                ellipfound=ellipfound+1;
                allellipses(:,:,ellipfound)=[p1;p2;p3;p4;p5];
                data(ellipfound,:) = [X0,Y0,rmax,rmin,theta];
                for i2i=1:S1
                    for j2j=1:S2
                        if(IB1(i2i,j2j)==ii)
                            IB(i2i,j2j)=0;
                        end
                    end
                end               
        end
    end
    clear vecbordXY
end

%After the procedure of the search through the entire number of labeled
%points is done, it must be applied the search to all the reamining edge points
while (finito==0)
    [vecbordXY(:,1),vecbordXY(:,2)]=find(IB);
    [searchspace,~]=size(vecbordXY);
    if(searchspace>=S1*S2*0.0001) %the process only is applied if the total points remaining is bigger enough
        [points,bestfit,funcalls]=devec3('funcDE',VTRe,5,[1 1 1 1 1],[searchspace searchspace searchspace searchspace searchspace],emptyIB,NP,Itera,0.8,0.8,8,1);

        points=round(points);
        p1bw=vecbordXY(points(1),:);
        p2bw=vecbordXY(points(2),:);
        p3bw=vecbordXY(points(3),:);
        p4bw=vecbordXY(points(4),:);
        p5bw=vecbordXY(points(5),:);
        p1=[p1bw(1,2),p1bw(1,1)];
        p2=[p2bw(1,2),p2bw(1,1)];
        p3=[p3bw(1,2),p3bw(1,1)];
        p4=[p4bw(1,2),p4bw(1,1)];
        p5=[p5bw(1,2),p5bw(1,1)];
% % % % % % % % % % % % % % % % % % % % % % % % 
        [SS1,SS2]=size(IB);
        [ellipse,rmin,rmax,X0,Y0,theta]=DEllipse(p1,p2,p3,p4,p5);
        ima2=zeros(SS1,SS2,3);
        ima2(:,:,1)=IB;
        [times,~]=size(ellipse);
        for ii=1:times
        	if~(ellipse(ii,1)<=0)||(ellipse(ii,2)<=0)||(ellipse(ii,1)>SS2)||(ellipse(ii,2)>=SS1)
                ima2(ellipse(ii,2),ellipse(ii,1),2)=1;
            end
        end
        subplot(121), imshow(ima2);
        use='';
        while(not(strcmp(use,'y'))&&not(strcmp(use,'Y'))&&not(strcmp(use,'n'))&&not(strcmp(use,'N')))
            use = input('Would you like to save this ellipse? [Y]/n: ', 's');
            if(isempty(use))
                use='Y';
            end
        end
        if ((strcmp(use,'y'))||(strcmp(use,'Y')))
            use=1;
        else
            use=0;
        end
%if the user decides to save an ellipse candidate, a procedure to erase the coincident points in the edge map ita realized        
        while(use==1)
            %the window size refers to how thick its necessary to draw a
            %virtual ellipse over the edge map and erase the coincident
            %points, this value must be typically 15
            wind_si = input('Define window size: ', 's');
            wind_si=str2double(wind_si);
            IBtemp=zeros(SS1,SS2,3);
            IBtemp(:,:,1)=IB;
            IBtemp(:,:,2)=ellipseremove(IB,ellipse,wind_si);
            imshow(IBtemp);
            finalch='';
            while (isempty(finalch))
                %This loop gives the user the opporunity to correct the
                %window value and erase less or more points in the edge-map
                finalch = input('Is the image correct? Y/N: ', 's');
            end
            if ((strcmp(finalch,'y'))||(strcmp(finalch,'Y')))
                IB=logical(IBtemp(:,:,2));
                use=0;
                ellipfound=ellipfound+1;
                data(ellipfound,:) = [X0,Y0,rmax,rmin,theta];
                allellipses(:,:,ellipfound)=[p1;p2;p3;p4;p5];
            else
                use=1;
            end
        end    
        clear vecbordXY;
    else
        bestfit=0;
    end
    if (searchspace==0)
        disp('The search has finished')
    end
    if (bestfit>=fitthresh)
        finito=1;
    end
end
%disp('Displaying results...')
%image=displayellipse(image,allellipses);
outellipses=allellipses;



%%%%%Some of the functions needed are listed here%%%%%%%%%%
function [result]=displayellipse(inputimage,vecEllipse)
%this function draws the final ellipses saved in "allellipses" 
%over the original image and in red color
image=cat(3,inputimage,inputimage,inputimage);
[~,~,tam3]=size(vecEllipse);
[Size1,Size2]=size(inputimage);
win=3;

tempimage=zeros(Size1,Size2);
for ii=1:tam3
	pt1=vecEllipse(1,1:2,ii);
	pt2=vecEllipse(2,1:2,ii);
	pt3=vecEllipse(3,1:2,ii);
	pt4=vecEllipse(4,1:2,ii);
	pt5=vecEllipse(5,1:2,ii);
	[elliptemp,~,~,~,~]=DEllipse(pt1,pt2,pt3,pt4,pt5);
	[veces,~]=size(elliptemp);
        for jj=1:veces
            if(elliptemp(jj,1)<=0)||(elliptemp(jj,2)<=0)||(elliptemp(jj,2)>Size1)||(elliptemp(jj,1)>=Size2)
            else
                tempimage(elliptemp(jj,2),elliptemp(jj,1),1)=1;
            end
        end
end

H=strel('disk',win,0);
tempimage=imdilate(tempimage,H);
imshow(tempimage);
for ii=1:Size1
    for jj=1:Size2
        if(tempimage(ii,jj)==1)
            image(ii,jj,1)=255; %red
            image(ii,jj,2)=0; %green
            image(ii,jj,3)=0; %blue
        end 
    end
end
figure, imshow(image);
result=image;

function [imabord]=Imaborders(I,winSize,thresMin,thresMax)
%Function that returns the edge map of original image, in this, 
%the structures that were founded to be too small has been erased
if sum(I(:)==1)+sum(I(:)==0)==length(I(:))
    IM = edge(I,'canny');
else
    IM=I<thresMax;
    IM2=I>thresMin;
    IM=(IM)&(IM2); %IM is the segmented image between the thresholds thresMin and thresMax
    IM = 1-IM; %for bright objects on dark bkgd
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here starts a process to intend minimize the noise and the non-relevant structures in image
H=strel('disk',winSize);
IE=imerode(IM,H);
IF=imfill(IE,'holes');
IInv=not(IF);
IBord=(IM)&(IInv);
IBord=bwmorph(IBord,'clean');
IBord2=bwlabel(IBord);
[Ssize1,Ssize2]=size(IBord);
IBord=zeros(Ssize1,Ssize2);
total=max(max(IBord2));
for ii=1:total
    [a(:,2),a(:,1)]=find(IBord2==ii);
    [Sizeaux,~]=size(a);
    if Sizeaux>=(Ssize1*Ssize2*0.001) %ratio of 0.001 of size to erase small structures;
        %If you decrease this ratio, only the smaller structures will be
        %considerate
        for jj=1:Sizeaux
            IBord(a(jj,2),a(jj,1))=1;
        end
    end
    clear a
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imabord=IBord;

function [resultado]=ellipseremove(IBord,Vecelip,wind)
%This function erase the trajectory if an ellipse in the edge-map
[tamxI,tamyI]=size(IBord);
[tam1,~]=size(Vecelip);
%Here, "wind" represents the value of thickness for a window to erase the
%coincident points in the edge-map
H=strel('disk',wind);
newim=zeros(tamxI,tamyI);

for ii=1:tam1
     if((Vecelip(ii,2)>0)&&(Vecelip(ii,1)>0)&&((Vecelip(ii,2)<tamxI))&&((Vecelip(ii,1)<tamyI)))
        newim(Vecelip(ii,2),Vecelip(ii,1))=1;
    end
end

newim=imdilate(newim,H);
newim=imfill(newim,'holes');
newim=not(newim);
newim=(IBord)&(newim);
imshow(newim);

resultado=newim;