function [costo]=funcDE(pts,imageB)
%function that calculates the required 
%points to form the perimeter of an ellipse starting to 5 initial points and 
%it's respective fitness value
p1=pts(1);
p2=pts(2);
p3=pts(3);
p4=pts(4);
p5=pts(5);

[border(:,1),border(:,2)]=find(imageB);
if p1<=0
    p1=1;
end
if p2<=0
    p2=1;
end
if p3<=0
    p3=1;
end
if p4<=0
    p4=1;
end
if p5<=0
    p5=1;
end
pt1=[border(p1,2),border(p1,1)];
pt2=[border(p2,2),border(p2,1)];
pt3=[border(p3,2),border(p3,1)];
pt4=[border(p4,2),border(p4,1)];
pt5=[border(p5,2),border(p5,1)];
[elipsita,~,~]=DEllipse(pt1,pt2,pt3,pt4,pt5);

costo=ellipcost(imageB,elipsita);

function cost=ellipcost(origimag,comparingvect)
H=strel('disk',2);
origimag=imdilate(origimag,H);
[pun1,pun2]=size(origimag);
[total,~]=size(comparingvect);
coincidence=0;

for ii=1:total
    if ((comparingvect(ii,2))>0 && (comparingvect(ii,1)>0) &&(comparingvect(ii,2))<=pun1 && (comparingvect(ii,1)<=pun2))
        if (origimag(comparingvect(ii,2),comparingvect(ii,1))==1)
            coincidence=coincidence+1;
        end
    end
end
if total==1
    cost=rand()*-0.001;
else
    cost=(coincidence/total)*(-1);
end
