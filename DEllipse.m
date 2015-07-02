%Function that draws an ellipse based on 5 points

function [salida,rmin,rmax,X0,Y0,theta]=DEllipse(p1,p2,p3,p4,p5)
%Returns a series of points that describes an ellipse trayectory through 
%the 5 points in analogical way

%The mathematical model of an ellipse is described in the paper
x=[p1(:,1);p2(:,1);p3(:,1);p4(:,1);p5(:,1)];
y=[p1(:,2);p2(:,2);p3(:,2);p4(:,2);p5(:,2)];


M = [2*x.*y y.^2 2*x 2*y ones(size(x))];

e = M\(-x.^2);

a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);

delta = b^2-a*c;

X0=round((c*d - b*f)/delta);
Y0=round((a*f - b*d)/delta);
phi = 0.5 * acot((c-a)/(2*b));

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = round(sqrt(nom/(delta* ( (c-a)*s -(c+a)))));

b_prime = round(sqrt(nom/(delta* ( (a-c)*s -(c+a)))));


rmax = round(max(a_prime, b_prime));
rmin = round(min(a_prime, b_prime));

if (a_prime < b_prime)
    phi = pi/2 - phi;
end

theta=(180*phi)/pi;
 if (isreal(rmin))&&(isreal(rmax))
    
    if (rmin>=35)&&(rmax>=35)&&(rmin<=75)&&(rmax<=75)
        %here rmin and rmax represent the possible rank for the semi-minor and 
        %semi-major axis of an ellipse
        [vecimelip,~]=ellipsimple(X0,Y0,rmax,rmin,theta);
    else
        vecimelip=[0,0];
    end
else
    vecimelip=[0,0];
end
salida=vecimelip;

function [point,s] = ellipsimple(x0, y0, rx, ry , degrees)
%function that calculates perimeter points of an ellipse, 
%through duplicate symmetric quadrants
x = 0;
y = ry;
pk = (ry^2)-(rx^2*ry)+(0.25*(rx^2));

count = 1;

points(count, 1) = x;
points(count, 2) = y;

count = count + 1;
%%%%%%%%%%%%%%%%%%% REGION 1 %%%%%%%%%%%%%%%%
while 2*ry^2*x <= 2*rx^2*y
    
    if pk < 0
        x = x + 1;
        points(count, 1) = x;
        points(count, 2) = y;
        pk = pk + (2*ry^2*x) + ry^2;
    else
        x = x + 1;
        y = y - 1;
        points(count, 1) = x;
        points(count, 2) = y;
        pk = pk + (2*ry^2*x) - (2*rx^2*y) + ry^2;
    end
    count = count + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pk = (ry^2*(x+0.5)^2)+(rx^2*(y-1)^2)-ry^2*rx^2;

%%%%%%%%%%%%%%%%%%% REGION 2 %%%%%%%%%%%%%%%%
while y > 0
    if pk >= 0
        y = y - 1;
        points(count, 1) = x;
        points(count, 2) = y;
        pk = pk - (2*rx^2*y) + rx^2;
    else
        x = x + 1;
        y = y - 1;
        points(count, 1) = x;
        points(count, 2) = y;
        pk = pk + (2*ry^2*x) - (2*rx^2*y) + rx^2;
    end
    count = count + 1;
end

pts = length(points);

for i = 1:pts
    points(i+pts, 1) = points(i,1) * -1;
    points(i+pts, 2) = points(i,2);
end

for i = 1:pts
    points(i+(pts*2), 1) = points(i,1) * -1;
    points(i+(pts*2), 2) = points(i,2) * -1;
end

for i = 1:pts
    points(i+(pts*3), 1) = points(i,1);
    points(i+(pts*3), 2) = points(i,2) * -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=points;
k=rotate(s,degrees);

k1(:,1)=[k(:,1)];
k1(:,2)=[k(:,2)];
[SS,~]=size(k1);
point=zeros(SS,2);
point(:,1) = k1(:,1) + x0;
point(:,2) = k1(:,2) + y0;

function [kl]=rotate(k,angle) 
%Function that rotates an ellipse the number of degrees specified in angle
kl(:,1)=round(k(:,1)*cos(angle) - k(:,2)*sin(angle));
kl(:,2)=round(k(:,1)*sin(angle) + k(:,2)*cos(angle));
