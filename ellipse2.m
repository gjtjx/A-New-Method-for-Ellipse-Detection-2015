function bw = ellipse2(sz,pos,major,minor,phi)  
%% Inputs
if nargin<5
    phi = 0;
end
if length(sz)==2
    m = sz(1); n = sz(2);
else
    m = sz; n = sz;
end
x0 = pos(1); y0 = pos(2);
clear sz pos
%% Set-Up
bw = zeros(m,n);
%% Circumference & Parameterisation in Theta
h = (major-minor)^2/(major+minor)^2;
circ = pi*(major+minor)*(1+((3*h)/(10+sqrt(4-(3*h)))));%Approximate (Ramanujan, 1914)
theta = linspace(0,2*pi,100*round(circ));
%% Calculate Ellipse Outline Pixels
x = x0 + 0.5*(major*cos(theta)*cosd(-phi) - minor*sin(theta)*sind(-phi));
y = y0 + 0.5*(major*cos(theta)*sind(-phi) + minor*sin(theta)*cosd(-phi));
%% Boundary Issues
x(x<1) = 1; y(y<1) = 1;
x(x>n) = n; y(y>m) = m;
%% Draw Outline
idx = sub2ind(size(bw),round(y),round(x));
bw(idx) = 1;
%% Fill Ellipse
%bw = bwfill(bw,'holes');%needed for octave
bw = imfill(bw,'holes');
end