function data = ellipsesFromTriangles(im,maxLength,resAngular)
% ellipsesFromTriangles, ellipse finding algorithm that uses point and
% tangent information. This code is a function wrapper for that found
% https://bitbucket.org/cicconet/triangles_matlab. This wrapper uses a
% set of default settings for running the code.
%
% data = ellipsesFromTriangles(im,maxLength,resAngular,estNumber)
%
% Details   Finds and fits ellipses by using the algorithm proposed in
%           "Ellipse from Triangles" by Cicconet, M., Gunsalus, K., Geiger, D.and
%           Werman, M. in 2014.
%
%           This code is a wrapper for the codes used in the aforementioned
%           paper. These required repositories can be found at:
%               https://bitbucket.org/cicconet/triangles_matlab
%               https://bitbucket.org/cicconet/pat
%           Please add these repositories to your MATLAB path.
%
%           There is a local copy of two files. The local copy of localmaxima.m
%           removes the randomisation, whilst this may increase running time
%           we found it to be much more stable. The local copy of intacc.m
%           removes the 'maximum number of objects' property.
%
%           Please cite the original works if you use this software in
%           your research:
%               Marcelo Cicconet, Davi Geiger, Kristin Gunsalus, and Michael Werman.
%               Mirror Symmetry Histograms for Capturing Geometric Properties in Images.
%               IEEE Conference on Computer Vision and Pattern Recognition. Columbus, Ohio. 2014
%
% Inputs    im - Edge image (binary or grayscale)
%           maxLength - Scale value denoting the maximum major/minor axis
%           length to consider (if a vector denotes [minimum,maximum]
%           major/minor axis length)
%           resAngular - Integer value between 1 and 180 denoting angular
%           resolution in degrees
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Based on  "Ellipse from Triangles" (Cicconet, M., Gunsalus, K., Geiger,
%           D.and Werman, M.; 2014).
%
% Examples:
% data = ellipsesFromTriangles(im,[10,20],45), detects ellipses with half-axes
% of lengths between 10 and 20 pixels at an angular resolution of 45 degrees.
%
% License   Redistribution and use in source and binary forms, with or
%           without modification, are permitted provided that the
%           following conditions are met:
%           * Redistributions of source code must retain the above copyright
%               notice, this list of conditions and the following disclaimer.
%           * Redistributions in binary form must reproduce the above
%               copyright notice, this list of conditions and the
%               following disclaimer in the documentation and/or other
%               materials provided with the distribution
%
%            THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
%            IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
%            THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%            PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
%            CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%            EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
%            PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
%            PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%            LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
%            NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%            SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% See also HEDAR, HYPERACCURATEELLIPSEFITTING, HEDAREXPERIMENTS
%% Inputs
% Set the estimated semi-minor and semi-major axis interval
if nargin<3 || isempty(resAngular); resAngular = 1; end
if nargin<2 || isempty(maxLength); maxLength = min(size(im)); end
if length(maxLength)>1
    radrange = [ceil(maxLength(1)/2) ceil(maxLength(2)/2)];
else
    radrange = [1 ceil(maxLength/2)];
end

im = mat2gray(im); % input should be double and in the range [0,1]

%% Points, Tangents, and Magnitudes
ignoredirections = 1;%treat 1-180 and 181-360 as same
stretch = 1;
scale = 1;
hopsize = 5;
halfwindowsize = 1;
magthreshold = 0.01;%set to ~0
[m,a,x,y] = coefficientslist(ignoredirections,im,180/resAngular,stretch,scale,hopsize,halfwindowsize,magthreshold);

%% Create Accumulator Array
nquadsfactor = 8;%IGNORE IF USING LOCAL intacc %set to an integer below length(m) to use randomisation
mindist = max(1,0.25*radrange(1));%minimum distance between points to consider them a useful set
maxdist = min(max(size(im)),2*radrange(2));%maximum distance between points to consider them a useful set
[A,pairs,iquads] = intacc(m,a,x,y,size(im,1),size(im,2),nquadsfactor,mindist,maxdist,radrange);

%% Locate Local Maxima Centers
hsize = 5;%Sigma for Gaussian blur
halfwindow = 8;%half width of neighbourhood to examine
mindistbetcent = 2*radrange(1);%minimum distance between maxima, set to smallest expected diameter
lowerbound = 0.25;%eps;%lower intensity bound, set to ~0
minarea = 0.1;%eps;%minimum maxima area, set to ~0
[centers,~,~] = localmaxima(A,hsize,halfwindow,lowerbound,minarea,mindistbetcent);

%% Clustering and Separating of Ellipses
proximitythreshold = 2;
ellipseindices = clusterbyellipse(pairs,iquads,centers,proximitythreshold);

%% Remaining Parameters
[~,~,data] = paintellipses(im,m,a,x,y,centers,pairs,iquads,ellipseindices);

if size(data,1)>0
    data(:,3:4) = 2*data(:,3:4);%change half-axis to axis length
    data(:,5) = rad2deg(data(:,5))-90;%set rotation from x axis not y axis
    data(:,1:2) = data(:,2:-1:1);%invert x and y
end
end
