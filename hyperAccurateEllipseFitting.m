function data = hyperAccurateEllipseFitting(im)
% hyperAccurateEllipseFitting, ellipse fitting algorithm that uses Kanatani's
% hyperaccurate fit method (2010). This code is a function wrapper for that
% found http://uk.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces
% This wrapper uses a set of default settings for running the code.
%
% data = ellipsesFromTriangles(im,maxLength,resAngular,estNumber)
%
% Details   Fits ellipses by using the algorithm proposed in "Hyperaccurate
%           Ellipse Fitting without Iterations" by Kanatani, K. and
%           Rangarajan, P. in 2010.
%
%           This code is a wrapper for the codes implemented by Levente
%           Hunyadi. The required codes can be found at:
%               http://uk.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces
%           Please add these codes to your MATLAB path.
%
%           Please cite the original works and codes if you use this software
%           in your research.
%
% Inputs    im - edge image (n-by-m array) or ellipse data points (p-by-2 array)
%
% Outputs   data - a matrix of ellipses. Each row contains five elements:
%           the center of the ellipse, its major and minor axes and
%           orientation of the major axis.
%
% Examples:
% data = ellipsesFromTriangles(im,[10,20],45), fits an ellipse with half-axes
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
% See also HEDAR, ELLIPSESFROMTRIANGLES, HEDAREXPERIMENTS
%% Inputs
if size(im,2)>2
    % assume edge image
    im = logical(im); % input should be logical
    [y,x] = find(im);
else
    % assume point data
    x = im(:,1);
    y = im(:,2);
end
data = zeros(1,5);

%% Fit Ellipse
p = quad2dfit_hyperaccurate(x,y);

%% Extract Ellipse Parameters
% adjust parameters in vector p
p(2) = 0.5 * p(2);
p(4:5) = 0.5 * p(4:5);

% use implicit equation to compute ellipse semi-axes
q = 2*(p(1)*p(5)^2+p(3)*p(4)^2+p(6)*p(2)^2-2*p(2)*p(4)*p(5)-p(1)*p(3)*p(6))/(p(2)^2-p(1)*p(3));
r = realsqrt((p(1)-p(3))^2+4*p(2)^2);
data(:,3) = 2*realsqrt(q/(r-(p(1)+p(3))));   % major axis
data(:,4) = 2*realsqrt(q/(-r-(p(1)+p(3))));  % minor axis

% check
if data(:,3)<data(:,4)
    data(:,3) = data(:,4); %major axis
    data(:,4) = 2*realsqrt(q/(r-(p(1)+p(3))));%minor axis
end

% readjust parameters in vector p
p(2) = 2 * p(2);
p(4:5) = 2 * p(4:5);

% calculate translation and rotation of ellipse
data(:,1:2) = imconictranslation(p);%position in space
R = imconicrotation(imconictranslate(p, -data(:,1:2)/2));%rotation matrix
data(:,5) = acosd(-R(1,1));%orientation of major axis compared to x axis

end
