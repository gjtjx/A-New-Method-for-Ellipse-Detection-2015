% Copyright � 2014 New York University.
% See notice at the end of this file.

function [locs,I,J] = localmaxima(IIn,hsize,halfwindow,lowerbound,minarea,mindistbetcent)
% This is a local copy of the localmaxima.m file found in the repository
% 'triangles_matlab'; this version ignores the maximum estimate property.
H = fspecial('gaussian',hsize,hsize/4);
I = conv2(double(IIn),H,'same');
I = I/max(max(I));

[nr,nc] = size(I);
J = zeros(nr,nc);
d = halfwindow; % default: 10
n = 8;
for i = d+1:nr-d
    for j = d+1:nc-d
        values1 = zeros(1,n);
        values2 = zeros(1,n);
        values3 = zeros(1,n);
        for k = 1:n
            ag = (k-1)/n*2*pi;
            vector = round(d*[cos(ag) sin(ag)]);
            values3(k) = I(i+vector(1),j+vector(2));
            vector = round(d/4*[cos(ag) sin(ag)]);
            values2(i) = I(i+vector(1),j+vector(2));
            vector = round(d/8*[cos(ag) sin(ag)]);
            values1(i) = I(i+vector(1),j+vector(2));
        end
        center = I(i,j);
        if center > max([values1 values2 values3]) && center > lowerbound
            J(i,j) = 1;
        end
    end
end

cc = bwconncomp(J,8);
stats = regionprops(cc,'Area','Centroid');
l = length(stats);
areas = zeros(1,l);
centers = zeros(2,l);
for i = 1:l
    areas(i) = stats(i).Area;
    centers(:,i) = round([stats(i).Centroid(2) stats(i).Centroid(1)]);
end

% get rid of blobs that are too small
mx = max(areas);
[mn,imn] = min(areas);
while mn < minarea*mx % default: 0.25
    areas(imn) = [];
    centers(:,imn) = [];
    [mn,imn] = min(areas);
end

% get rid of blobs for which there's a strong blob nearby
D = distbetcent(centers);
[mn,imn] = min(D);
[mn1,imn1] = min(mn);
while mn1 < mindistbetcent
    c = imn1;
    r = imn(c);
    if areas(c) < areas(r)
        idx = c;
    else
        idx = r;
    end
    areas(idx) = [];
    centers(:,idx) = [];

    D = distbetcent(centers);
    [mn,imn] = min(D);
    [mn1,imn1] = min(mn);
end

locs = centers;

end

function D = distbetcent(centers)
    d = size(centers,2);
    D = zeros(d,d);
    for i = 1:d
        for j = 1:d
            if i == j
                D(i,j) = Inf;
            else
                D(i,j) = norm(centers(:,i)-centers(:,j));
            end
        end
    end
end

% Copyright � 2014 New York University.
%
% All Rights Reserved. A license to use and copy this software and its documentation
% solely for your internal research and evaluation purposes, without fee and without a signed licensing agreement,
% is hereby granted upon your download of the software, through which you agree to the following:
% 1) the above copyright notice, this paragraph and the following paragraphs
% will prominently appear in all internal copies;
% 2) no rights to sublicense or further distribute this software are granted;
% 3) no rights to modify this software are granted; and
% 4) no rights to assign this license are granted.
% Please Contact The Office of Industrial Liaison, New York University, One Park Avenue, 6th Floor,
% New York, NY 10016 (212) 263-8178, for commercial licensing opportunities,
% or for further distribution, modification or license rights.
%
% Created by Marcelo Cicconet.
%
% IN NO EVENT SHALL NYU, OR ITS EMPLOYEES, OFFICERS, AGENTS OR TRUSTEES (?COLLECTIVELY ?NYU PARTIES?)
% BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY KIND ,
% INCLUDING LOST PROFITS, ARISING OUT OF ANY CLAIM RESULTING FROM YOUR USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
% EVEN IF ANY OF NYU PARTIES HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH CLAIM OR DAMAGE.
%
% NYU SPECIFICALLY DISCLAIMS ANY WARRANTIES OF ANY KIND REGARDING THE SOFTWARE,
% INCLUDING, BUT NOT LIMITED TO, NON-INFRINGEMENT, THE IMPLIED WARRANTIES OF  MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE, OR THE ACCURACY OR USEFULNESS,
% OR COMPLETENESS OF THE SOFTWARE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
% IF ANY, PROVIDED HEREUNDER IS PROVIDED COMPLETELY "AS IS".
% NYU HAS NO OBLIGATION TO PROVIDE FURTHER DOCUMENTATION, MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
%
% Please cite the following reference if you use this software in your research:
%
% Marcelo Cicconet, Davi Geiger, Kristin Gunsalus, and Michael Werman.
% Mirror Symmetry Histograms for Capturing Geometric Properties in Images.
% IEEE Conference on Computer Vision and Pattern Recognition. Columbus, Ohio. 2014.
