% Copyright � 2014 New York University.
% See notice at the end of this file.

function [A,pairs,iquads] = intacc(m,a,x,y,nr,nc,nquadsfactor,mindist,maxdist,radrange)
%%Replaces the default intacc and removes randomisation
%disp('Not using randomisation!')

% intersections accumulator

% find triangles
s = length(m);
npairs = s*(s-1)/2;
pairs = zeros(npairs,7); % used, index 1, index 2, midp, dirv
pairindex = 0;
rp = randperm(s);
for jj = 1:s-1
    % fprintf('pairs of points: %f\n', jj/s);
    for kk = jj+1:s
        j = rp(jj);
        k = rp(kk);
        pairindex = pairindex+1;
        p = [x(j); y(j)];
        q = [x(k); y(k)];
        taup = [cos(a(j)); sin(a(j))];
        tauq = [cos(a(k)); sin(a(k))];
        if norm(p-q) > mindist && norm(p-q) < maxdist
            midp = round((p+q)/2); % midpoint
            [intp,denom] = intertanlin(p,q,taup,tauq);
            if denom ~= 0 % tangent lines intersect
                pn = p-intp;
                pn = pn/norm(pn);
                qn = q-intp;
                qn = qn/norm(qn);
                if abs(dot(pn,qn)) < 0.95 % triangle not too thin or fat
                    dirv = midp-intp; dirv = dirv/norm(dirv);
                    pairs(pairindex,:) = [1 j k midp' dirv'];
                end
            end
        end
    end
end

pairs = pairs(pairs(:,1) == 1,:); % throwing away unused pairs (triangles)
npairs = size(pairs,1);

% accumulate for pairs of triangles
nquads = npairs*(npairs-1)/2;
iquads = zeros(nquads,2);
A = zeros(nr,nc);
count = 0;

% non random
for i=1:npairs
    for j=1:(i-1)
            midp1 = pairs(i,4:5)';
            dirv1 = pairs(i,6:7)';
            midp2 = pairs(j,4:5)';
            dirv2 = pairs(j,6:7)';
            [intp,~] = intertanlin(midp1,midp2,dirv1,dirv2);

            indj1 = pairs(i,2);
            indk1 = pairs(i,3);
            indj2 = pairs(j,2);
            indk2 = pairs(j,3);
            p1 = [x(indj1); y(indj1)];
            q1 = [x(indk1); y(indk1)];
            p2 = [x(indj2); y(indj2)];
            q2 = [x(indk2); y(indk2)];

            if norm(intp-size(A)'/2) < min(size(A))/2 && ...
                    norm(intp-p1) >= radrange(1) && norm(intp-p1) < radrange(2) && ...
                    norm(intp-q1) >= radrange(1) && norm(intp-q1) < radrange(2) && ...
                    norm(intp-p2) >= radrange(1) && norm(intp-p2) < radrange(2) && ...
                    norm(intp-q2) >= radrange(1) && norm(intp-q2) < radrange(2)
                A(intp(1),intp(2)) = A(intp(1),intp(2))+1;
            end
            count = count+1;
            iquads(count,:) = [i j];
    end
end
A = A/max(max(A));

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
