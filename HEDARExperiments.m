function HEDARExperiments(en,N)
%% Experiments for 'A New Method for Ellipse Detection' by Nelson et al. (2015)
%
% HEDARExperiments(en,N)
%
% Details   This function contains all the experiments implemented for
%           the paper 'A New Method for Ellipse Detection' by Carl J.
%           Nelson, Philip T. G. Jackson and Boguslaw Obara in 2015.
% Inputs    en - Experiment number, between 1 and 6
%           N - number or repeat trials (only reqiured for certain
%           experiments; default is 1)
% Outputs   N/A - all experiments save their data to file
%
% Examples:
% HEDARExperiments(4), runs the experiments used to create Figure 7
% HEDARExperiments(2,10), runs the experiments used to create Figure 6b ten times
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
%
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%
% See also HEDAR, ELLIPTICALHOUGH, ELLIPSESFROMTRIANGLES, ACCURACYMAPSCRIPTS

switch en
%% Experiment 1 - Change of image size, set kernel size range (binary image)
    case 1
        % Inputs
        if nargin==1; disp('Running only once'); N = 1; end;
        % Create File for Results
        headers = {'N','m',...
            'HEDAR:time','HEDAR:Jaccard',...
            'Hough:time','Hough:Jaccard',...
            'EFT:time','EFT:Jaccard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment1.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        % Set-Up
        major=randi(23,1)+7; minor=randi(major-7,1)+7;
        rot=randi(180,1); ms = [64,128,256,512];
        for m=1:length(ms)
            % Data
            data = cell(N,8);
            data(:,2) = cellstr(num2str(repmat(ms(m),N,1)));
            % Create Image
            bw = ellipse2(ms(m),[ceil((ms(m)+1)/2),ceil((ms(m)+1)/2)],major,minor,rot);
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run HEDAR
                tic, results = hedar(bw,32);
                data{rn,3} = num2str(toc);
                % HEDAR Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(results,1));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                tic, edges = edge(bw,'canny');
                results = ellipticalHough(edges,32);
                data{rn,5} = num2str(toc);
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Ellipses From Triangles (EFT)
                tic, edges = imgradient(bw);
                results = ellipsesFromTriangles(edges,32);
                data{rn,7} = num2str(toc);
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
            end
            % Append Data to File
            fid = fopen('experiment1.dat','a+');
            for o=1:N
                dato = strjoin(data(o,:),',');
                fprintf(fid,'%s\r\n',dato);
            end
            fclose(fid); clear data dato fid bw
        end
%% Experiment 2 - Change of kernel size (max), set image size (binary image)
    case 2
        % Inputs
        if nargin==1; disp('Running only once'); N = 1; end;
        % Create File for Results
        headers = {'N','k',...
            'HEDAR:time','HEDAR:Jaccard',...
            'Hough:time','Hough:Jaccard',...
            'EFT:time','EFT:Jaccard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment2.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        % Set-Up
        major=randi(7,1)+7; minor=randi(major-7,1)+7;
        rot=randi(180,1);
        ks = 15:5:65;
        % Create Image
        bw = ellipse2(256,[ceil((256+1)/2),ceil((256+1)/2)],major,minor,rot);
        for k=1:length(ks)
            % Data
            data = cell(N,8);
            data(:,2) = cellstr(num2str(repmat(ks(k),N,1)));
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run HEDAR
                tic, results = hedar (bw,ks(k));
                data{rn,3} = num2str(toc);
                % HEDAR Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(results,1));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                tic, edges = edge(bw,'canny');
                results = ellipticalHough(edges,ks(k));
                data{rn,5} = num2str(toc);
                % Hough Jaccard
                bwo = zeros(size(bw));bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Ellipses From Triangles (EFT)
                tic, edges = imgradient(bw);
                results = ellipsesFromTriangles(edges,ks(k));
                data{rn,7} = num2str(toc);
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
            end
            % Append Data to File
            fid = fopen('experiment2.dat','a+');
            for o=1:N
                dato = strjoin(data(o,:),',');
                fprintf(fid,'%s\r\n',dato);
            end
            fclose(fid); clear data dato fid
        end
%% Experiment 3 - Different numbers of objects (non-overlapping)
    case 3
        % Inputs
        if nargin==1; disp('Running only once'); N = 1; end;
        % Create File for Results
        headers = {'N','number',...
            'HEDAR:n','HEDAR:time','HEDAR:Jaccard',...
            'Hough:n','Hough:time','Hough:Jaccard',...
            'EFT:n','EFT:time','EFT:Jaccard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment3.dat','w');
        fprintf(fid,'%s\r\n',headers);
        fclose(fid);
        %% Set-up Figure
        % Create Figure
        bw = zeros(256,256);
        xposs = cell(200,1); yposs = cell(200,1); majors = cell(200,1);
        for l=1:200;
            % Add Ellipse to Figure
            major=randi(7,1)+7; minor=randi(major-7,1)+7; rot=randi(180,1);
            xpos = randi(226)+15; ypos = randi(226)+15;
            if l>1
                overlapping = true;
                while overlapping
                    overlapping = false;
                    for ps=1:l
                        distmap = zeros(256,256);
                        distmap(xposs{ps},yposs{ps})=1;
                        distmap = bwdist(distmap);
                        if distmap(xpos,ypos)<=(major+majors{ps})
                            xpos = randi(226)+15; ypos = randi(226)+15;
                            overlapping = true;
                            break
                        end
                    end
                end
            end
            clear overlapping ps distmap
            % disp([xpos,ypos,major,minor,rot])
            bwadd = ellipse2(256,[xpos,ypos],major,minor,rot);
            xposs{l} = xpos; yposs{l} = ypos; majors{l} = major;
            bw = double((bw+bwadd)>0);
            % figure(1), subplot(221), imshow(bw)
            clear major minor rot xpos ypos overlapping distmap
            % Data
            data = cell(N,11);
            data(:,2) = cellstr(num2str(repmat(l,N,1)));
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run HEDAR
                tic, results = hedar (bw,32);
                data{rn,4} = num2str(toc);
                data{rn,3} = num2str(size(results,1));
                % disp(results)
                % HEDAR Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(results,1));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                % figure(1), subplot(222), imshow(bwo)
                data{rn,5} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                tic, edges = edge(bw,'canny');
                results = ellipticalHough(edges,15);
                data{rn,7} = num2str(toc);
                data{rn,6} = num2str(size(results,1));
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                % figure(1), subplot(223), imshow(bwo)
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Ellipses From Triangles (EFT)
                tic, edges = imgradient(bw);
                results = ellipsesFromTriangles(edges,15);
                data{rn,10} = num2str(toc);
                data{rn,9} = num2str(size(results,1));
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3); clear o results
                % figure(1), subplot(224), imshow(bwo), drawnow
                data{rn,11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
            end
            % Append Data to File
            fid = fopen('experiment3.dat','a+');
            for o=1:N
                dato = strjoin(data(o,:),',');
                fprintf(fid,'%s\r\n',dato);
            end
            fclose(fid);
            clear data dato fid
        end
%% Experiment 4 - Clustered and Overlapping Objects
    case 4
        % Create File for Results
        headers = {'distx','disty',...
            'HEDAR:n','HEDAR:Jaccard',...
            'Hough:n','Hough:Jaccard',...
            'EFT:n','EFT:Jaccard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment4.dat','w');
        fprintf(fid,'%s\r\n',headers);
        fclose(fid);
        distmax = 105;
        distmin = -105;
        major=round(30); minor=round(20);
        xpos = 128; ypos = 128;
        parfor l=distmin:distmax
            for k=distmin:distmax
                % Create Figure
                sepx = (major+minor) * l/100;
                sepy = (major+minor) * k/100;
                idx1 = ellipse2(256,[xpos,ypos],major,minor,0);
                idx2 = ellipse2(256,[xpos+sepx,ypos+sepy],major,minor,90);
                bw = max(idx1,idx2);
                % Data
                data = cell(1,8);
                data{1} = num2str(l);
                data{2} = num2str(k);
                % Run HEDAR
                results = hedar (bw,[10,40],1,90);
                data{3} = num2str(size(results,1));
                % HEDAR Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(results,1));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Hough
                edges = edge(bw,'canny');
                results = ellipticalHough(edges,[10,40]);
                data{5}= num2str(size(results,1));
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Ellipses From Triangles (EFT)
                edges = imgradient(bw);
                results = ellipsesFromTriangles(edges,[10,40]);
                data{7}= num2str(size(results,1));
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Append Data to File
                fid = fopen('experiment4.dat','a+');
                data = strjoin(data,',');
                fprintf(fid,'%s\r\n',data);
                fclose(fid);
            end
        end
%% Experiment 5 - Accuracy over lengths and rotations
    case 5
        % Data File Set-Up
        headers = {'X','Y','Major','Minor','Rotation',...
            'HEDAR:n','HEDAR:Jaccard',...
            'Hough:n','Hough:Jaccard',...
            'EFT:n','EFT:Jaccard'};
        headers = strjoin(headers,','); fid = fopen('experiment5.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        major = 30;%Only use one major axis value for accuracy maps
        parfor rot=0:179
            % for major=3:30
                for minor=1:major
                    % Data Set-Up
                    data = cell(1,11);
                    data{1} = num2str(ceil((64+1)/2));
                    data{2} = num2str(ceil((64+1)/2));
                    data{3} = num2str(major);
                    data{4} = num2str(minor);
                    data{5} = num2str(rot);
                    % Create Image
                    bw = ellipse2(64,[ceil((64+1)/2),ceil((64+1)/2)],major,minor,rot);
                    % Run HEDAR
                    results = hedar (bw);
                    l = size(results,1); data{6} = num2str(l);
                    % HEDAR Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{7} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Run Hough
                    edges = edge(bw,'canny'); results = ellipticalHough(edges);
                    l = size(results,1); data{8} = num2str(l);
                    % Hough Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,l));
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{9} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Run Ellipses From Triangles (EFT)
                    edges = imgradient(bw);
                    try
                        results = ellipsesFromTriangles(edges);
                        l = size(results,1);
                    catch
                        l=0;
                    end
                    data{10} = num2str(l);
                    % EFT Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,l));
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Append Data to File
                    fid = fopen('experiment5.dat','a+');
                    data = strjoin(data,',');
                    fprintf(fid,'%s\r\n',data);
                    fclose(fid);
                end
            % end
        end
%% Experiment 6 - Change SNR (binary)
    case 6
        % Inputs
        if nargin==1; disp('Running only once'); N = 1; end;
        % Create File for Results
        headers = {'N','SNR-Theoretical','SNR-Real',...
            'HEDAR:Jaccard','HEDAR:Time',...
            'Hough:Jaccard','Hough:Time',...
            'EFT:Jaccard','EFT:Time'};
        headers = strjoin(headers,',');
        fid = fopen('experiment6.dat','w');
        fprintf(fid,'%s\r\n',headers);
        fclose(fid);
        clear headers fid
        % Set-Up
        major=randi(20,1)+10; minor=randi(major-10,1)+10;
        rot=randi(180,1); m=128; k=35;
        for snr=34:-1:0
            % Data
            data = cell(N,9);
            data(:,2) = cellstr(num2str(repmat(snr,N,1)));
            % Create Image
            bw = ellipse2(m,[ceil((m+1)/2),ceil((m+1)/2)],major,minor,rot);
            for rn = 1:N
                % Apply Noise
                noise = (10^(-(snr-5)/20) * (randn(size(bw))));
                bw1 = bw+noise;
                bw1 = imadjust(bw1);
                [x,y] = ndgrid(1:30,1:30);
                RMS_noise = sqrt(mean(var(bw1(x,y))));
                actualSNR = 20 * log10(1/RMS_noise);
                % Blur (to remove noise)
                bw1 = imgaussfilt(bw1,3);
                % Data
                data{rn,1} = num2str(rn);
                data{rn,3} = num2str(actualSNR);
                % Run HEDAR
                tic; results = hedar (bw1,[5,k],1,1);
                data{rn,5} = num2str(toc);
                % HEDAR Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(results,1));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{rn,4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Hough
                tic; edges = edge(bw1,'canny');
                results = ellipticalHough (edges,[5,k]);
                data{rn,7} = num2str(toc);
                clear edges
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{rn,6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Ellipses From Triangles (EFT)
                tic; edges = imgradient(bw);
                results = ellipsesFromTriangles(edges);
                data{rn,9} = num2str(toc);
                clear edges
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,max(1,size(results,1)));
                for o=1:size(results,1)
                    bwo(:,:,o) = ellipse2(size(bw),[results(o,1),results(o,2)],results(o,3),results(o,4),results(o,5));
                end
                bwo = max(bwo,[],3);
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            end
            % Append Data to File
            fid = fopen('experiment6.dat','a+');
            for o=1:N
                dato = strjoin(data(o,:),',');
                fprintf(fid,'%s\r\n',dato);
            end
            fclose(fid);
        end
end
end
