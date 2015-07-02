function PODExperiments(en,N)
%% Experiments for 'A New Method for Ellipse Detection' by Nelson et al. (2015)
%
% PODExperiments(en,N)
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
% PODExperiments(4), runs the experiments used to create Figure 7
% PODExperiments(2,10), runs the experiments used to create Figure 6b ten times
%
% Copyright 2015 Carl J. Nelson, Durham University, UK
% 
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="https://github.com/ChasNelson1990/...
%              A-New-Method-for-Ellipse-Detection-2015/">The GitHub
%              Repository</a>
%
% See also POD2, PODH, ELLIPSEDETECTION, ACCURACYMAPSCRIPTS

%% Set-Up
profile -memory
switch en
%% Experiment 1 - Change of image size, set kernel size range (binary image)
    case 1
        % Inputs
        if nargin==1; N = 1; end;
        % Create File for Results
        headers = {'N','m','POD:time','POD:memory','POD:Jaqqard','PODH:time',...
            'PODH:memory','PODH:Jaqqard','Hough:time','Hough:memory','Hough:Jaqqard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment1.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        % Set-Up
        major=randi(27,1)+3; minor=randi(major-3,1)+3;
        rot=randi(180,1); ms = [64,128,256,512,1024];
        for m=1:length(ms)
            % Data
            data = cell(N,11);
            data(:,2) = cellstr(num2str(repmat(ms(m),N,1)));
            % Create Image
            bw = ellipse2(ms(m),[ceil((ms(m)+1)/2),ceil((ms(m)+1)/2)],major,minor,rot);
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run POD
                profile on, pod = pod2 (bw,(2*15)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'pod2'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,3} = num2str(time); data{rn,4} = num2str(memory);
                profile off, clear stats funct time memory
                % POD Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod;
                data{rn,5} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run PODH
                profile on, pod = podh (bw,(2*15)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'podh'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,6} = num2str(time); data{rn,7} = num2str(memory);
                profile off, clear stats funct time memory
                % PODH Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                profile on, edges = edge(bw,'canny');
                hough = ellipseDetection(edges,(2*15)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'ellipseDetection'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,9} = num2str(time); data{rn,10} = num2str(memory);
                profile off, clear stats funct time memory edges
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(hough,1));
                for o=1:size(hough,1)
                    bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
                end
                bwo = max(bwo,[],3); clear o hough
                data{rn,11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
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
        headers = {'N','k','POD:time','POD:memory','POD:Jaqqard','PODH:time',...
            'PODH:memory','PODH:Jaqqard','Hough:time','Hough:memory','Hough:Jaqqard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment2.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        % Set-Up
        ks = 15:5:65;
        for k=1:length(ks)
            % Set-Up
            major=randi(7,1)+7; minor=randi(major-1,1)+1;
            rot=randi(180,1); 
            % Create Image
            bw = ellipse2(256,[ceil((256+1)/2),ceil((256+1)/2)],major,minor,rot);
            % Data
            data = cell(N,11);
            data(:,2) = cellstr(num2str(repmat(ks(k),N,1)));
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run POD
                profile on, pod = pod2 (bw,ks(k));
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'pod2'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,3} = num2str(time); data{rn,4} = num2str(memory);
                profile off, clear stats funct time memory
                % POD Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod;
                data{rn,5} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run PODH
                profile on, pod = podh (bw,ks(k));
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'podh'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,6} = num2str(time); data{rn,7} = num2str(memory);
                profile off, clear stats funct time memory
                % PODH Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod
                data{rn,8} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                profile on, edges = edge(bw,'canny');
                hough = ellipseDetection(edges,ks(k));
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'ellipseDetection'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,9} = num2str(time); data{rn,10} = num2str(memory);
                profile off, clear stats funct time memory edges
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(hough,1));
                for o=1:size(hough,1)
                    bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
                end
                bwo = max(bwo,[],3); clear o hough
                data{rn,11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
            end
            % Append Data to File
            fid = fopen('experiment2.dat','a+');
            for o=1:N
                dato = strjoin(data(o,:),',');
                fprintf(fid,'%s\r\n',dato);
            end
            fclose(fid); clear data dato fid bw
        end
%% Experiment 3 - Different numbers of objects (non-overlapping)
    case 3
        % Inputs
        if nargin==1; N = 1; end
        % Create File for Results
        headers = {'N','number','POD:n','POD:time','POD:memory',...
            'POD:Jaqqard','PODH:n','PODH:time','PODH:memory','PODH:Jaqqard',...
            'Hough: n','Hough:time','Hough:memory','Hough:Jaqqard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment3.dat','w');
        fprintf(fid,'%s\r\n',headers);
        fclose(fid);
        for l=1:200;
            % Create Figure
            bw = zeros(256,256,l);
            xposs = cell(l,1); yposs = cell(l,1); majors = cell(l,1);
            for o=1:l
                major=randi(11,1)+3; minor=randi(major-3,1)+3; rot=randi(180,1);
                xpos = randi(226)+15; ypos = randi(226)+15;
                if o>1
                    overlapping = true;
                    while overlapping
                        overlapping = false;
                        for ps=2:o
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
                bw(:,:,o) = ellipse2(256,[xpos,ypos],major,minor,rot);
                xposs{o} = xpos; yposs{o} = ypos; majors{o} = major;
            end
            bw = max(bw,[],3);
            clear o major minor rot xpos ypos xposs yposs majors
            % Data
            data = cell(N,14);
            data(:,2) = cellstr(num2str(repmat(l,N,1)));
            for rn=1:N
                % Data
                data{rn,1} = num2str(rn);
                % Run POD
                profile on, pod = pod2 (bw,(2*7)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'pod2'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,3} = num2str(size(pod,1));
                data{rn,4} = num2str(time); data{rn,5} = num2str(memory);
                profile off, clear stats funct time memory
                % POD Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod;
                data{rn,6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run PODH
                profile on, pod = podh (bw,(2*7)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'podh'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,7} = num2str(size(pod,1));
                data{rn,8} = num2str(time); data{rn,9} = num2str(memory);
                profile off, clear stats funct time memory
                % PODH Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3); clear o pod
                data{rn,10} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                clear bwo
                % Run Hough
                profile on, edges = edge(bw,'canny');
                hough = ellipseDetection(edges,(2*7)+1);
                stats = profile('info');
                funct = find(cellfun(@(x)isequal(x,'ellipseDetection'),{stats.FunctionTable.FunctionName}));
                time = stats.FunctionTable(funct).TotalTime;
                memory = stats.FunctionTable(funct).TotalMemAllocated;
                data{rn,11} = num2str(size(hough,1));
                data{rn,12} = num2str(time); data{rn,13} = num2str(memory);
                profile off, clear stats funct time memory edges
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(hough,1));
                for o=1:size(hough,1)
                    bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
                end
                bwo = max(bwo,[],3); clear o hough
                data{rn,14} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
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
        headers = {'distx','disty','POD:n','POD:Jaqqard','PODH:n','PODH:Jaqqard','Hough: n','Hough:Jaqqard'};
        headers = strjoin(headers,',');
        fid = fopen('experiment4.dat','w');
        fprintf(fid,'%s\r\n',headers);
        fclose(fid);
%         [L,K] = readexp4();
        distmax = 150;
        distmin = -150;
        major=round(2*50/3); minor=round(50/3);
        xpos = 257; ypos = 257;
        parfor l=distmin:distmax
            for k=distmin:distmax
%                 test = sum(K(L==l)==k);
%                 if test>0; continue; end;
                % Create Figure
                sepx = (major+minor) * l/100;
                sepy = (major+minor) * k/100;
                idx1 = ellipse2(512,[xpos,ypos],major,minor,0);
                idx2 = ellipse2(512,[xpos+sepx,ypos+sepy],major,minor,90);
                bw = max(idx1,idx2);
                % Data
                data = cell(1,8);
                data{1} = num2str(l);
                data{2} = num2str(k);
                % Data
                % Run POD
                pod = pod2 (bw,[15,35],1,45);
                data{3} = num2str(size(pod,1));
                % POD Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3);
                data{4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run PODH
                pod = podh (bw,[15,40],1,1);
                data{5} = num2str(size(pod,1));
                % PODH Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3);
                data{6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Hough
                edges = edge(bw,'canny');
                hough = ellipseDetection(edges,[15,40]);
                data{7}= num2str(size(hough,1));
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(hough,1));
                for o=1:size(hough,1)
                    bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
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
        headers = {'X','Y','Major','Minor','Rotation','POD:n',...
            'POD:Jaqqard','PODH:n','PODH:Jaqqard','Hough:n','Hough:Jaqqard'};
        headers = strjoin(headers,','); fid = fopen('experiment5.dat','w');
        fprintf(fid,'%s\r\n',headers); fclose(fid);
        parfor rot=0:179
            for major=3:30
                for minor=3:major
                    % Data Set-Up
                    data = cell(1,11);
                    data{1} = num2str(ceil((64+1)/2));
                    data{2} = num2str(ceil((64+1)/2));
                    data{3} = num2str(major);
                    data{4} = num2str(minor);
                    data{5} = num2str(rot);
                    % Create Image
                    bw = ellipse2(64,[ceil((64+1)/2),ceil((64+1)/2)],major,minor,rot);
                    % Run POD
                    pod = pod2 (bw);
                    l = size(pod,1); data{6} = num2str(l);
                    % POD Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{7} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Run PODH
                    pod = podh (bw);
                    l = size(pod,1); data{8} = num2str(l);
                    % PODH Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{9} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Run Hough
                    edges = edge(bw,'canny'); hough = ellipseDetection(edges);
                    l = size(hough,1); data{10} = num2str(l);
                    % Hough Jaccard
                    bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
                    for o=1:l
                        bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
                    end
                    bwo = max(bwo,[],3);
                    data{11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                    % Append Data to File
                    fid = fopen('experiment5.dat','a+');
                    data = strjoin(data,',');
                    fprintf(fid,'%s\r\n',data);
                    fclose(fid);
                end
            end
        end
%% Experiment 6 - Change SNR (binary)
    case 6
        % Inputs
        if nargin==1; N = 1; end
        % Create File for Results
        headers = {'N','SNR-Theoretical','SNR-Real','POD:Jaqqard','POD:Time','PODH:Jaqqard','PODH:Time','Hough:Jaqqard','Hough:Time'};
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
                % Run POD
                tic; pod = pod2 (bw1,[5,k],1,1);
                data{rn,5} = num2str(toc);
                % POD Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3);
                data{rn,4} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run PODH
                tic; pod = podh (bw1,[5,k],1,1);
                data{rn,7} = num2str(toc);
                % PODH Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(pod,1));
                for o=1:size(pod,1)
                    bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                end
                bwo = max(bwo,[],3);
                data{rn,6} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
                % Run Hough
                tic; edges = edge(bw1,'canny');
                hough = ellipseDetection (edges,[5,k]);
                data{rn,9} = num2str(toc);
                % Hough Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,size(hough,1));
                for o=1:size(hough,1)
                    bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
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