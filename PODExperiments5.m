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
            disp(exist('bw','var'))
            % Run POD
            pod = pod2 (bw);
            disp(exist('bw','var'))
            l = size(pod,1); data{6} = num2str(l);
            % POD Jaccard
            bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
            disp(exist('bw','var'))
            for o=1:l
                bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                disp(exist('bw','var'))
            end
            bwo = max(bwo,[],3);
            data{7} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            disp(exist('bw','var'))
            % Run PODH
            pod = podh (bw);
            disp(exist('bw','var'))
            l = size(pod,1); data{8} = num2str(l);
            % PODH Jaccard
            bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
            disp(exist('bw','var'))
            for o=1:l
                bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
                disp(exist('bw','var'))
            end
            bwo = max(bwo,[],3);
            data{9} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            disp(exist('bw','var'))
            % Run Hough
            edges = edge(bw,'canny'); hough = ellipseDetection(edges);
            disp(exist('bw','var'))
            l = size(hough,1); data{10} = num2str(l);
            % Hough Jaccard
            bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
            disp(exist('bw','var'))
            for o=1:l
                bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
                disp(exist('bw','var'))
            end
            bwo = max(bwo,[],3);
            data{11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            disp(exist('bw','var'))
            % Append Data to File
            fid = fopen('experiment5.dat','a+');
            data = strjoin(data,',');
            fprintf(fid,'%s\r\n',data);
            fclose(fid);
        end
    end
end