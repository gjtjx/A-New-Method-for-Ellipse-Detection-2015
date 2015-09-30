% Data File Set-Up
headers = {'X','Y','Major','Minor','Rotation',...%'POD:n','POD:Jaqqard','PODH:n','PODH:Jaqqard','Hough:n','Hough:Jaqqard',...
    'EFT:n','EFT:Jaqqard'};
headers = strjoin(headers,','); fid = fopen('experiment5c.dat','w');
fprintf(fid,'%s\r\n',headers); fclose(fid);
parfor rot=0:179
    for major=30
        for minor=3:major
            % Data Set-Up
            data = cell(1,7);%11);
            data{1} = num2str(ceil((64+1)/2));
            data{2} = num2str(ceil((64+1)/2));
            data{3} = num2str(major);
            data{4} = num2str(minor);
            data{5} = num2str(rot);
            % Create Image
            bw = ellipse2(64,[ceil((64+1)/2),ceil((64+1)/2)],major,minor,rot);
%             % Run POD
%             pod = pod2 (bw);
%             l = size(pod,1); data{6} = num2str(l);
%             % POD Jaccard
%             bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
%             for o=1:l
%                 bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
%             end
%             bwo = max(bwo,[],3);
%             data{7} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
%             % Run PODH
%             pod = podh (bw);
%             l = size(pod,1); data{8} = num2str(l);
%             % PODH Jaccard
%             bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
%             for o=1:l
%                 bwo(:,:,o) = ellipse2(size(bw),[pod(o,1),pod(o,2)],pod(o,3),pod(o,4),pod(o,5));
%             end
%             bwo = max(bwo,[],3);
%             data{9} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
%             % Run Hough
%             edges = edge(bw,'canny'); hough = ellipticalHough(edges);
%             l = size(hough,1); data{10} = num2str(l);
%             % Hough Jaccard
%             bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
%             for o=1:l
%                 bwo(:,:,o) = ellipse2(size(bw),[hough(o,1),hough(o,2)],hough(o,3),hough(o,4),hough(o,5));
%             end
%             bwo = max(bwo,[],3);
%             data{11} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            % Run EFT
            edges = edge(bw,'canny'); parameters = ellipsesFromTriangles(edges);
            l = size(parameters,1); data{6} = num2str(l);
            if l~=0
                % EFT Jaccard
                bwo = zeros(size(bw)); bwo = repmat(bwo,1,1,l);
                for o=1:l
                    bwo(:,:,o) = ellipse2(size(bw),[parameters(o,1),parameters(o,2)],parameters(o,3),parameters(o,4),parameters(o,5));
                end
                bwo = max(bwo,[],3);
                data{7} = num2str(sum(sum(imabsdiff(bw,bwo)))/sum(bw(:)|bwo(:)));
            else
                data{7} = num2str(1);
            end
            % Append Data to File
            fid = fopen('experiment5c.dat','a+');
            data = strjoin(data,',');
            fprintf(fid,'%s\r\n',data);
            fclose(fid);
        end
    end
end