% Written by Seth Konig June 23, 2015
% code analyzes viewing behavior from VR rigs when monkeys viewed SCM sets.

%---[1] Find ROIs---%
% run scmgui_RM to find ROIs
% code has been slighlty modified to accommadate different name schemes

%---[2] Import Eye Data---%
all_data_dir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';
all_data_files = {'15_06_16_13_50','15_06_16_14_32','15_06_17_13_57','15_06_17_12_55',...
    '15_06_18_14_38','15_06_22_14_53','15_06_23_14_09','15_06_24_13_01',...
    '15_06_29_13_27','15_06_29_14_12'};
monk = 'JN';


% all_data_dir = 'R:\Buffalo Lab\VR Task Data UW\Gromit\panda data\';
% all_data_files = {'15_06_19_11_24','15_06_22_10_28','15_06_23_10_48',...
%     '15_06_24_09_21','15_06_29_10_19'};
% monk = 'GR';

% for file = length(all_data_files)
%     data_dir = [all_data_dir  monk '_' all_data_files{file}(1:8) '\'];
%     ImportVREyeData(all_data_files{file},data_dir)
% end
%%

imageX = 800;
imageY = 600;

img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Figures\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Eye Data\';

allnumFixationsInROI = cell(length(all_data_files),4); %number of fixations in ROI by session by type
num10FixationsInROI = cell(length(all_data_files),4); %number of fixations in ROI by session by type for 1st 10 fixations

all_sets = NaN(1,length(all_data_files));

all_time_by_time = cell(length(all_data_files),5);

allTime_ROI_timeWindow =  cell(length(all_data_files),4);%amount of time in ROI by session by type
Time3_ROI_timeWindow =  cell(length(all_data_files),4);%amount of time in ROI by session by type in 1st 3 secs
area = cell(length(all_data_files),4);

for file = 1:length(all_data_files)
    load([data_dir all_data_files{file} '-fixation.mat']);
    load([img_dir 'VRset' num2str(setnum) '_ROIs.mat']);
    figure_dir2 = [figure_dir 'VRSet' num2str(setnum) '\'];
    
    all_sets(file) = setnum;
    
    all_time = zeros(5,5000);%novel, familiar, replaced, moved,new_moved
    all_time_points = zeros(5,5000);
    
    for img = 1:size(pairings,2)
        if ~isnan(pairings(1,img))
            nov_ind = pairings(1,img);%1st presentation (nov)
            rep_ind = pairings(2,img);%2nd presentation  (rep)
            
            rep_name = image_name{rep_ind};
            if strcmpi(rep_name(8),'p')
                trialtype = 2;
            elseif strcmpi(rep_name(8),'r')
                trialtype = 3;
            elseif strcmpi(rep_name(8),'m')
                trialtype = 4;
            else
                error('image type unkown')
            end
            
            nov_time_out = zeros(1,5000); %time outside image or crosshair fixation
            novx = fixationstats{nov_ind}.XY(1,1:5000);
            novy = fixationstats{nov_ind}.XY(2,1:5000);
            novfixations = fixationstats{nov_ind}.fixations;
            novtimes =  fixationstats{nov_ind}.fixationtimes;
            novtimes= ceil(novtimes/5);
            if (novfixations(1,1) > 300 && novfixations(1,1) < 500) && ...
                    (novfixations(2,1) > 200 && novfixations(2,1) < 400)
                nov_time_out(1:novtimes(2,1)) = 1;
                novx(1:novtimes(2,1)) = [];
                novy(1:novtimes(2,1)) = [];
                novfixations(:,1) = [];
                novtimes(:,1) = [];
            end
            nov_outside = find(novx > imageX | novx < 1 | novy > imageX | novy < 1);
            nov_time_out(nov_outside) = 1;
            novtimes = novtimes*5;
            
            rep_time_out = zeros(1,5000);
            if length(fixationstats{rep_ind}.XY) > 5000
                repx = fixationstats{rep_ind}.XY(1,1:5000);
                repy = fixationstats{rep_ind}.XY(2,1:5000);
            end
            repfixations = fixationstats{rep_ind}.fixations;
            repfixations(1,:) =  repfixations(1,:);
            repfixations(2,:) = repfixations(2,:);
            reptimes =  fixationstats{rep_ind}.fixationtimes;
            reptimes= ceil(reptimes/5);
            if (repfixations(1,1) > 300 && repfixations(1,1) < 500) && ...
                    (repfixations(2,1) > 200 && repfixations(2,1) < 400)
                rep_time_out(1:reptimes(2,1)) = 1;
                repfixations(:,1) = [];
                reptimes(:,1) = [];
                repx(1:reptimes(2,1)) = [];
                repy(1:reptimes(2,1)) = [];
            end
            rep_outside = find(repx > imageX | repx < 1 | repy > imageX | repy < 1);
            rep_time_out(rep_outside) = 1;
            reptimes = reptimes*5;
            
            novy(novx > imageX) = [];
            novx(novx > imageX) = [];
            novx(novy > imageY) = [];
            novy(novy > imageY) = [];
            novy(novx < 1) = [];
            novx(novx < 1) = [];
            novx(novy < 1) = [];
            novy(novy < 1) = [];
            
            repy(repx > imageX) = [];
            repx(repx > imageX) = [];
            repx(repy > imageY) = [];
            repy(repy > imageY) = [];
            repy(repx < 1) = [];
            repx(repx < 1) = [];
            repx(repy < 1) = [];
            repy(repy < 1) = [];
            
            %ignore anything after 5000 ms
            late = find(novtimes(1,:) > 5000);
            novtimes(:,late) = [];
            novfixations(:,late) = [];
            novtimes(novtimes > 5000) = 5000;
            
            late = find(reptimes(1,:) > 5000);
            reptimes(:,late) = [];
            repfixations(:,late) = [];
            reptimes(reptimes > 5000) = 5000;
            
            if trialtype == 2 || trialtype == 3 %familiar and replaced, respectively
                ROI1 = ROIs{img}{1};
                ROI1(1) = ROI1(1)-36;
                ROI1(2) = ROI1(2)+36;
                ROI1(3) = ROI1(3)-36;
                ROI1(4) = ROI1(4)+36;
                ROI1(ROI1 < 1) = 1;
                ROI1(ROI1 > imageX) = imageX;
                if ROI1(3) > imageY; ROI1(3) = imageY;end
                if ROI1(4) > imageY; ROI1(4) = imageY;end
                ROI2 =ROI1;
            elseif trialtype == 4 %moved
                ROI1 = ROIs{img}{1};%original location
                ROI1(1) = ROI1(1)-36;
                ROI1(2) = ROI1(2)+36;
                ROI1(3) = ROI1(3)-36;
                ROI1(4) = ROI1(4)+36;
                ROI1(ROI1 < 1) = 1;
                ROI1(ROI1 > imageX) = imageX;
                if ROI1(3) > imageY; ROI1(3) = imageY;end
                if ROI1(4) > imageY; ROI1(4) = imageY;end
                ROI2 = ROIs{img}{2};%new location after move
                ROI2(1) = ROI2(1)-36;
                ROI2(2) = ROI2(2)+36;
                ROI2(3) = ROI2(3)-36;
                ROI2(4) = ROI2(4)+36;
                ROI2(ROI2 < 1) = 1;
                ROI2(ROI2 > imageX) = imageX;
                if ROI2(3) > imageY; ROI2(3) = imageY;end
                if ROI2(4) > imageY; ROI2(4) = imageY;end
            end
            
            
            rep_trialType = trialtype;
            
            ROIfix1 = find(...
                (novfixations(1,:) > ROI1(1) & novfixations(1,:) < ROI1(2)) & ...
                (novfixations(2,:) > ROI1(3) & novfixations(2,:) < ROI1(4))); %ROIfix1 contains indices of novfixations corresponding to fixations that occur in ROI of novel image.
            ROIfix2 = find(...
                repfixations(1,:) > ROI2(1) & repfixations(1,:) < ROI2(2) & ...
                repfixations(2,:) > ROI2(3) & repfixations(2,:) < ROI2(4));%ROIfix2 contains indices of repfixations corresponding to fixations that occur in ROI of the manipulated image
            
            if trialtype ==4 %moved image. Data from the old ROI (where the object was before being moved)
                ROIfix2_oldROI = find(...
                    (repfixations(1,:) > ROI1(1) & repfixations(1,:) < ROI1(2) & ...
                    repfixations(2,:) > ROI1(3) & repfixations(2,:) < ROI1(4)));%ROIfix2_oldROI contains indices of repfixations corresponding to fixations that occur in the 'old' ROI of the 'moved' image
            else
                ROIfix2_oldROI = [];
            end
            
            %added commented out code below to plot eye traces on images
            
            %             if ~isempty(ROIfix2_oldROI)
            %                 disp('found fixations')
            %             end
            figure
            subplot(2,2,1)
            hold on
            image(flipud(imread([img_dir 'VRset' num2str(setnum) '\' image_name{nov_ind}])));
            plot(novx,novy);
            plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
                [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
            plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
                [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
            hold off
            box off
            xlim([0 imageX])
            ylim([0 imageY])
            axis off
            axis equal
            
            subplot(2,2,3)
            hold on
            image(flipud(imread([img_dir 'VRset' num2str(setnum) '\' image_name{rep_ind}])));
            plot(repx,repy);
            plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
                [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
            plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
                [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
            hold off
            box off
            xlim([0 imageX])
            ylim([0 imageY])
            axis off
            axis equal
            
            if isempty(ROIfix1)%didnt look at original region, moving on to next image
                subplot(2,2,2)
                title('No fixations in Novel ROI')
                save_and_close_fig(figure_dir2,[monk num2str(img)]);
                continue
            end
            
            
            
            %%%%Code for calculating avg time in ROI for each image type
            tempvec1=zeros(1,5000);%novel image
            for j=1:length(ROIfix1)
                tempvec1(novtimes(1,ROIfix1(j)):novtimes(2,ROIfix1(j)))=1;
            end
            tempvec1(find(nov_time_out)) = NaN;
            
            tempvec2=zeros(1,5000);%manipulated image new ROI
            for j=1:length(ROIfix2)
                tempvec2(reptimes(1,ROIfix2(j)):reptimes(2,ROIfix2(j)))=1;
            end
            tempvec2(find(rep_time_out)) = NaN;
            
            tempvec3=zeros(1,5000);%manipulated imaged old moved ROI
            for j=1:length(ROIfix2_oldROI)
                tempvec3(reptimes(1,ROIfix2_oldROI(j)):reptimes(2,ROIfix2_oldROI(j)))=1;
            end
            tempvec3(find(rep_time_out)) = NaN;
            
            area_ROI1 = abs(ROI1(1)-ROI1(2))*abs(ROI1(3)-ROI1(4)); %area of ROI for novel image
            area_ROI2 = abs(ROI2(1)-ROI2(2))*abs(ROI2(3)-ROI2(4)); %for replaced, moved or familiar
            
            %for novel presentation
            area{file,1} = [area{file,1} area_ROI1];
            
            propall_fixations = length(ROIfix1)/length(novfixations);
            prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
            propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
            prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
            
            allnumFixationsInROI{file,1} =[allnumFixationsInROI{file,1} propall_fixations];
            num10FixationsInROI{file,1} = [num10FixationsInROI{file,1} prop10_fixaitons];
            allTime_ROI_timeWindow{file,1} = [allTime_ROI_timeWindow{file,1} propall_time];
            Time3_ROI_timeWindow{file,1} = [Time3_ROI_timeWindow{file,1} prop3_time];
            
            s(1) = subplot(2,2,2);
            bar([propall_fixations propall_time prop10_fixaitons prop3_time])
            set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
            ylabel('Proportion in ROI')
            title(['First Fixation in ROI is # ' num2str(ROIfix1(1))])
            
            
            %for second presentation
            if trialtype == 4
                area{file,trialtype} = [area{file,trialtype} area_ROI2+area_ROI1];
                
                
                propall_fixations = (length(ROIfix2)+length(ROIfix2_oldROI))/length(repfixations);
                prop10_fixaitons = (length(ROIfix2(ROIfix2 <= 10))+...
                    length(ROIfix2_oldROI(ROIfix2_oldROI <= 10)))./10;
                propall_time = (nansum(tempvec2)+nansum(tempvec3))/sum(~isnan(tempvec2));
                prop3_time = (nansum(tempvec2(1:3000))+nansum(tempvec3(1:3000)))...
                    /sum(~isnan(tempvec2(1:3000)));
                
                allnumFixationsInROI{file,trialtype} =[allnumFixationsInROI{file,trialtype} propall_fixations];
                num10FixationsInROI{file,trialtype} = [num10FixationsInROI{file,trialtype} prop10_fixaitons];
                allTime_ROI_timeWindow{file,trialtype} = [allTime_ROI_timeWindow{file,trialtype} propall_time];
                Time3_ROI_timeWindow{file,trialtype} = [Time3_ROI_timeWindow{file,trialtype} prop3_time];
                
                s(2) = subplot(2,2,4);
                bar([propall_fixations propall_time prop10_fixaitons prop3_time])
                set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
                ylabel('Proportion in ROI')
                if ~isempty(ROIfix2)
                    title(['First Fixation in ROI is # ' num2str(ROIfix2(1))])
                end
                
            else %for replaced or repeat
                area{file,trialtype} = [area{file,trialtype} area_ROI2];
                
                propall_fixations = length(ROIfix2)/length(repfixations);
                prop10_fixaitons = length(ROIfix2(ROIfix2 <= 10))/10;
                propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
                prop3_time = nansum(tempvec2(1:3000))/sum(~isnan(tempvec2(1:3000)));
                
                allnumFixationsInROI{file,trialtype} =[allnumFixationsInROI{file,trialtype} propall_fixations];
                num10FixationsInROI{file,trialtype} = [num10FixationsInROI{file,trialtype} prop10_fixaitons];
                allTime_ROI_timeWindow{file,trialtype} = [allTime_ROI_timeWindow{file,trialtype} propall_time];
                Time3_ROI_timeWindow{file,trialtype} = [Time3_ROI_timeWindow{file,trialtype} prop3_time];
                
                s(2) = subplot(2,2,4);
                bar([propall_fixations propall_time prop10_fixaitons prop3_time])
                set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
                ylabel('Proportion in ROI')
                if ~isempty(ROIfix2)
                    title(['First Fixation in ROI is # ' num2str(ROIfix2(1))])
                end
            end
            yl1 = get(s(1),'ylim');
            yl2 = get(s(2),'ylim');
            ymax = max([yl1 yl2]);
            set(s(1),'ylim',[0 ymax]);
            set(s(2),'ylim',[0 ymax]);
            
            save_and_close_fig(figure_dir2,[monk num2str(img)]);
            
            
            all_time_points(1,~isnan(tempvec1)) = all_time_points(1,~isnan(tempvec1))+1;
            tempvec1(isnan(tempvec1)) = 0;
            all_time(1,:) = all_time(1,:) + tempvec1;
            
            all_time_points(trialtype,~isnan(tempvec2)) = ...
                all_time_points(trialtype,~isnan(tempvec2))+1;
            tempvec2(isnan(tempvec2)) = 0;
            all_time(trialtype,:) =  all_time(trialtype,:)+tempvec2;
            
            if trialtype == 4 %for moved only
                all_time_points(5,~isnan(tempvec3)) = ...
                    all_time_points(5,~isnan(tempvec3))+1;
                tempvec3(isnan(tempvec3)) = 0;
                all_time(5,:) =  all_time(5,:)+tempvec3+tempvec2;
            end
        end
    end
    
    all_time_points(all_time_points == 0) = 1; %for division purposes
    all_time = all_time./all_time_points; %normalize
    
    for t=1:5
        all_time_by_time{file,t} = all_time(t,:);
    end
    %
    %     clrs = 'brmgg';
    %     line = {'-','-','-','-','--'};
    %     image_count = cellfun(@length,allnumFixationsInROI(file,:));%num images in each set
    %     image_count = [image_count image_count(end)];
    %     figure
    %     hold on
    %     for r = 1:size(all_time,1)
    %         plot(100*filtfilt(1/200*ones(1,200),1,all_time(r,:)),...
    %             [clrs(r) line{r}]);
    %     end
    %     hold off
    %     xlabel('Time (ms)')
    %     ylabel('% of time in ROI')
    %     legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
    %     title(['Set ' num2str(setnum)])
end
%%
clrs = 'brmgg';
line = {'-','-','-','-','--'};
all_sets_avg = zeros(5,5000);
for a = min(all_sets):max(all_sets);
    this_set = find(all_sets == a);
    set_avg = zeros(5,5000);
    for t=1:length(this_set)
        for tp = 1:5
            all_sets_avg(tp,:) =  all_sets_avg(tp,:)+ all_time_by_time{this_set(t),tp};
            set_avg(tp,:) = set_avg(tp,:)+all_time_by_time{this_set(t),tp};
        end
    end
    set_avg = set_avg/length(this_set);
    
    
    image_count = cellfun(@length,allnumFixationsInROI(file,:));%num images in each set
    image_count = [image_count image_count(end)];
    figure
    hold on
    for r = 1:size(set_avg,1)
        plot(100*filtfilt(1/200*ones(1,200),1,set_avg (r,:)),...
            [clrs(r) line{r}]);
    end
    hold off
    xlabel('Time (ms)')
    ylabel('% of time in ROI')
    legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
    title(['Set ' num2str(a)])
end
all_sets_avg = all_sets_avg/length(all_sets);

figure
hold on
for r = 1:size(set_avg,1)
    plot(100*filtfilt(1/200*ones(1,200),1,all_sets_avg (r,:)),...
        [clrs(r) line{r}]);
end
hold off
xlabel('Time (ms)')
ylabel('% of time in ROI')
legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
title(['Average across all sets, n = ' num2str(max(all_sets))])
%%
sess_num_fix =NaN(max(all_sets),4);
sess_num10fix = NaN(max(all_sets),4);
sess_allTime = NaN(max(all_sets),4);
sess_Time3 = NaN(max(all_sets),4);
for a = min(all_sets):max(all_sets);
    this_set = find(all_sets == a);
    set_avg_nf = [];
    set_avg_nf10 = [];
    set_avg_at = [];
    set_avg_at3 = [];
    for t=1:length(this_set)
        for tp = 1:4
            set_avg_nf(t,tp) = mean(allnumFixationsInROI{this_set(t),tp});
            set_avg_nf10(t,tp) = mean(num10FixationsInROI{this_set(t),tp});
            set_avg_at(t,tp) = mean(num10FixationsInROI{this_set(t),tp});
            set_avg_at3(t,tp) = mean(num10FixationsInROI{this_set(t),tp});
        end
    end
    
    for tp = 1:4
        sess_num_fix(a,tp) = 100*mean(set_avg_nf(:,tp));
        sess_num10fix(a,tp) = 100*mean(set_avg_nf10(:,tp));
        sess_allTime(a,tp)= 100*mean(set_avg_at(:,tp));
        sess_Time3(a,tp)= 100*mean(set_avg_at3(:,tp));
    end
end

figure
subplot(2,2,1)
hold on
bar(nanmean(sess_num10fix));
errorb(nanmean(sess_num10fix),nanstd(sess_num10fix)./sqrt(sum(~isnan(sess_num10fix(:,1)))))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of 1st 10 Fixations in ROI')

subplot(2,2,2)
hold on
bar(nanmean(sess_Time3));
errorb(nanmean(sess_Time3),nanstd(sess_Time3)./sqrt(sum(~isnan(sess_Time3(:,1)))))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of Time in 1st 3 seconds in ROI')

subplot(2,2,3)
hold on
bar(nanmean(sess_num_fix));
errorb(nanmean(sess_num_fix),nanstd(sess_num_fix)./sqrt(sum(~isnan(sess_num_fix(:,1)))))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Fixations in ROI')

subplot(2,2,4)
hold on
bar(nanmean(sess_allTime));
errorb(nanmean(sess_allTime),nanstd(sess_allTime)./sqrt(sum(~isnan(sess_allTime(:,1)))))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Time in ROI')

subtitle(monk)