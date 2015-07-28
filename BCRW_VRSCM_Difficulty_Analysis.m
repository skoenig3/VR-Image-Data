% Code written By Seth Konig to determine if certain maniuplations are
% harder than others to spot. Written July 2015

%%
%---[1] Make Salience maps for manipulated imagse---%
% image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';
%
% for SET =1:8
%     cd([image_dir 'VRset' num2str(SET) '\']);
%     list = ls;
%     for img = 1:size(list,1)
%         if ~isempty(strfind(list(img,:),'.bmp'))
%             getSalienceMap(list(img,:))
%         end
%     end
% end
% emailme('Finished creating saliencemaps')
%%
% %---[2] Run simulations for manipulated images---%
% image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';
% tags = {'MP','TT','JN','IW'};
% Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets'...
%     '\CombinedViewingBehavior.mat'];
% use behavioral files from the 4 monkeys we have lots of behavior from not
% the 1 or 2 we have little from
%
% load(Combinedbehaviorfile,'allview')
% imageX = 800; imageY = 600;
% plotoptions.runs = 'none'; %all/none
% plotoptions.probdens = 'none';
% plotoptions.type = 'sal'; %sal/image
% IOR_tau = [1/17];
% % %same name for everyset
% %
% for SET = 1:8
%    emailme(['Running Simulations for Set ' num2str(SET)])
%     cd([image_dir 'VRset' num2str(SET) '\']);
%     a = what;
%     saliencemapfiles = a.mat;
%
%     for i = 1:size(saliencemapfiles,1)
%         for t = 1:length(tags)
%             disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' num2str(SET)])
%             run_BCRWCF(allview{t},saliencemapfiles{i},tags{t},imageX,imageY,plotoptions,IOR_tau)
%         end
%     end
% end
%%
%---[3] Determine the number of fixations to Manipulated ROI---%

image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';

tags = {'MP','TT','JN','IW'};

imageY = 600;
imageX = 800;

replaced_difficulty = NaN(8,36);
moved_difficulty = NaN(8,36);

replaced_total = NaN(8,36);
moved_total = NaN(8,36);
replaced_area = NaN(8,36);
moved_area = NaN(8,36);
%area should be the same since the same object is being manipulated
for SET =1:8;
    cd([image_dir 'VRset' num2str(SET) '\']);
    load([image_dir 'VRset' num2str(SET) '_ROIs.mat'])
    
    for index = 1:36
        if exist(['VRS' num2str(SET) 'I' num2str(index) 'r-saliencemap.mat'],'file')
            s =1; %for replaced
            zero = [];
            char = 'r';
        elseif exist(['VRS' num2str(SET) 'I' num2str(index) 'm-saliencemap.mat'],'file')
            s = 2;%for moved
            zero = [];
            char = 'm';
        elseif exist(['VRS' num2str(SET) 'I0' num2str(index) 'r-saliencemap.mat'],'file')  %some sest have 0's for images < 10
            s = 1; %for replaced
            zero = '0';
            char = 'r';
        elseif exist(['VRS' num2str(SET) 'I0' num2str(index) 'm-saliencemap.mat'],'file') %some sest have 0's for images < 10
            s = 2; %for moved
            zero = '0';
            char = 'm';
        else %not a manipulated image, either repeat or novel
            continue
        end
        
        if  s == 1  %eplaced
            ROI = ROIs{index}{1};
            ROI(1) = ROI(1)-36;
            ROI(2) = ROI(2)+36;
            ROI(3) = ROI(3)-36;
            ROI(4) = ROI(4)+36;
            ROI(ROI < 1) = 1;
            ROI(ROI > imageX) = imageX;
            if ROI(3) > imageY; ROI(3) = imageY;end
            if ROI(4) > imageY; ROI(4) = imageY;end
            ROI2 =ROI;
        elseif s == 2 %moved
            ROI = ROIs{index}{1};%original location
            ROI(1) = ROI(1)-36;
            ROI(2) = ROI(2)+36;
            ROI(3) = ROI(3)-36;
            ROI(4) = ROI(4)+36;
            ROI(ROI < 1) = 1;
            ROI(ROI > imageX) = imageX;
            if ROI(3) > imageY; ROI(3) = imageY;end
            if ROI(4) > imageY; ROI(4) = imageY;end
        end
        
        for t = 1:length(tags)
            load(['BCRW IOR TAU 17\' tags{t} '-' 'VRS' num2str(SET) 'I' zero num2str(index) char '-BCRW.mat'],'fixationtimes');
            for i = 1:size(fixationtimes,1);
                fixations{t}{i} = NaN(2,60);
            end
            for i = 1:size(fixationtimes,1);
                tind = find(fixationtimes(i,:,1) > 0);
                if length(tind) > 60
                    tind = tind(1:60);
                end
                for ii = 1:length(tind)
                    x = fixationtimes(i,tind(ii),1);
                    y = imageY-fixationtimes(i,tind(ii),2);%ROI is inverted from BCRW coordinates
                    fixations{t}{i}(:,ii) = [x;y];
                end
            end
        end
        
        difficulty = NaN(4,length(fixations{1}));
        total = NaN(4,length(fixations{1}));
        for t = 1:length(fixations)
            for i = 1:length(fixations{t})
                %add 1 dva or 24 pixel buffer around ROI
                ind = find(fixations{t}{i}(1,:) > ROI(1)...
                    & fixations{t}{i}(1,:) < ROI(2) ....
                    & fixations{t}{i}(2,:) > ROI(3) ...
                    & fixations{t}{i}(2,:) < ROI(4));
                if ~isempty(ind)
                    difficulty(t,i) = min(ind);
                    total(t,i) = length(ind);
                else
                    difficulty(t,i) = NaN; %if not found cap
                    total(t,i) = 0;
                end
            end
        end
        
        %         img = zeros(600,800);
        %         for t = 1:length(fixations)
        %             for i = 1:length(fixations{t})
        %                 for ii = 1:size(fixations{t}{i},2)
        %                     if ~isnan(fixations{t}{i}(1,ii))
        %                         img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii)) =  img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii))+1;
        %                     end
        %                 end
        %             end
        %         end
        
        if s == 1
            replaced_difficulty(SET,index) = nanmean(difficulty(1:end));
            replaced_total(SET,index) = nanmean(total(1:end));
            replaced_area(SET,index) = (ROI(2)-ROI(1))*(ROI(4)-ROI(3));
        else
            moved_difficulty(SET,index) = nanmean(difficulty(1:end));
            moved_total(SET,index) = nanmean(total(1:end));
            moved_area(SET,index) =  (ROI(2)-ROI(1))*(ROI(4)-ROI(3));
        end
    end
end

nanmean(nanmean(replaced_difficulty))
nanmean(nanmean(moved_difficulty))
nanmean(nanmean(replaced_total))
nanmean(nanmean(moved_total))
nanmean(nanmean(replaced_area))
nanmean(nanmean(moved_area))

[~,tpd] = ttest2(replaced_difficulty(1:end),moved_difficulty(1:end))
[~,kpd] = kstest2(replaced_difficulty(1:end),moved_difficulty(1:end))
[~,tpt] = ttest2(replaced_total(1:end),moved_total(1:end))
[~,kpt] = kstest2(replaced_total(1:end),moved_total(1:end))
[~,tpa] = ttest2(replaced_area(1:end),moved_area(1:end))
[~,kpa] = kstest2(replaced_area(1:end),moved_area(1:end))
%%
all_data_dir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';
all_data_files = {'15_06_16_13_50','15_06_16_14_32','15_06_17_13_57','15_06_17_12_55',...
    '15_06_18_14_38','15_06_22_14_53','15_06_23_14_09','15_06_24_13_01',...
    '15_06_29_13_27','15_06_29_14_12','15_07_16_14_10',...
    '15_07_17_12_31','15_07_23_12_28','15_07_23_12_58'};
monk = 'JN';

imageX = 800;
imageY = 600;

img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Figures\';
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Eye Data\';

allnumFixationsInROI = cell(length(all_data_files),4); %number of fixations in ROI by session by type
num10FixationsInROI = cell(length(all_data_files),4); %number of fixations in ROI by session by type for 1st 10 fixations

all_sets = NaN(1,length(all_data_files));
all_img_nums =  cell(length(all_data_files),4);

all_time_by_time = cell(length(all_data_files),5);

allTime_ROI_timeWindow =  cell(length(all_data_files),4);%amount of time in ROI by session by type
Time3_ROI_timeWindow =  cell(length(all_data_files),4);%amount of time in ROI by session by type in 1st 3 secs
Time15_ROI_timeWindow =  cell(length(all_data_files),4);%amount of time in ROI by session by type in 0.5-1.5 secs window

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
            imgnum = str2double(rep_name(6:7));
            
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
            
            
            %for novel presentation
            
            propall_fixations = length(ROIfix1)/length(novfixations);
            prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
            propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
            prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
            prop15_time = nansum(tempvec1(500:1500))/sum(~isnan(tempvec1(500:1500)));
            
            allnumFixationsInROI{file,1} =[allnumFixationsInROI{file,1} propall_fixations];
            num10FixationsInROI{file,1} = [num10FixationsInROI{file,1} prop10_fixaitons];
            allTime_ROI_timeWindow{file,1} = [allTime_ROI_timeWindow{file,1} propall_time];
            Time3_ROI_timeWindow{file,1} = [Time3_ROI_timeWindow{file,1} prop3_time];
            Time15_ROI_timeWindow{file,1}= [ Time15_ROI_timeWindow{file,1} prop15_time];
            
            %for second presentation
            if trialtype == 4
                
                propall_fixations = (length(ROIfix2)+length(ROIfix2_oldROI))/length(repfixations);
                prop10_fixaitons = (length(ROIfix2(ROIfix2 <= 10))+...
                    length(ROIfix2_oldROI(ROIfix2_oldROI <= 10)))./10;
                propall_time = (nansum(tempvec2)+nansum(tempvec3))/sum(~isnan(tempvec2));
                prop3_time = (nansum(tempvec2(1:3000))+nansum(tempvec3(1:3000)))...
                    /sum(~isnan(tempvec2(1:3000)));
                prop15_time = (nansum(tempvec2(500:1500))+nansum(tempvec3(500:1500)))...
                    /sum(~isnan(tempvec2(500:1500)));
                
                allnumFixationsInROI{file,trialtype} =[allnumFixationsInROI{file,trialtype} propall_fixations];
                num10FixationsInROI{file,trialtype} = [num10FixationsInROI{file,trialtype} prop10_fixaitons];
                allTime_ROI_timeWindow{file,trialtype} = [allTime_ROI_timeWindow{file,trialtype} propall_time];
                Time3_ROI_timeWindow{file,trialtype} = [Time3_ROI_timeWindow{file,trialtype} prop3_time];
                Time15_ROI_timeWindow{file,trialtype} = [Time15_ROI_timeWindow{file,trialtype} prop15_time];
                all_img_nums{file,trialtype} = [ all_img_nums{file,trialtype} imgnum];
                
            else %for replaced or repeat
                
                propall_fixations = length(ROIfix2)/length(repfixations);
                prop10_fixaitons = length(ROIfix2(ROIfix2 <= 10))/10;
                propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
                prop3_time = nansum(tempvec2(1:3000))/sum(~isnan(tempvec2(1:3000)));
                prop15_time = nansum(tempvec2(500:1500))/sum(~isnan(tempvec2(500:1500)));
                
                allnumFixationsInROI{file,trialtype} =[allnumFixationsInROI{file,trialtype} propall_fixations];
                num10FixationsInROI{file,trialtype} = [num10FixationsInROI{file,trialtype} prop10_fixaitons];
                allTime_ROI_timeWindow{file,trialtype} = [allTime_ROI_timeWindow{file,trialtype} propall_time];
                Time3_ROI_timeWindow{file,trialtype} = [Time3_ROI_timeWindow{file,trialtype} prop3_time];
                Time15_ROI_timeWindow{file,trialtype} = [Time15_ROI_timeWindow{file,trialtype} prop15_time];
                all_img_nums{file,trialtype} = [ all_img_nums{file,trialtype} imgnum];
                
                
            end
            
            
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
end
%%
replaced_allfix = [];
replaced_10fix = [];
replaced_alltim = [];
replaced_t3 = [];
replaced_t15 = [];

moved_ar = [];
moved_allfix = [];
moved_10fix = [];
moved_alltim = [];
moved_t3 = [];
moved_t15 = [];

replaced_ar = [];
replaced_dif = [];
replaced_tot = [];
moved_dif = [];
moved_tot = [];

for file = 1:size(allnumFixationsInROI,1)
    for type = 3:4 %replaced, moved
        set = all_sets(file); %should be all same set
        imgnums = all_img_nums{file,type};
        if type == 4
            moved_ar = [moved_ar moved_area(set,imgnums)];
            moved_dif = [moved_dif moved_difficulty(set,imgnums)];
            moved_tot = [moved_tot  moved_total(set,imgnums)];
            moved_allfix = [moved_allfix allnumFixationsInROI{file,type}];
            moved_10fix = [moved_10fix num10FixationsInROI{file,type}];
            moved_alltim = [moved_alltim allTime_ROI_timeWindow{file,type}];
            moved_t3 = [moved_t3 Time15_ROI_timeWindow{file,type}];
            moved_t15 = [moved_t15  Time3_ROI_timeWindow{file,type}];
        else
            replaced_ar = [replaced_ar replaced_area(set,imgnums)];
            replaced_dif = [replaced_dif replaced_difficulty(set,imgnums)];
            replaced_tot = [replaced_tot replaced_total(set,imgnums)];
            replaced_allfix = [replaced_allfix allnumFixationsInROI{file,type}];
            replaced_10fix = [replaced_10fix num10FixationsInROI{file,type}];
            replaced_alltim = [replaced_alltim allTime_ROI_timeWindow{file,type}];
            replaced_t3 = [replaced_t3  Time3_ROI_timeWindow{file,type}];
            replaced_t15 = [replaced_t15 Time15_ROI_timeWindow{file,type}];
        end
    end
end
%%
figure

subplot(3,5,1)
plot(replaced_dif,replaced_allfix,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(replaced_dif,replaced_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,2)
plot(replaced_dif,replaced_10fix,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(replaced_dif,replaced_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,3)
plot(replaced_dif,replaced_alltim,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(replaced_dif,replaced_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,4)
plot(replaced_dif,replaced_t3,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(replaced_dif,replaced_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,5)
plot(replaced_dif,replaced_t15,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(replaced_dif,replaced_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,6)
plot(replaced_tot,replaced_allfix,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(replaced_tot,replaced_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,7)
plot(replaced_tot,replaced_10fix,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(replaced_tot,replaced_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,8)
plot(replaced_tot,replaced_alltim,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(replaced_tot,replaced_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,9)
plot(replaced_tot,replaced_t3,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(replaced_tot,replaced_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,10)
plot(replaced_tot,replaced_t15,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(replaced_tot,replaced_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,11)
plot(replaced_ar,replaced_allfix,'k.')
xlabel('ROI Area')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(replaced_ar,replaced_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,12)
plot(replaced_ar,replaced_10fix,'k.')
xlabel('ROI Area')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(replaced_ar,replaced_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,13)
plot(replaced_ar,replaced_alltim,'k.')
xlabel('ROI Area')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(replaced_ar,replaced_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,14)
plot(replaced_ar,replaced_t3,'k.')
xlabel('ROI Area')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(replaced_ar,replaced_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,15)
plot(replaced_ar,replaced_t15,'k.')
xlabel('ROI Area')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(replaced_ar,replaced_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

%subtitle('Replaced')
%%
figure

subplot(3,5,1)
plot(moved_dif,moved_allfix,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(moved_dif,moved_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,2)
plot(moved_dif,moved_10fix,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(moved_dif,moved_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,3)
plot(moved_dif,moved_alltim,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(moved_dif,moved_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,4)
plot(moved_dif,moved_t3,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(moved_dif,moved_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,5)
plot(moved_dif,moved_t15,'k.')
xlabel('BCRW first fixation')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(moved_dif,moved_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,6)
plot(moved_tot,moved_allfix,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(moved_tot,moved_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end


subplot(3,5,7)
plot(moved_tot,moved_10fix,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(moved_tot,moved_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,8)
plot(moved_tot,moved_alltim,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(moved_tot,moved_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,9)
plot(moved_tot,moved_t3,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(moved_tot,moved_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,10)
plot(moved_tot,moved_t15,'k.')
xlabel('BCRW total fixations')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(moved_tot,moved_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,11)
plot(moved_ar,moved_allfix,'k.')
xlabel('ROI Area')
ylabel('Proportion of all fix in ROI')
[r,p] = corrcoef(moved_ar,moved_allfix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,12)
plot(moved_ar,moved_10fix,'k.')
xlabel('ROI Area')
ylabel('Proportion of 1st 10 fix in ROI')
[r,p] = corrcoef(moved_ar,moved_10fix);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,13)
plot(moved_ar,moved_alltim,'k.')
xlabel('ROI Area')
ylabel('Proportion of All time in ROI')
[r,p] = corrcoef(moved_ar,moved_alltim);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,14)
plot(moved_ar,moved_t3,'k.')
xlabel('ROI Area')
ylabel('Proportion of 1st 3 seconds in ROI')
[r,p] = corrcoef(moved_ar,moved_t3);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

subplot(3,5,15)
plot(moved_ar,moved_t15,'k.')
xlabel('ROI Area')
ylabel('Proportion of 0.5-1.5 secs in ROI')
[r,p] = corrcoef(moved_ar,moved_t15);
if p(2) < 0.05/15
    title(num2str(r(2)))
end

% subtitle('moved')