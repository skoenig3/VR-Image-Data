% written by Seth Konig 5/19/2015
% Code calculates fixation durations from images viewed on the VR rig

eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Eye Data\';

data_files = {'15_01_06_14_28','15_01_06_15_28','15_01_06_16_34',...
    '15_01_07_13_15','15_01_07_14_07','15_01_08_14_40','15_01_09_12_23',...
    '15_01_09_13_36','15_01_12_12_26','15_01_12_13_22','15_01_13_15_25',...
    '15_01_13_17_13','15_01_14_14_52','15_01_16_13_04','15_01_16_14_06',...
    '15_01_16_15_06','15_01_20_14_42','15_01_20_15_47','15_01_21_16_50',...
    '15_01_22_12_25','15_01_22_13_49'};

all_fix_nov = NaN(length(data_files),50);
all_fix_rep = NaN(length(data_files),50);


for file = 1:length(data_files)
    load([eyedata_dir data_files{file} '-fixation.mat'])
    
    fix_nov = NaN(16,50);
    fix_rep = NaN(16,50);
    for img = 1:size(pairings,2)
        
        %         nov_fixations = fixationstats{pairings(1,img)}.fixations;
        nov_fixationtimes = fixationstats{pairings(1,img)}.fixationtimes;
        
        if size(nov_fixationtimes,2) < 3
            continue;
        end
        
        nov_fixationtimes = nov_fixationtimes(:,2:end);
        
        %         %remove 1st and 2nd fixation if they're within fixation window
        %         if (nov_fixations(1,1) < 100 || nov_fixations(1,1) > -100) && ...
        %                 (nov_fixations(2,1) < 100 || nov_fixations(2,1) > -100)
        %             %central fixation
        %             nov_fixationtimes(:,1) = [];
        %         end
        %
        %         if (nov_fixations(1,1) < 100 || nov_fixations(1,1) > -100) && ...
        %                 (nov_fixations(2,1) < 100 || nov_fixations(2,1) > -100)
        %             %central fixation
        %             dist12 = sqrt(sum((nov_fixations(:,2)-nov_fixations(:,1)).^2));
        %             if dist12 < 25
        %                 nov_fixationtimes(:,1) = [];
        %             end
        %         end
        
        
        %         rep_fixations = fixationstats{pairings(2,img)}.fixations;
        rep_fixationtimes = fixationstats{pairings(2,img)}.fixationtimes;
        
        if size(rep_fixationtimes,2) < 3
            continue;
        end
        
        rep_fixationtimes = rep_fixationtimes(:,2:end);
        
        %         if (rep_fixations(1,1) < 100 || rep_fixations(1,1) > -100) && ...
        %                 (rep_fixations(2,1) < 100 || rep_fixations(2,1) > -100)
        %             %central fixation
        %             rep_fixationtimes(:,1) = [];
        %             rep_fixations(:,1) = [];
        %         end
        %
        %         if (rep_fixations(1,1) < 100 || rep_fixations(1,1) > -100) && ...
        %                 (rep_fixations(2,1) < 100 || rep_fixations(2,1) > -100)
        %             %central fixation
        %             dist12 = sqrt(sum((rep_fixations(:,2)-rep_fixations(:,1)).^2));
        %             if dist12 < 25
        %                 rep_fixationtimes(:,1) = [];
        %             end
        %         end
        
        nov_durs = nov_fixationtimes(2,:)-nov_fixationtimes(1,:)+1;
        if length(nov_durs) > 50, nov_durs = nov_durs(1:50);end
        rep_durs = rep_fixationtimes(2,:)-rep_fixationtimes(1,:)+1;
        if length(rep_durs) > 50, rep_durs = rep_durs(1:50);end
        
        fix_nov(img,1:length(nov_durs)) = nov_durs;
        fix_rep(img,1:length(rep_durs)) = rep_durs;
    end
    
    median_nov_fix = sum(~isnan(fix_nov'));
    median_nov_fix(median_nov_fix == 0) = [];
    median_nov_fix = median(median_nov_fix);
    
    median_rep_fix = sum(~isnan(fix_rep'));
    median_rep_fix(median_rep_fix == 0) = [];
    median_rep_fix = median(median_rep_fix);
    
    fix_nov(:,median_nov_fix+1:end) = [];
    fix_rep(:,median_rep_fix+1:end) = [];
    
    all_fix_nov(file,1:size(fix_nov,2)) = nanmean(fix_nov);
    all_fix_rep(file,1:size(fix_rep,2)) = nanmean(fix_rep);
end

median_nov_fix = sum(~isnan(all_fix_nov'));
median_nov_fix(median_nov_fix == 0) = [];
median_nov_fix = median(median_nov_fix);

median_rep_fix = sum(~isnan(all_fix_rep'));
median_rep_fix(median_rep_fix == 0) = [];
median_rep_fix = median(median_rep_fix);

all_fix_nov(:,median_nov_fix+1:end) = [];
all_fix_rep(:,median_rep_fix+1:end) = [];

figure
hold on
errorbar(nanmean(all_fix_nov),nanstd(all_fix_nov)./sqrt(sum(~isnan(all_fix_nov))))
errorbar(nanmean(all_fix_rep),nanstd(all_fix_rep)./sqrt(sum(~isnan(all_fix_rep))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration')
legend('Novel','Repeat')
title('JN data for VR Calibration Pics')
%%
[~,p] = ttest2(nanmean(all_fix_nov'),nanmean(all_fix_rep'))

figure
hold on
bar([1 2],[mean(nanmean(all_fix_nov')) mean(nanmean(all_fix_rep'))])
errorb(1,mean(nanmean(all_fix_nov')),std(nanmean(all_fix_nov'))./sqrt(length(data_files)))
errorb(2,mean(nanmean(all_fix_rep')),std(nanmean(all_fix_rep'))./sqrt(length(data_files)),'r')
if p < 0.05
    plot(1.5,375,'*k')
end
hold off
set(gca,'Xtick',[1:2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Fixation Duration')