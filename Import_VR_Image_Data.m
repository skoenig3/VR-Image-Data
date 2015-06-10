%modified from GoGo_import_eye_cal.m by Seth Konig 5/19/15



data_dir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\eye_calibration\';
img_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\ListSQ\Image Sets\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Figures\';
eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\Eye Data\';

data_files = {'15_01_06_14_28','15_01_06_15_28','15_01_06_16_34',...
    '15_01_07_13_15','15_01_07_14_07','15_01_08_14_40','15_01_09_12_23',...
    '15_01_09_13_36','15_01_12_12_26','15_01_12_13_22','15_01_13_15_25',...
    '15_01_13_17_13','15_01_14_14_52','15_01_16_13_04','15_01_16_14_06',...
    '15_01_16_15_06','15_01_20_14_42','15_01_20_15_47','15_01_21_16_50',...
    '15_01_22_12_25','15_01_22_13_49'}; 

for file = 19:length(data_files)
    eye_cal_file = ['eye_cal2_' data_files{file}]; %contains time, x, and y data
    time_cal_file =['time_cal2_' data_files{file}]; %contains event and tieyeing data
    
    
    %---Read in Eye and Event Data--%
    [eye1 eye2 eye3] = textread([data_dir eye_cal_file],'%f%f%f','headerlines',1,'delimiter',',','whitespace','\n');
    
    try %sometimes eye data file gets cutoff at end?
        eye = [eye1 eye2 eye3];
    catch
        eye = [eye1(1:end-1) eye2 eye3];
    end
    
    t0 = eye(1,1); %time 0
    eye(:,1) = eye(:,1)-t0; %make time relative to 0 otherwise relative to 1970
    
    time = NaN(1,1000); %prealocate space for time of events
    square_pos = NaN(1000,2); %prealocate space for square_position of calibration squares
    events = NaN(1,1000); %prealocate space for events
    % 1: square square_position, 2: square on, 3:reward, 4:break fixation
    % 5: clrchng, 6: square off
    
    fid = fopen([data_dir time_cal_file]);
    linecount = 1;
    tline = fgets(fid); %skip line 1: timestaeyep, events
    tline = fgets(fid); %skip line 2: [timestaeyep1] start collecting eye data
    image_name = {};
    img_count = 1;
    while ischar(tline)
        tline = fgets(fid); %get next line
        if tline ~= -1; %end of file is noted by -1
            C = textscan(tline,'%s'); %parse line by spaces and comas
            time(linecount) = str2double(C{1}{1});
            if strcmpi(C{1}{2},'Square')
                if strcmpi(C{1}{3},'Position,') %gives position of square in pixels
                    events(linecount) = 1;
                    square_pos(linecount,1) = str2double(C{1}{4});
                    square_pos(linecount,2) = str2double(C{1}{5});
                elseif strcmpi(C{1}{3},'moved') || strcmpi(C{1}{3},'on') %square moved to different location and displayed
                    events(linecount) = 2;
                elseif strcmpi(C{1}{3},'dims') %color change occured
                    events(linecount) = 5;
                elseif strcmpi(C{1}{3},'off') %suqre diseappears
                    events(linecount) = 6;
                else
                    disp('unknown events parameter in column 3 of event array') %so we known we didn't miss anything
                end
            elseif strcmpi(C{1}{2},'reward') %monkey received reward
                events(linecount) = 3;
            elseif strcmpi(C{1}{2},'no') %trial aborted due to break fixation, etc.
                events(linecount) = 4;
            elseif ~isempty(strfind(lower(C{1}{2}),'photo')) %image trials
                if length(C{1}) == 3
                    if  strcmpi(C{1}{3},'on') %image on
                        events(linecount) = 7;
                    elseif strcmpi(C{1}{3},'off') %image off
                        events(linecount) = 8;
                    end
                else
                    if~isempty(strfind(C{1}{2},'/')) %image name
                        events(linecount) = 9;
                        slash = strfind(C{1}{2},'/');
                        image_name{img_count} = C{1}{2}(slash+1:end);
                        img_count = img_count+1;
                    else
                        disp('Unknown parameter for column 2 for picture trials')  %so we known we didn't miss anything
                    end
                end
            end
            linecount = linecount+1; %to go to the next line
        end
    end
    fclose(fid);
    
    %---Get Set number for images...for later use---%
    setnum = str2double(image_name{1}(2:3)); %get the set number, assumes set number is the same for all images,
    % may be incorrect
    set_dir = [figure_dir 'Set' num2str(setnum) '\'];
    mkdir(set_dir);  %make directory if it doesn't already
    %exist for individual sets for organization purpose
    
    %---Get the calibration Function---% 
    %doing for individual sessions of calibration. Currently assumes stable
    %cailbration within short session. Does not track drift
    if linecount < 1000 %if you didn't use the whole prealocated space remove extra rows
        time = time(1:linecount-1);
        events = events(1:linecount-1);
        square_pos = square_pos(1:linecount-1,:);
    end
    time = time-t0; %use t0 to ensure time is the saeyee but should be
    
    % get unique calibration points
    xsquare_poses = unique(square_pos(:,1)); xsquare_poses(isnan(xsquare_poses)) = [];
    ysquare_poses = unique(square_pos(:,2)); ysquare_poses(isnan(ysquare_poses)) = [];
    
    
    rewards = find(events == 3); %find where a reward occured
    caldat_x = cell(length(xsquare_poses),length(ysquare_poses));
    caldat_y = cell(length(xsquare_poses),length(ysquare_poses));
    square_dims = find(events == 5);
    square_pos_events = find(events == 1);
    for rew = 1:length(rewards);
        rewardind = rewards(rew); %index of event array in which reward was delivered
        
        %trialstarttieye = time(rewardind-3); %square on
        square_change = square_dims(square_dims < rewardind); %find the last time the color change
        square_change = square_change(end);
        clrchngtieye = time(square_change); %look at time in which square changes color (dims)
        start_data_collection = clrchngtieye - 0.5; %500 eyes
        
        %find indices associated with times between 500 ms before color change
        %and color change
        eyeind = find(eye(:,1) > start_data_collection & eye(:,1) <= clrchngtieye);
        
        %determine x and y eye position for those 500 ms
        xsquare_pos = eye(eyeind,2);
        ysquare_pos = eye(eyeind,3);
        
        %position of square
        sq_pos = square_pos_events(square_pos_events < rewardind);
        sq_pos = sq_pos(end);
        square_x = square_pos(sq_pos,1);
        square_y = square_pos(sq_pos,2);
        
        xind = find(xsquare_poses == square_x);
        yind = find(ysquare_poses == square_y);
        
        caldat_x{xind,yind} = [caldat_x{xind,yind} mean(xsquare_pos)];
        caldat_y{xind,yind} = [caldat_y{xind,yind} mean(ysquare_pos)];
    end
    
    % control are unique combinations x and y position of square in pixel coordinates
    control_x = [];
    control_y = [];
    for xi = 1:length(xsquare_poses);
        for yi = 1:length(ysquare_poses);
            control_x = [control_x; xsquare_poses(xi)];
            control_y = [control_y; ysquare_poses(yi)];
        end
    end
    
    % input are estimated position of the color change squares in mV (the ouptu
    % of eyescan).
    input_x = [];
    input_y = [];
    for xi = 1:length(xsquare_poses);
        for yi = 1:length(ysquare_poses);
            input_x = [input_x; mean(caldat_x{xi,yi})];
            input_y = [input_y; mean(caldat_y{xi,yi})];
        end
    end
    
    %tform automatically transforms estimated eye positions in mV into pixel
    %coordinates in  a "magical" way.
    tform = cp2tform([control_x control_y], [input_x input_y],'polynomial',4);
    tform.forward_fcn = tform.inverse_fcn;
    
    
    %---Create Plot of Calibration Quality---%
    %Visualize position of estimated square position (blue *) compared to
    %actual position of color change squares (red +)
    [xp,yp] = tformfwd(tform,input_x,input_y);%transform inputs from mV to pixels to plot
    figure
    hold on
    for xi = 1:length(xsquare_poses);
        for yi = 1:length(ysquare_poses);
            plot(xsquare_poses(xi),ysquare_poses(yi),'r+')
        end
    end
    for xi = 1:length(xp);
        plot(xp(xi),yp(xi),'b*')
    end
    hold off
    save_and_close_fig(set_dir,['Calibration for ' data_files{file}])
   
    
    %---Recalibrate Collected Eye Data---%
    eyedata = {};
    img_on = find(events == 7);
    img_off = find(events == 8);
    for img = 1:length(image_name);
        eye_ind = find(eye(:,1) > time(img_on(img)) & eye(:,1) <= time(img_off(img)));
        
        [xp,yp] = tformfwd(tform,eye(eye_ind,2)',eye(eye_ind,3)');%transform inputs from mV to pixels to plot
        eyedat{img} = [xp;yp];
    end
    
    %---Get image numbers and the order they appeared in---%
    imgnum = NaN(1,32);
    for img = 1:length(image_name);
        imgnum(img) = str2double(image_name{img}(5:6));
    end
    
    pairings = NaN(2,16);
    for img = min(imgnum):max(imgnum) %assumes images don't cross over sets
        ind = find(imgnum == img);
        pairings(1,img-min(imgnum)+1) = ind(1); %novel
        pairings(2,img-min(imgnum)+1) = ind(2); %repeat
    end
    
    %---Plot Eye Data on Corresponding Figure---%
    for img = 1:size(pairings,2);
        im1 = imread([img_dir 'LSQ' image_name{pairings(1,img)}(2:3)...
            '\' image_name{pairings(1,img)}]);
        im2 = imread([img_dir 'LSQ' image_name{pairings(2,img)}(2:3)...
            '\' image_name{pairings(2,img)}]);
        %im1 should be the same as im2 but want to make sure that
        %everythign works so going to load 2x
        
        figure
        subplot(1,2,1)
        imshow(im1)
        hold on
        plot(eyedat{pairings(1,img)}(1,:)+400,...
            600-(eyedat{pairings(1,img)}(2,:)+300))
        hold off
        title('Novel Image')
        
        subplot(1,2,2)
        imshow(im2)
        hold on
        plot(eyedat{pairings(2,img)}(1,:)+400,...
            600-(eyedat{pairings(2,img)}(2,:)+300))
        hold off
        title('Repeat Image')
        
        subtitle(['Set ' num2str(setnum) ' Image ' num2str(pairings(1,img)) ...
            ' Spacing ' num2str(pairings(2,img)-pairings(1,img))]);
        
        save_and_close_fig(set_dir,['Set_' num2str(setnum) '_img_' ...
            num2str(pairings(1,img)) '_' data_files{file}]);
    end
    
    fixationstats = ClusterFixation_Short(eyedat);
    save([eyedata_dir data_files{file} '-fixation.mat'],'fixationstats',...
        'image_name','pairings','imgnum','setnum')
end