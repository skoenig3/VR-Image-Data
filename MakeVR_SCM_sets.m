%make VR image sets


original_loc = 'C:\Users\seth.koenig\Desktop\Images\';
setnum = 3;
folder_loc = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\';
set_folder = ['VRset' num2str(setnum)];
img_tag = ['VRS' num2str(setnum) 'I'];

new_set_folder = [folder_loc set_folder '\'];
mkdir(new_set_folder);

image_list = ls(original_loc);
%%
%find novel images
novel_images_index = NaN(1,36);
img_num = 1;
for l = 1:size(image_list);
    period = strfind(image_list(l,:),'.bmp');
    if ~isempty(period) %so it's an image
        if ~isnan(str2double(image_list(l,period-1)))
            novel_images_index(img_num) = l;
            img_num = img_num +1;
        end
    end
end

if any(isnan(novel_images_index))
    error('Not enought novel images present')
end
%%
image_type2 = NaN(1,36);

for img = 1:36
    period = strfind(image_list(novel_images_index(img),:),'.bmp');
    name = image_list(novel_images_index(img),1:period-1);
    if exist([original_loc name 'p.bmp'])
        image_type2(img) = 2;
    elseif exist([original_loc name 'r.bmp'])
        image_type2(img)=  3;
    elseif exist([original_loc name 'm.bmp'])
        image_type2(img)= 4;
    else
        error('2nd image not found')
    end
end

order = randperm(36);
novel_images_index = novel_images_index(order);
image_type2 = image_type2(order);

order = randperm(12);
repeat_order = [ones(1,6) 2*ones(1,6)];
repeat_order = repeat_order(order);

order = randperm(12);
replaced_order = [ones(1,6) 2*ones(1,6)];
replaced_order = replaced_order(order);

order = randperm(12);
moved_order = [ones(1,6) 2*ones(1,6)];
moved_order = replaced_order(order);

repeat_count = 1;
replaced_count = 1;
moved_count = 1;

for img = 1:36
    period = strfind(image_list(novel_images_index(img),:),'.bmp');
    name = image_list(novel_images_index(img),1:period-1);
    if image_type2(img) == 2
        if repeat_order(repeat_count) == 1
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name 'p.bmp'],[new_set_folder img_tag num2str(img) 'p.bmp'])
        else
            copyfile([original_loc name 'p.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) 'p.bmp'])
        end
        repeat_count = repeat_count + 1; 
    elseif  image_type2(img) == 3
        if replaced_order(replaced_count) == 1
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name 'r.bmp'],[new_set_folder img_tag num2str(img) 'r.bmp'])
        else
            copyfile([original_loc name 'r.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) 'r.bmp'])
        end
         replaced_count = replaced_count + 1; 
    elseif  image_type2(img) == 4
        if moved_order(moved_count) == 1
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name 'm.bmp'],[new_set_folder img_tag num2str(img) 'm.bmp'])
        else
            copyfile([original_loc name 'm.bmp'],[new_set_folder img_tag num2str(img) '.bmp'])
            copyfile([original_loc name '.bmp'],[new_set_folder img_tag num2str(img) 'm.bmp'])
        end
         moved_count = moved_count + 1; 
    end
end
