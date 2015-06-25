%remove images already organized into sets
VRsets = 1:4;

VRset_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\VRset';

values = NaN(72*4,3);
count = 1;
for v = VRsets
    cd([VRset_dir num2str(VRsets(v)) '\'])
    list = ls;
    for l = 1:size(list,1)
        if ~isempty(strfind(list(l,:),'.bmp'))
            img = imread(list(l,:));
            img = rgb2gray(img);
            values(count,1) = mean2(img);
            values(count,2) = std2(img);
            values(count,3)=  entropy(img);%pixel intesnity entropy
            count = count+1;
        end
    end
end

cd('C:\Users\seth.koenig\Desktop\all repeats\');
list = ls;

for l = 1:size(list,1)
    if ~isempty(strfind(list(l,:),'.bmp'))
        img = imread(list(l,:));
        img = rgb2gray(img);
        values2(1) = mean2(img);
        values2(2) = std2(img);
        values2(3)=  entropy(img);%pixel intesnity entropy
        
        D = pdist2(values,values2);
        if any(D == 0);
            delete(list(l,:))
            disp('img deleted')
        end
    end
end

cd('C:\Users\seth.koenig\Desktop\all replaced\');
list = ls;

for l = 1:size(list,1)
    if ~isempty(strfind(list(l,:),'.bmp'))
        img = imread(list(l,:));
        img = rgb2gray(img);
        values2(1) = mean2(img);
        values2(2) = std2(img);
        values2(3)=  entropy(img);%pixel intesnity entropy
        
        D = pdist2(values,values2);
        if any(D == 0);
            delete(list(l,:))
            disp('img deleted')
        end
    end
end


cd('C:\Users\seth.koenig\Desktop\all moved\');
list = ls;

for l = 1:size(list,1)
    if ~isempty(strfind(list(l,:),'.bmp'))
        img = imread(list(l,:));
        img = rgb2gray(img);
        values2(1) = mean2(img);
        values2(2) = std2(img);
        values2(3)=  entropy(img);%pixel intesnity entropy
        
        D = pdist2(values,values2);
        if any(D == 0);
            delete(list(l,:))
            disp('img deleted')
        end
    end
end

