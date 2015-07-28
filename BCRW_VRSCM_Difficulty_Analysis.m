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
%---[2] Run simulations for manipulated images---%
image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VR Image Data\VRSCM Images\';
tags = {'MP','TT','JN','IW'};
Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets'...
    '\CombinedViewingBehavior.mat'];

load(Combinedbehaviorfile,'allview')
imageX = 800; imageY = 600;
plotoptions.runs = 'none'; %all/none
plotoptions.probdens = 'none';
plotoptions.type = 'sal'; %sal/image
IOR_tau = [1/17];
% %same name for everyset
%
for SET = 1:8
   emailme(['Running Simulations for Set ' num2str(SET)])
    cd([image_dir 'VRset' num2str(SET) '\']);
    a = what;
    saliencemapfiles = a.mat; 

    for i = 1:size(saliencemapfiles,1)
        for t = 1:length(tags)
            disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' num2str(SET)])
            run_BCRWCF(allview{t},saliencemapfiles{i},tags{t},imageX,imageY,plotoptions,IOR_tau)
        end
    end
end
