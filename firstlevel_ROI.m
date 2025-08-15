function [mean_firstlevel]=firstlevel_ROI_AW(atlas, cii)
%%
% This function carries out the extraction of first-level ROIs of
% functional data

% based on scripts provided by Alina Podschun

% functional data is stored in cii
% dlabel.nii file for ROI info is stored in atlas

%%


try
    cii = cii.cdata; %extract dense timeseries

    for k=1:max(atlas.cdata) % loop across ROIs
            sysmat(k,1)=mean(mean(cii(atlas.cdata==k,:))); %extract mean activation
            sysmat(k,2)=k;
    end
        mean_firstlevel = sysmat; % assign mat
        clear sysmat % not needed  any more
        clear cii % not needed  any more
catch ME
    disp(['Failed: ', ME.message])
end
    
