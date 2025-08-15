function [data_reord]=region_reordering_AW(goal, atlas, inp)
% function reorderes HCP data according to order of regions as in CATO DWI 
% output; Alina Podschun and Sebastian Markett (2022)
% adapted by Aradia Wilms for HCP MEG data
%
% needs variables:
%   - goal -> this has the regions as named in CATO DWI output
%   - atlas -> ciftified HCP parcellation
%   - inp -> match-ordered input (cope file) that needs to be reordered to goal-order
%
% output:
%   - vector containing CATO-ordered region label and corresponding value
%       of functional data of cope file of interest

%% create list of labels

% test lists all region names that contain the ctx-lh prefix (instead of
% rh)
test = goal(contains(goal,'ctx-lh')==1);

% check which part is the individual names; this is about finding out where
% the prefix ends, in this case through trial we know its 8; change number
% if needed
test{1}(8:numel(test{1}));

% make boolean list lh which contains 1 for every left hemisphere region
lh = contains(goal,'ctx-lh');
rh = contains(goal,'ctx-rh');

% add L and R to pure label name extractions; result is vector "label" that
% contains labels for all 219 regions -> names are freesurfer style but
% order is CATO style still
for i=1:length(goal)
    if lh(i) == 1
        label{i,1} = ['L' goal{i}(8:numel(goal{i}))];
    elseif rh(i) == 1
        label{i,1} = ['R' goal{i}(8:numel(goal{i}))];
    else 
        % mark subcortical areas with an X so they can be deleted
        % the CATO Desikan Killiany contains subcortical ROIs as well
        label{i,1} = ['X'];
    end
    label{i,2} = i;
end


% delete subcortical areas 
rows_to_delete = strcmp(label(:, 1), 'X');
label(rows_to_delete, :) = [];



% now the label template has 68 regions!


%% build matching HCP ordered list and inp

% extract relevant content from parc

% result here is vector "match" with freesurfer-style names in
% freesurfer-style order
inmed = {atlas.diminfo{1,2}.maps.table.name};
for k = 1:71
    match{k,1} = inmed{1,k};
end


%%
% from match, drop unknown area and corpus callosum
match=match(~ismember(match(:,1),'???'),:);
match=match(~ismember(match(:,1),'Lcorpuscallosum'),:);
match=match(~ismember(match(:,1),'Rcorpuscallosum'),:);


%%
% delete corpus callosum from functional input
inp(4,:)=[];

% note: this index needs to be row -1 because the table will shift from
% removing row 4
inp(38,:)=[];

%%  reorder functional data in a reordered way (matching CATO DWI output)
data_reord = zeros(68,1);
for i=1:68
    inter=label{i};
    find = strcmp(inter,match);
    for k=1:68
        if find(k)==1
            data_reord(i)=inp(k);
        else
            continue
        end
    end
end
