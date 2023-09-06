clear;clc;close all;
plot_options = 1;%1:6
%1: hist;      2:futu;       3: bias(futu-hist)
%4: hist-std;  5:futu-std;   6: bias-std;
frostdoy = ncread('multicrops/xfrostdoy.nc','FROSTDOY');
yield_hist_ncfile = 'multicrops/multi_crop_yield_mpi_2000_2014_r7.nc';
yield_futu_ncfile = 'multicrops/multi_crop_yield_mpi_2036_2050_r7.nc';

hardiness_hist_ncfile = 'hardiness/hardiness_CWRF_2000_2014_zones_5a_8a.nc';
hardiness_futu_ncfile = 'hardiness/hardiness_CWRF_2036_2050_zones_5a_8a.nc';
%these two files do not matter here. It's just for reading in ONLY!!
crop_allocation_hist = 'crop_allocation_index/max_index_yield_based_hist_r3.mat';
crop_allocation_futu = 'crop_allocation_index/max_index_yield_based_futu_r3.mat';
%
for k = 1:length(plot_options)
    plot_option = plot_options(k);
    if plot_option==1
        yield_all = ncread(yield_hist_ncfile,'yield');
        wue_all   = ncread(yield_hist_ncfile,'wue');
        hardiness = ncread(hardiness_hist_ncfile,'hardiness');
        miscanthus = yield_all(:,:,1);
        miscanthus(hardiness~=1) = nan; %1 is miscanthus
        yield_all(:,:,1) = miscanthus;
        loadfile  = load(crop_allocation_hist);
        max_index = loadfile.max_index1;
        outfile_namekey = 'hist';
        cmins = [0,0,0];
        cmaxs = [50,20,50];
        cmaps = {'colormaps/Blues_6.gpl','colormaps/Greens_6.gpl','colormaps/RdPu_6.gpl'};
    elseif plot_option==2
        yield_all = ncread(yield_futu_ncfile,'yield');
        wue_all   = ncread(yield_futu_ncfile,'wue');
        hardiness = ncread(hardiness_futu_ncfile,'hardiness');
        miscanthus = yield_all(:,:,1);
        miscanthus(hardiness~=1) = nan; %1 is miscanthus
        yield_all(:,:,1) = miscanthus;
        loadfile = load(crop_allocation_futu);
        max_index = loadfile.max_index1;
        outfile_namekey = 'futu';
        cmins = [0,0,10];
        cmaxs = [30,10,60];
        cmaps = {'colormaps/Blues_6.gpl','colormaps/Greens_6.gpl','colormaps/RdPu_6.gpl'};
    elseif plot_option==3
        outfile_namekey = 'bias';
        yield_hist = ncread(yield_hist_ncfile,'yield');
        hardiness_hist = ncread(hardiness_hist_ncfile,'hardiness');
        yield_futu = ncread(yield_futu_ncfile,'yield');
        hardiness_futu = ncread(hardiness_futu_ncfile,'hardiness');
        %this difference is the difference between the same specie!!!
        yield_diff = yield_futu - yield_hist;        
        loadfile  = load(crop_allocation_hist);
        max_index_hist = loadfile.max_index1;
        loadfile  = load(crop_allocation_futu);
        max_index_futu = loadfile.max_index1;
        %compare futu and hist Index
        max_index_diff = max_index_futu - max_index_hist;
        %find locations that are changed in crop types
        changed_points = find(max_index_diff~=0 & ~isnan(max_index_diff));
        
        yield_hist_1d=zeros(length(changed_points),3);
        yield_futu_1d=zeros(length(changed_points),3);
        for i=1:3
            tmp = yield_hist(:,:,i);
            yield_hist_1d(:,i)=tmp(changed_points);
            tmp = yield_futu(:,:,i);
            yield_futu_1d(:,i)=tmp(changed_points);
        end
        %this is the actual difference between species 
        yield_diff_actual = yield_diff;
        yield_changed_hist = indexing_yield(yield_hist_1d,max_index_hist(changed_points));
        yield_changed_futu = indexing_yield(yield_futu_1d,max_index_futu(changed_points));
        for i=1:3
            tmp = yield_diff_actual(:,:,i);
            tmp(changed_points) = yield_changed_futu - yield_changed_hist;
            yield_diff_actual(:,:,i) = tmp;
        end
        yield_all = yield_diff_actual;
        max_index = max_index_futu;
%         %assign a number to these points to remove.
%         %show only the same specie.
%         max_index(changed_points) = 99; 
        cmins = -10;
        cmaxs = 10;
        cmaps = {'my_cmap'};
    elseif plot_option==4
        yield_std = ncread(yield_hist_ncfile,'yield_std');
        yield = ncread(yield_hist_ncfile,'yield');
        yield_all = yield_std./yield*100;
        hardiness = ncread(hardiness_hist_ncfile,'hardiness');
        loadfile  = load(crop_allocation_hist);
        max_index = loadfile.max_index1;
        outfile_namekey = 'hist_std_percent';
        cmins = 0; cmaxs = 30;
        cmaps = {parula};
    elseif plot_option==5
        yield_std = ncread(yield_futu_ncfile,'yield_std');
        yield = ncread(yield_futu_ncfile,'yield');
        yield_all = yield_std./yield*100;
        hardiness = ncread(hardiness_futu_ncfile,'hardiness');
        loadfile  = load(crop_allocation_futu);
        max_index = loadfile.max_index1;
        outfile_namekey = 'futu_std_percent';
        cmins = 0; cmaxs = 30;
        cmaps = {parula};
    elseif plot_option==6
        outfile_namekey = 'bias_std_percent';
        yield_std = ncread(yield_hist_ncfile,'yield_std');
        yield = ncread(yield_hist_ncfile,'yield');
        yield_hist = yield_std./yield*100;
        hardiness_hist = ncread(hardiness_hist_ncfile,'hardiness');
       
        yield_std = ncread(yield_futu_ncfile,'yield_std');
        yield = ncread(yield_futu_ncfile,'yield');
        yield_futu = yield_std./yield*100;
    %     wue_all   = ncread('multicrops/multi_crop_yield_mpi_2036_2049.nc','wue');
        hardiness_futu = ncread(hardiness_futu_ncfile,'hardiness');
        
        yield_all = yield_futu - yield_hist;
        loadfile  = load(crop_allocation_hist);
        max_index_hist = loadfile.max_index1;
        loadfile  = load(crop_allocation_futu);
        max_index_futu = loadfile.max_index1;
        %compare futu and hist Index
        max_index_diff = max_index_futu - max_index_hist;
        %find locations that are changed in crop types
        changed_points = find(max_index_diff~=0 & ~isnan(max_index_diff));
        max_index = max_index_futu;
        max_index(changed_points) = 99; %assign a number to these points
        cmins = -10;
        cmaxs = 10;
        cmaps = {'my_cmap'};
    end

% %maximum index
% %this part is used to generate the index map, saved as .mat
    max_index1 = zeros(195,138)*nan; %yield-based index
    max_index2 = zeros(195,138)*nan; %wue-based index
    for i =1:195
        for j=1:138
            values1 = yield_all(i,j,:);
            values2 = wue_all(i,j,:);
            
%             %make switchgrass to be zero to exclude it
%             values1(2) = 0;
            
            [maxV,maxInd]=max(values1,[],'omitnan');
            max_index1(i,j) = maxInd;
            
            if(hardiness(i,j)==2)  %for energycane area, compare switchgrass & energycane
            %here in the hardiness file, 2 is energycane, but it's 3 in the
            %index file
                values2 = wue_all(i,j,2:3);
                [maxV,maxInd]=max(values2,[],'omitnan');
                max_index2(i,j) = maxInd+1;
            else
                [maxV,maxInd]=max(values2,[],'omitnan');
                max_index2(i,j) = maxInd;
            end
        end
    end
    % max_index2(~isnan(hardiness))=nan;
    max_index1(isnan(wue_all(:,:,1))) = nan;
    max_index2(isnan(wue_all(:,:,1))) = nan;
    save(['crop_allocation_index/max_index_yield_based_',outfile_namekey,'_r7_new.mat'],'max_index1')
    save(['crop_allocation_index/max_index_wue_based_',outfile_namekey,'_r7_new.mat'],'max_index2')
end

function y=indexing_yield(A,index)
    [m,n] = size(A);
    B=A';
    index_1d = index+[0:1:(m-1)]'*n;
    y = B(index_1d);
end
