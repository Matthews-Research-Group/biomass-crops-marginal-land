function subplots_one_multiple_colormaps2(dataset,fig_name,cmins,cmaxs,cmaps,max_index,ml,titles,unit)
screensize = get(0, 'Screensize');
customsize = screensize;
customsize(3) = customsize(3)/2;%width
customsize(4) = customsize(4);%height
load('conus_masks/conus_mask.mat') %#ok<LOAD>
FigH = figure('Position', customsize);
fs = 20;
lw = 2;
nrows = 1;ncols=1;
% xlimit = [-2500,2200];
xlimit = [-1100,2200];
ylimit = [-1600,1600];
xy_km = csvread('CWRF_coor_in_km.csv',1,0);
x=xy_km(:,1);
y=xy_km(1:138,2);
total_figs = 3;
remove_western=true;

for i=1:total_figs
    ax{i}=axes('xlim',xlimit,'ylim',ylimit,'Color','none');
    plot(gca,2000,1000,'wo')
    hold on;
    if(length(cmaps)>1)
        my_cmap = load(cmaps{i});
        my_cmap = my_cmap/255;
    else
        if(strcmp(cmaps{1},'my_cmap'))
            my_cmap = load('colormaps/BlueRed10.txt');
        else
            my_cmap = cmaps{1};
        end
    end
    if length(cmins)>1 
        cmin = cmins(i);
        cmax = cmaxs(i);
    else
        cmin = cmins;
        cmax = cmaxs;
    end
    if total_figs>1 
     z = dataset(:,:,i);
    else
     z = dataset;
    end
    z(z<cmin) = cmin;
    z(z>cmax) = cmax;
    z( conus < 0.5) = 9999; %Mask out the non-CONUS area
    z( ml < 0.05 ) = 9999; %mask out the non-ML grids.
    %keep the current crop type grid points only
    z(max_index ~= i) = 9999;
    
    z = z'; %REQUIRE a Transpose here!!!!!!
    non_blanks = z~=9999 & ~isnan(z);
    im=imagesc(ax{i},x,y,z,'AlphaData', non_blanks);
%         imagesc(x,y,z,'AlphaData', zeros(size(z)))
%         s =  pcolor(x,y,z);
%         s.EdgeColor = 'none';
    set(gca,'fontsize',fs,'linewidth',0.1,'xtick',[],'ytick',[])
%          set(gca, 'visible', 'off');
    % Set colormap and color limits for all subplots
%         set(ax, 'Colormap', jet, 'CLim', [cmin cmax])
    set(gca, 'Colormap', my_cmap, 'CLim', [cmin cmax])
    %         title(titles{i},'FontSize',16)
    %Add some texts at the bot-left corner
    pos = get(gca, 'position');
%     dim = pos.*[1 1 0.2 0.2];
%     annotation(gcf, 'textbox', dim, 'String', titles,'FontSize',14,'LineStyle','none')
    %add colorbar at the end
%             cb=colorbar('southoutside','XTick', cmin:cmax);
    cb=colorbar('southoutside');
    %[left bottom width height]
    %deal with position and width of colorbars
    if(length(cmaps)>1) %only needed if we have multiple color maps
        left_indent = 0.05;
        width_old = cb.Position(3);
        cb.Position(3) = width_old*1/3 - 0.1;
        cb.Position(2) = cb.Position(2);% + 0.1;
        cb.Position(1) = cb.Position(1) + left_indent;
        if(i==2)
            cb.Position(1) = cb.Position(1)+width_old/3;
        elseif(i==3)
    %         disp(cb.Position)
            cb.Position(1) = cb.Position(1)+width_old*2/3;
    %             annotation(gcf,'textbox',...
    %             [cb.Position(1)+cb.Position(3)+0.01 0.045 0.018 0.053],...
    %             'LineStyle','none',...
    %             'FontSize',14,...
    %             'String',{unit});
        end
    else
        left_indent = 0.1;
        width_old = cb.Position(3);
        cb.Position(3) = width_old - left_indent*2;
        cb.Position(2) = cb.Position(2);% + 0.1;
        cb.Position(1) = cb.Position(1) + left_indent;
    end
    hold off
    % Add CONUS Map
    add_conus_map(ax{i},xlimit,ylimit,remove_western)
end
    linkaxes([ax{1}, ax{2}, ax{3}])
    ax{3}.Visible = 'off';
    ax{3}.XTick = [];
    ax{3}.YTick = [];
    ax{2}.Visible = 'off';
    ax{2}.XTick = [];
    ax{2}.YTick = [];
    print(fig_name,'-dpng','-r300')
    % close;
end

function add_conus_map(current_ax,xlimit,ylimit,remove_western)
        us=shaperead('./USshp/us_lcc_YH/us_lcc.shp');

        if(remove_western)
            states_midwest  = ["North Dakota","South Dakota","Nebraska","Kansas","Minnesota","Iowa","Missouri",...
               "Wisconsin","Michigan","Illinois","Indiana","Ohio","Pennsylvania","West Virginia","Virginia",...
               "Maryland","New Jersey","New York","New Jersey","Vermont","New Hampshire","Maine",...
               "Massachusetts","Connecticut","Rhode Island","Delaware"]; 
            states_southern  = ["Texas","Oklahoma","Arkansas","Louisiana",...           
                       "Mississippi","Kentucky","Tennessee","Alabama","Florida",...
                       "North Carolina","South Carolina","Georgia"];
            states_all = extractfield(us,'NAME');
            states_all = string(states_all);

            ind = ismember(states_all,[states_midwest,states_southern]);
            ind_out = ismember(states_all,[states_midwest,states_southern])==0;
            us_sub_keep = us(ind);
            us_sub_remove = us(ind_out);

            map=mapshow(current_ax,us_sub_keep);
    %         set(map.Children,'FaceAlpha',0);
            set(map.Children,'FaceColor','none');
            set(map.Children,'EdgeColor',[0.85,0.85,0.85]);
            set(map.Children,'EdgeAlpha',1);
            set(map.Children,'LineWidth',1);
            xlim(xlimit)
            ylim(ylimit)

            hold on;
            map=mapshow(current_ax,us_sub_remove);
    %         set(map.Children,'FaceAlpha',0);
            set(map.Children,'FaceColor','white');
            set(map.Children,'EdgeColor',[0.85,0.85,0.85]);
            set(map.Children,'EdgeAlpha',1);
            set(map.Children,'LineWidth',1);
            xlim(xlimit)
            ylim(ylimit)
        else
            map=mapshow(current_ax,us);
            set(map.Children,'FaceColor','none');
            set(map.Children,'EdgeColor',[0.85,0.85,0.85]);
            set(map.Children,'EdgeAlpha',1);
            set(map.Children,'LineWidth',1);
            xlim(xlimit)
            ylim(ylimit)
        end

        hold on;
        conus_clean=shaperead('USshp/conus_shp_clean_lcc/conus.shp');
        poly1 = polyshape(conus_clean.X,conus_clean.Y);
%         plot(poly1);
%         map=mapshow(conus_clean);
%         set(map.Children,'FaceColor','none');
%         set(map.Children,'EdgeColor','b');
%         set(map.Children,'EdgeAlpha',1);
%         set(map.Children,'LineWidth',1);
%         xlim(xlimit)
%         ylim(ylimit)
%        set(gca,'Color','black');%Make all color outside CONUS white, including lakes
%         setm(gca,'ffacecolor','white')
        X=poly1.Vertices;
        fillout(X(:,1),X(:,2),[xlimit,ylimit],'white');
end