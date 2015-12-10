function CLImATHET_plot_normalized_results(Datafile,resultsfile,plotsdir,barcode)
% 21/05/2015 by Zhenhua
% This function is used to plot figures given the data, all 
% these figures are stored in a specified directory

%----------------------read results  ------------------------%
fid = fopen(resultsfile,'r');
if fid == -1
    error(['Can not open result file: ' resultsfile]);
end

%get estimated parameters from the first row of the result file
beta = [];
lambda = [];

while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'StartPos')),break,end
%     % alpha
%     str = 'Frequencies of subclonal population in the sample: ';
%     indx = strfind(tline, str);
%     if ~isempty(indx)
%         result1 = regexp(tline(indx+length(str):end),'\s','split');
%         if ~isempty(result1)
%             alpha = [alpha str2double(result1)];
%         end
%     end
%     result1 = regexp(tline,'Proportion in the sample:\s*(\S+)','tokens','once');
%     if ~isempty(result1)
%         alpha = [alpha str2double(result1{1})];
%     end
    % beta
    str = 'Subclonal populations in the sample: ';
    indx = strfind(tline, str);
    if ~isempty(indx)
        result1 = regexp(tline(indx+length(str):end),'\s','split');
        if ~isempty(result1)
            beta = [beta str2double(result1)];
        end
    end
%     % p
%     str = 'Parameter p of NB distributions: ';
%     indx = strfind(tline, str);
%     if ~isempty(indx)
%         result1 = regexp(tline(indx+length(str):end),'\s','split');
%         if ~isempty(result1)
%             p = [p str2double(result1)];
%         end
%     end
    % lambda
    result1 = regexp(tline,'Copy neutral read counts:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        lambda = [lambda str2double(result1{1})];
    end
end
%report errors if these values are not parsed successfully
if isempty(beta)
    error(['Can not read estimated subclonal populations from file ',resultsfile]);
end
% if isempty(p)
%     error(['Can not read estimated parameter p of NB distributions from file ',resultsfile]);
% end
if isempty(lambda)
    error(['Can not read estimated read counts for normal copy from file ',resultsfile]);
end

beta = [0 beta];

% then read the results
% '%d\t%d\t%d\t%d\t%1.3f\t%d\t%f\t%3.1f\t%d\n'
results = textscan(fid,'%f %f %f %f %f %*f %f %*f','treatAsEmpty', {'NA', 'na'});
fclose(fid);
chr_seg = results{1};
pstart_seg = results{2};
pend_seg = results{3};
cn_seg = results{4};
mcn_seg = results{5};
% AI_seg = results{6};
purity_seg = results{6};
% score_seg = results{8};
clear results;

eval(['load ' Datafile]);

tv = data_chr_all == 24;
if sum(tv) > 0
    gender = 1;
else
    gender = 0;
end

h=figure(1);
set(h,'visible','off');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 13 8])
FontSize = 15;

%--------------- plot figures ---------------------%
rc_colors = [0.5 0.5 0.5;
             0 0.9 0;
             0 0 0.9;
             0.9 0 0];
baf_colors = [0 250 0;
             0 0 250;
             250 0 0]./255;

for i = reshape(intersect(unique(data_chr_all),1:30),1,[])
    tv = ismember(data_chr_all,i);
    data_rc = data_rc_all(tv);
    data_bd = data_bd_all(tv);
    data_td = data_td_all(tv);
    data_baf = data_bd./data_td;
    data_pos = data_pos_all(tv);
    min_pos = min(data_pos)-1000;
    max_pos = max(data_pos)+1000;
    
    indx1 = find(chr_seg==i);
    
    %plot
    clf;
    marker_size = 3;
    %---plot CN---
    subplot(4,1,1)
    hold on
    set(gca,'YGrid','on')
    set(gca,'FontSize',FontSize);
    %     axis ([-Inf Inf -0.05 7.05])
%     axis ([-Inf Inf -0.1 7.1])
    axis ([min_pos max_pos -0.1 7.1])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:1:7],'Box','on')
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('Copy number');
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        mCN = mcn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN > 7
            CN = 7;
        end
        line_style = 'r-';
        plot([data_pos(indx(1)) data_pos(indx(end))],[CN CN], line_style, 'LineWidth',2.5);
        line_style = 'b-';
        if CN == mCN
            plot([data_pos(indx(1)) data_pos(indx(end))],[mCN-0.25 mCN-0.25], line_style, 'LineWidth',2.5);
        else
            plot([data_pos(indx(1)) data_pos(indx(end))],[mCN mCN], line_style, 'LineWidth',2.5);
        end
    end
    %replace '_' with '-' in barcode to display it correctly
	barcode_m = barcode;
    tmp = strfind(barcode_m,'_');
    barcode_m(tmp) = '-';
    set(gca,'XTick',[])
    title (['Chromosome ' num2str(i) ', ' barcode_m])
    
    %---plot BAF---
    subplot(4,1,2)
    set(gca,'FontSize',FontSize);
    hold on
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        if CN == 0
            Muc = 0.5;
        else
            Muc = mcn_seg(j)/cn_seg(j);
        end
        
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 2
            k = 1; % Del
        else
            if Muc == 1
                k = 3; % LOH
            else
                k = 2; % Het
            end
        end
        plot(data_pos(tv),data_baf(tv),'.','MarkerSize',marker_size, 'Color', baf_colors(k,:));
    end
%     plot(data_pos,data_baf,'b.', 'MarkerSize',marker_size)
    for j = 0:0.25:1
        plot([data_pos(1) data_pos(end)],[j j],'k-','LineWidth',0.5)
    end
    %plot expected BAF mean values
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        mCN = mcn_seg(j);
        w = purity_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN == 0
            baf_mean = 0.5;
        else
            if gender == 1 && (i == 23 || i == 24)
                temp1 = w*mCN+(1-w);
                temp2 = w*CN+(1-w);
            else
                temp1 = w*mCN+(1-w);
                temp2 = w*CN+(1-w)*2;
            end
            baf_mean = temp1/temp2;
        end
        plot([data_pos(indx(1)) data_pos(indx(end))],[baf_mean baf_mean],'k-','LineWidth',1.2);
        plot([data_pos(indx(1)) data_pos(indx(end))],[1-baf_mean 1-baf_mean],'k-','LineWidth',1.2);
    end
    ylabel('BAF');
    axis ([min_pos max_pos -0.03 1.03])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0 0.5 1],'Box','on')
    set(gca,'XTick',[])
    
    %---plot Read count---
    subplot(4,1,3);
    set(gca,'FontSize',FontSize);
    hold on
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 1
            k = 1;
        else
            k = CN+1;
        end
        if k > 4
            k = 4;
        end
        plot(data_pos(tv),data_rc(tv),'.','MarkerSize',marker_size, 'Color', rc_colors(k,:));
    end
    % plot expected RC mean values
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        w = purity_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN == 0 %total deletion
            rc_median = lambda*(1-w);
        else
            if gender == 1 && (i == 23 || i == 24)
                rc_median = lambda*(w*CN+(1-w))/2; 
            else
                rc_median = lambda*(w*CN+(1-w)*2)/2; 
            end
        end
        plot([data_pos(indx(1)) data_pos(indx(end))],[rc_median rc_median],'k-','LineWidth',1.2);
    end
    
    ylabel('Read counts');
%     set(gca,'XTick',[]);
    set(gca,'Box','on')
%     axis ([min_pos max_pos min_rd max_rd])
    axis ([min_pos max_pos -Inf Inf])
    set(gca,'XTick',[])
%     axis([-Inf Inf -Inf Inf])

    %---plot Tumor fraction---
    subplot(4,1,4);
    set(gca,'FontSize',FontSize);
    hold on
    for j=reshape(indx1,1,[])
        line_style = 'r-';
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([data_pos(indx(1)) data_pos(indx(end))],[purity_seg(j) purity_seg(j)], line_style, 'LineWidth',2.5);
    end
    axis ([min_pos max_pos -0.03 1.03])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:0.5:1],'Box','on')
    ylabel('Cellularity');
    
%     %---plot Score---
%     subplot(5,1,5);
%     hold on
%     for j=reshape(indx1,1,[])
%         line_style = 'r-';
%         indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
%         if isempty(indx)
%             continue;
%         end
%         plot([data_pos(indx(1)) data_pos(indx(end))],[score_seg(j) score_seg(j)], line_style, 'LineWidth',1.5);
%     end
%     axis ([min_pos max_pos -10 110])
% %     set(gca,'XTick',[]);
%     set(gca,'YTick',[0:50:100],'Box','on')
%     ylabel('Score');

    %save figure
%     figpath = [plotsdir '\Chr_' num2str(i) '_' barcode];
    figpath = [plotsdir '/Chr_' num2str(i) '_' barcode '.png'];
    %     eval(['print -djpeg -r600 ' figpath ])
    eval(['print -dpng -r400 ' figpath ])

end