function CLImATHET_plot_results(Datafile,resultsfile,plotsdir,barcode)
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
%--------------- plot figures ---------------------%
rc_colors = [0.5 0.5 0.5;
             0 0.9 0;
             0 0 0.9;
             0.9 0 0];
baf_colors = [0 250 0;
             0 0 250;
             250 0 0]./255;

Chromosomes = intersect(unique(data_chr_all),1:30);
max_pos = zeros(1,length(Chromosomes));

for i = 1:length(Chromosomes)
    tv1 = data_chr_all == Chromosomes(i);
    max_pos(i) = max(data_pos_all(tv1));
end

ratio = max_pos/sum(max_pos);
xtick = cumsum([0 ratio(1:end-1)])+ratio/2;

clf;
line_style = 'k-';
LineWidth = 1.0;
MarkerSize = 2;
%CN
subplot(4,1,1);
hold on
set(gca,'YGrid','on');
set(gca,'FontSize',FontSize);
set(gca,'YTick',[0:1:7],'Box','on');
set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
ylabel('Copy number');
pre_x = 0;
chr_epos = zeros(length(Chromosomes),1);
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_pos_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
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
        plot([x(indx(1)) x(indx(end))],[CN CN], 'r-', 'LineWidth',2.5);
        if CN == mCN
            plot([x(indx(1)) x(indx(end))],[mCN-0.25 mCN-0.25], 'b-', 'LineWidth',2.5);
        else
            plot([x(indx(1)) x(indx(end))],[mCN mCN], 'b-', 'LineWidth',2.5);
        end
    end
    chr_epos(i) = max(x);
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-0.1 7.1], line_style, 'LineWidth',LineWidth)
end
set(gca,'XTick',[])
axis([0 1 -0.1 7.1]);
barcode_m = barcode;
tmp = strfind(barcode_m,'_');
barcode_m(tmp) = '-';
title(barcode_m);

%BAF
subplot(4,1,2);
set(gca,'FontSize',FontSize);
hold on
pre_x = 0;
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_pos_all(tv);
    data_bd = data_bd_all(tv);
    data_td = data_td_all(tv);
    data_baf = data_bd./data_td;
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
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
        plot(x(tv),data_baf(tv),'.','MarkerSize',MarkerSize, 'Color', baf_colors(k,:));
    end
    for j = 0:0.25:1
        plot([x(1) x(end)],[j j],'k-','LineWidth',0.5)
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
        plot([x(indx(1)) x(indx(end))],[baf_mean baf_mean],'k-','LineWidth',1.2);
        plot([x(indx(1)) x(indx(end))],[1-baf_mean 1-baf_mean],'k-','LineWidth',1.2);
    end
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-0.03 1.03], line_style, 'LineWidth',LineWidth)
end
ylabel('BAF');
axis([0 1 -0.03 1.03])
set(gca,'YTick',[0 0.5 1],'Box','on')
set(gca,'XTick',[])

max_rc = max(data_rc_all);
subplot(4,1,3);
set(gca,'FontSize',FontSize);
hold on
pre_x = 0;
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_pos_all(tv);
    data_rc = data_rc_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
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
        plot(x(tv),data_rc(tv),'.','MarkerSize',MarkerSize, 'Color', rc_colors(k,:));
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
        plot([x(indx(1)) x(indx(end))],[rc_median rc_median],'k-','LineWidth',1.2);
    end
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [0 max_rc+10], line_style, 'LineWidth',LineWidth)
end
ylabel('Read counts');
axis ([0 1 0 max_rc+10])
set(gca,'Box','on')
set(gca,'XTick',[])
% axis([0 1 -Inf Inf])

subplot(4,1,4);
set(gca,'FontSize',FontSize);
hold on
pre_x = 0;
for i = 1:length(Chromosomes)
    tv = data_chr_all == Chromosomes(i);
    data_pos = data_pos_all(tv);
    x = data_pos*ratio(i)/max_pos(i)+pre_x;
    indx1 = find(chr_seg == Chromosomes(i));
    for j = reshape(indx1,1,[])
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([x(indx(1)) x(indx(end))],[purity_seg(j) purity_seg(j)], line_style, 'LineWidth',2.5);
    end
    pre_x = pre_x+ratio(i);  
end
for i = 1:length(Chromosomes)-1
    plot([chr_epos(i) chr_epos(i)], [-0.03 1.03], line_style, 'LineWidth',LineWidth)
end
axis ([0 1 -0.03 1.03])
set(gca,'YTick',[0:0.5:1],'Box','on')
ylabel('Cellularity');
xlabel('Chromosome');
set(gca,'Box','on')
set(gca,'XTick',xtick);
set(gca,'XTickLabel',mat2cell(reshape(Chromosomes,1,[]),1,length(Chromosomes)));

%save figure
figpath = [plotsdir '/' barcode '.png'];
eval(['print -dpng -r400 ' figpath ])

end