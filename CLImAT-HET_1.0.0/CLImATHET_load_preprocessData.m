function [data_chr_all, data_pos_all, data_bd_all, data_td_all, data_rc_all] = CLImATHET_load_preprocessData(df_file)
% 10/12/2015 by Zhenhua

global window
global gender

global td_th
global gc_th
global map_th

% parameters to filter data
rc_th = [0.1 99.5];
hom_th = 0.01;

% load read counts data
fid = fopen(df_file, 'r');
if fid == -1
    error(['Can not open file ' df_file]);
end
fline = fgetl(fid);
result = regexp(fline,'#window\s*size\s*=\s*(\S+)','tokens','once');
window = str2double(result{1});
results = textscan(fid, '%f%f%f%f%f%f%f', 'HeaderLines', 1, 'treatAsEmpty', {'NA', 'na'});% only autosomes are considered here
data_chr_all = results{1};
data_pos_all = results{2};
data_bd_all = results{3};% B alellic count
data_td_all = results{4};% alellic total count
data_rc_all = results{5};% read count
data_gc_all = results{6};% GC-content
data_map_all = results{7};% GC-content
clear results;
fclose(fid);

tv = data_chr_all == 24;
if sum(tv) > 0
    gender = 1;
else
    gender = 0;
end
clear tv;

fprintf(1,'Total %d windows are loaded from file "%s".\n',length(data_chr_all),df_file);
if length(data_chr_all) < 100000
    fprintf(1,'Warning: the number of windows loaded is too small, check the format of data file and whether data are completely loaded!\n');
end

% initial filtering
Chromosomes = intersect(unique(data_chr_all),1:24); % only use autosome
low_p = prctile(data_rc_all,rc_th(1));
high_p = prctile(data_rc_all,rc_th(2));
tv = ismember(data_chr_all,Chromosomes) & (data_gc_all > gc_th(1) & data_gc_all < gc_th(2)) & (data_map_all > map_th(1) & data_map_all < map_th(2))...
    & (data_rc_all > low_p & data_rc_all < high_p) & (data_td_all >= td_th(1) & data_td_all <= td_th(2));
data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bd_all = data_bd_all(tv);
data_td_all = data_td_all(tv);
data_rc_all = data_rc_all(tv);
data_gc_all = data_gc_all(tv);
data_map_all = data_map_all(tv);
clear tv;

% sort by position
for i = 1:length(Chromosomes)
    tv = ismember(data_chr_all, Chromosomes(i));
    temp = sortrows([data_pos_all(tv) data_bd_all(tv) data_td_all(tv) data_rc_all(tv) data_gc_all(tv) data_map_all(tv)], 1);
    data_pos_all(tv) = temp(:,1);
    data_bd_all(tv) = temp(:,2);
    data_td_all(tv) = temp(:,3);
    data_rc_all(tv) = temp(:,4);
    data_gc_all(tv) = temp(:,5);
    data_map_all(tv) = temp(:,6);
end
clear tv temp;

data_rc_all = CLImATHET_GC_MAP_correction(data_rc_all,data_gc_all,data_map_all);

clear data_map_all data_gc_all;

high_p = prctile(data_rc_all,99.8);
tv = data_rc_all < high_p;
data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bd_all = data_bd_all(tv);
data_td_all = data_td_all(tv);
data_rc_all = data_rc_all(tv);
clear tv;

% Quantile normalization of BAF
[data_bd_all, data_td_all] = CLImATHET_baf_tQN(data_bd_all, data_td_all);

% make sure that the number of heterozygous positions is comparative with the number of homozygous positions
tv = binopdf(data_td_all-data_bd_all,data_td_all,0.01) < hom_th & binopdf(data_bd_all,data_td_all,0.01) < hom_th;
if sum(tv) < sum(~tv)
    indx = find(~tv); 
    step = max(1,floor(length(indx)/(0.5*sum(tv))));
%     step = max(1,floor(length(indx)/sum(tv)));
    ds = 1:step:length(indx);
    tv(indx(ds)) = 1;
    clear indx ds
end
data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bd_all = data_bd_all(tv);
data_td_all = data_td_all(tv);
data_rc_all = data_rc_all(tv);
clear tv;

end