function [data_bd, data_td] = CLImATHET_baf_tQN(data_bd, data_td)
% Quantile normalization of BAF
% 03/23/2014 by Zhenhua

flag1 = 0;
flag2 = 0;
if size(data_bd,2) > 1
    data_bd = data_bd';
    flag1 = 1;
end
if size(data_td,2) > 1
    data_td = data_td';
    flag2 = 1;
end

data_ad = data_td-data_bd;
a_fre = data_ad./data_td;
b_fre = data_bd./data_td;

% average distance before tQN
tv = b_fre < 0.9 & b_fre > 0.1;
dist1 = sum(abs(sort(a_fre(tv))-sort(b_fre(tv))))/sum(tv);

[temp1,indx1] = sort(a_fre);
[temp2,indx2] = sort(b_fre);

temp = [temp1 temp2];
mean_values = mean(temp,2);
temp1 = zeros(length(a_fre),1);
temp2 = zeros(length(b_fre),1);
temp1(indx1) = mean_values;
temp2(indx2) = mean_values;

% thre = 0.9;
thres = 0.5:0.05:1.5;
dists = zeros(1,length(thres));
for i = 1:length(thres)
    t1 = temp1;
    t2 = temp2;
    tv = t1./(a_fre+eps) > thres(i);
    t1(tv) = thres(i)*a_fre(tv);
    tv = t2./(b_fre+eps) > thres(i);
    t2(tv) = thres(i)*b_fre(tv);

    baf = t2./(t1+t2);
    aaf = 1-baf;

    tv = baf < 0.9 & baf > 0.1;
    x = sort(aaf(tv));
    y = sort(baf(tv));
%     figure(i);
%     clf;
%     plot(x, y, '.');
    dists(i) = sum(abs(x-y))/length(x);
end

[min_dist, I] = min(dists);
best_thre = thres(I);
if dist1 > min_dist
    tv = temp1./(a_fre+eps) > best_thre;
    temp1(tv) = best_thre*a_fre(tv);
    tv = temp2./(b_fre+eps) > best_thre;
    temp2(tv) = best_thre*b_fre(tv);
    a_fre = temp1;
    b_fre = temp2;
    baf = b_fre./(a_fre+b_fre);
    tv = baf > 0 & baf < 1;
    d = round((data_td(tv).*baf(tv)-data_bd(tv))./(1-baf(tv)));
    data_bd(tv) = data_bd(tv)+d;
    data_td(tv) = data_td(tv)+d;
end

if flag1
    data_bd = data_bd';
end
if flag2
    data_td = data_td';
end

end