function het_seq = CLImATHET_het_estimate(chr, data_baf, cn, mcn, beta)

global gender

if gender == 1 && (chr == 23 || chr == 24)
    Y = (1-beta)+beta*cn;
    Z = (1-beta)+beta*mcn;
else
    Y = (1-beta)*2+beta*cn;
    Z = (1-beta)+beta*mcn;
end
obslik_bd_Homo = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),0.997);% homo
if cn > 0 && mcn/cn == 0.5
    obslik_bd_Het = 1.9*CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),Z/Y); % het
else
    obslik_bd_Het = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),Z/Y); % het
end
het_seq = obslik_bd_Het > obslik_bd_Homo;

end