function [f_CN, f_mCN, f_sp] = CLImATHET_segment_annotation(data_rc, data_baf, beta, p, lambda)
% 21/05/2015 by Zhenhua
% This function is used to estimate copy number, major copy number and
% subclonal population for highly amplified regions

global Het_prior
global max_copy

prior_w = [1 1.7 1.7 1.4 1.3 1 0.8 0.7];% 0,1,2,3,4,5,6,7 copy

Hom_prior = 1-Het_prior;
rc_median = median(data_rc);
obslik_bd_Homo = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),0.997);% homo
max_LL = -Inf;
f_CN = -1;
f_mCN = -1;
f_sp = -1;
for sp = 1:length(beta)
    CN = round(2*rc_median/(lambda*beta(sp)+eps)-2*(1-beta(sp))/(beta(sp)+eps));
    if CN < 0
        CN = 0;
    end
    if CN > max_copy
        j = max_copy+1;
    else
        j = CN+1;
    end
    if CN == 0
        Y = 0.001*beta(sp)+2*(1-beta(sp));
    else
        Y = CN*beta(sp)+2*(1-beta(sp));
    end  
    obslik_rc = CLImATHET_eval_pdf_RC(data_rc,lambda*Y/2,p(j));
    M = ceil(CN/2);
    LL = zeros(1,CN-M+1);
    for mCN = M:CN
        Z = mCN*beta(sp)+(1-beta(sp));
        obslik_bd_Het = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),Z/Y); % het
        if mCN/CN == 0.5 && CN > 1
            obslik_bd_Het = 1.7*obslik_bd_Het;
        end
        obslik_bd = Hom_prior*obslik_bd_Homo+Het_prior*obslik_bd_Het;
        if CN > 7
            k = 8;
        else
            k = CN+1;
        end
        obslik = prior_w(k)*obslik_rc.*obslik_bd;
        LL(mCN-M+1) = sum(log(obslik));
    end 
    [Y, I] = max(LL);
    temp = M:CN;
    if Y > max_LL
        max_LL = Y;
        f_CN = CN;
        f_mCN = temp(I);
        f_sp = sp;
    end
end

if f_CN == 2 && f_mCN == 1
    f_sp = 0;
end

end