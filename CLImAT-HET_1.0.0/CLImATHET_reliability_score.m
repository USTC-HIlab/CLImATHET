function score = CLImATHET_reliability_score(chr,data_rc,data_baf,het_tv,beta,p,lambda,cn,mcn,sp)
% 21/05/2015 by Zhenhua
% This function is used to evaluate the reliability of CLImATHET results

global gender
global max_copy

if sp == 0
    sp = 1;
end
if gender == 1 && (chr == 23 || chr == 24)
    Y = cn*beta(sp)+(1-beta(sp));
    Z = mcn*beta(sp)+(1-beta(sp));
else
    Y = cn*beta(sp)+2*(1-beta(sp));
    Z = mcn*beta(sp)+(1-beta(sp));
end
p_het = Z/Y;
lambda_c = lambda*Y/2;
if cn > max_copy
    j = max_copy+1;
else
    j = cn+1;
end
e_rc = lambda_c-p(j)/(1-p(j));
if e_rc < 0
    e_rc = 0;
end

if sum(het_tv) > 0
    score_bd = CLImATHET_eval_pdf_BD(data_baf(1,het_tv),data_baf(2,het_tv),p_het)./CLImATHET_eval_pdf_BD(round(data_baf(2,het_tv)*p_het),data_baf(2,het_tv),p_het);
    score_bd(score_bd > 1) = 1;
    score_rc = CLImATHET_eval_pdf_RC(data_rc(het_tv),lambda_c,p(j))/CLImATHET_eval_pdf_RC(e_rc,lambda_c,p(j));
    score_rc(score_rc > 1) = 1; 
    score = mean(score_bd.*score_rc);
else
    score_bd = CLImATHET_eval_pdf_BD(data_baf(1,~het_tv),data_baf(2,~het_tv),0.997)./CLImATHET_eval_pdf_BD(round(data_baf(2,~het_tv)*0.997),data_baf(2,~het_tv),0.997);
    score_bd(score_bd > 1) = 1;
    score_rc = CLImATHET_eval_pdf_RC(data_rc(~het_tv),lambda_c,p(j))/CLImATHET_eval_pdf_RC(e_rc,lambda_c,p(j));
    score_rc(score_rc > 1) = 1; 
    score = mean(score_bd.*score_rc);
end

if cn > 0 && mcn/cn == 0.5
    score = score*1.2;
end
