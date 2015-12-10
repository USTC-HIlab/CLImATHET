function [obslik, condi_probs, condi_probs_fluct] = CLImATHET_get_obslik(chr,data_baf,data_rc,beta,p,lambda,depend_table,normal_prior)
% 04/02/2014 by Zhenhua

global gender
global Het_prior

max_rc = max(data_rc);
N = length(data_rc); %number of data points

ns = 2; %copy number of stromal cells
mus = 0.5;% BAF mean of stromal cells
tv_S = depend_table(:,2)~=0;
Nc = depend_table(tv_S,3); %vector of copy numbers of different entries
Muc = depend_table(tv_S,4); %vector of BAF means of different entries
cn_u = unique(Nc);% unique copy number

%----------------------------------------------------------------------%
%-----------------calculate all of the probablities--------------------%
Hom_prior = 1-Het_prior;
fluct_prob = 0.001;% 0.001

temp1 = repmat(Nc,1,length(beta));
temp2 = repmat(Muc,1,length(beta));
temp3 = repmat(beta',length(Nc),1);
if gender == 1 && (chr == 23 || chr == 24)
    Y = temp1.*temp3+(1-temp3);
    Z = temp1.*temp2.*temp3+(1-temp3);
else
    Y = temp1.*temp3+ns*(1-temp3);
    Z = temp1.*temp2.*temp3+ns*mus*(1-temp3);
end
p_het = Z./Y;
lambda_c = lambda*Y/2;

US_indx = depend_table(tv_S,1);
normal_indx = find(Nc == 2 & Muc == 0.5);

S = sum(tv_S);
K = length(beta);

obslik = zeros((S-1)*K+1,N);
condi_probs = zeros((S-1)*K+1,N);
condi_probs_fluct = zeros((S-1)*K+1,N);

obslik_bd_Homo = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),0.997);%0.997
obslik_fluct = 1./((data_baf(2,:)+1)*max_rc);

% normal state
k = cn_u == Nc(normal_indx);
obslik_rc = CLImATHET_eval_pdf_RC(data_rc,lambda,p(k));
obslik_bd_Het = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),0.5);

temp = find([1 3 8 15] == US_indx(normal_indx));
if ~isempty(temp)
    obslik_bd_Het = normal_prior(temp)*obslik_bd_Het; 
end
obslik_bd = Hom_prior*obslik_bd_Homo+Het_prior*obslik_bd_Het;

obslik((S-1)*K+1,:) = (1-fluct_prob)*obslik_rc.*obslik_bd+fluct_prob*obslik_fluct;
        
condi_probs((S-1)*K+1,:) = (1-fluct_prob)*obslik_rc.*(Het_prior*obslik_bd_Het)./obslik((S-1)*K+1,:);
condi_probs_fluct((S-1)*K+1,:) = fluct_prob*obslik_fluct./obslik((S-1)*K+1,:);

prior_w = [1 1 1 1 1 0.9 0.8 0.7];% 0,1,2,3,4,5,6,7 copy


for j = 1:K
    for i = 2:size(Y,1)
        if Nc(i) < 1
            fluct_prob = 0.001;%0.001
        else
            fluct_prob = 0.001;%0.001
        end
        if j < K
            k = floor(Nc(i))+1;
            wg = prior_w(k);
        else
            wg = 1;
        end
        if gender == 1 && (chr == 23 || chr == 24) && Muc(i) ~= 1
            obslik((i-2)*K+j,:) = (1-fluct_prob)*eps+fluct_prob*obslik_fluct;
            condi_probs((i-2)*K+j,:) = (1-fluct_prob)*eps./obslik((i-2)*K+j,:);
            condi_probs_fluct((i-2)*K+j,:) = (fluct_prob*obslik_fluct)./obslik(i(i-2)*K+j,:);
        else
            k = cn_u == Nc(i);
            obslik_rc = CLImATHET_eval_pdf_RC(data_rc,lambda_c(i,j),p(k));
            obslik_bd_Het = CLImATHET_eval_pdf_BD(data_baf(1,:),data_baf(2,:),p_het(i,j));
            temp = find([1 3 8 15] == US_indx(i));
            if ~isempty(temp)
                obslik_bd_Het = normal_prior(temp)*obslik_bd_Het; 
            end
            obslik_bd = Hom_prior*obslik_bd_Homo+Het_prior*obslik_bd_Het;

            obslik((i-2)*K+j,:) = (1-fluct_prob)*wg*obslik_rc.*obslik_bd+fluct_prob*obslik_fluct;

            condi_probs((i-2)*K+j,:) = (1-fluct_prob)*wg*obslik_rc.*(Het_prior*obslik_bd_Het)./obslik((i-2)*K+j,:);
            condi_probs_fluct((i-2)*K+j,:) = fluct_prob*obslik_fluct./obslik((i-2)*K+j,:);
        end
    end
end

end