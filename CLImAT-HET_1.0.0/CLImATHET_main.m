function CLImATHET_paras = CLImATHET_main(init_CLImATHET_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose1,fn_nosuffix)
%10/12/2015 by Zhenhua
%CLImATHET main function,basically everything is done here
%--------------------------- screening -------------------------
global clamp_thres
global NoSolutionFlag
global mc_w
global pre_best_paras
clamp_thres = 1-1e-10;
mc_w = 0.8;
max_cp_num = 10;

pre_LL = -Inf;
pre_best_indx = [];
pre_CLImATHET_paras = [];
pre_best_paras = cell(6,1);
normal_pro = zeros(max_cp_num,1);

% bic_fid = fopen('INFO.bic', 'a+');
% fprintf(bic_fid,'----------------------------%s----------------------------\n',fn_nosuffix);
%-------------------------------------------------------------------
% Iteratively increase the number of subclonal populations to find 
% the optimal model using Bayesian information criterion (BIC)
%-------------------------------------------------------------------
for cp_num = 1:max_cp_num
    cp_flag = 0;
    %-------------------------------------------------------------------
    %               ---> d-sampling screening <---
    %-------------------------------------------------------------------
    %initialize parameters
    thres_del = 0.02;

    init_CLImATHET_paras = CLImATHET_Init_paras(init_CLImATHET_paras,depend_table,cp_num,0);
    [LL,CLImATHET_paras,p_states,num_loci,aCN] = CLImATHET_screening...
            (stepsize_ds,init_CLImATHET_paras,depend_table,thres_EM,max_iter,verbose1);
    %-------------------------------------------------------------------
    %               ---> summarize ds results <---
    %-------------------------------------------------------------------
    %get detailed info for further ananlysis
    %CN           0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7
    N_genotype = [0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4];
    tv = depend_table(:,2) ~= 0;
    tv_del = depend_table(tv,3) < 1;
    US_indx = depend_table(tv,1);
    tv_normal = depend_table(tv,3) == 2 & depend_table(tv,4) == 0.5;
    
    beta_all = cell2mat(CLImATHET_paras{3});
    p_all = cell2mat(CLImATHET_paras{4});
    lambda_all = cell2mat(CLImATHET_paras{5});

    p_total_del = sum(p_states(tv_del,:),1);
    p_total_normal = sum(p_states(tv_normal,:),1);
    DF_bias = (1-1./(N_genotype(US_indx)+1))*p_states;

    tv = (p_total_del<thres_del) & (aCN<4.5);
    if cp_num == 1
        NoSolutionFlag = false;
    end
    if ~any(tv)
        tv = ~tv;
        if cp_num > 1
            break;
        else
            warning('Can not find a feasible solution with pre-defined criteria!');
            NoSolutionFlag = true;
        end
    end
    sel_indx = find(tv);

    [temp,I] = max(LL(tv));
    best_indx = sel_indx(I);

    thres1_ratio = 0.005;%

    %correction for possible signal noise
    ratio = zeros(1,length(LL));
    for i = 1:length(LL)
        temp1 = LL(best_indx)-LL(i);
        temp2 = (log(num_loci(i))/2).*(DF_bias(best_indx)-DF_bias(i));
        if abs(temp1) <= 1 && abs(temp2) <= 0.01
            ratio(i) = 0;
        else
            ratio(i) = temp2/temp1;
        end
    end
    
    candi_best = best_indx;
    ds_candi1 = (ratio>=thres1_ratio) & tv;
    if any(ds_candi1)
        [temp,I] = max(ratio(ds_candi1));
%         [temp,I] = max(LL(ds_candi1));
        indx = find(ds_candi1);
        best_indx = indx(I);
    end
    normal_pro(cp_num) = p_total_normal(best_indx);
    if cp_num > 1
        if (normal_pro(cp_num-1)-normal_pro(cp_num))/(normal_pro(cp_num-1)+eps) > 0.3
            break;
        end
    end
    
%     % for test, save intermediate results
%     %-----------------------------------------------------------------------------------
%     fid = fopen('Model.details','a+');
%     fprintf(fid,'---------------------------%s----------------------------\n', fn_nosuffix);
%     fprintf(fid,'sp number:%d \n',cp_num);
%     for i = 1:length(ratio)
%         if i == best_indx
%             fprintf(fid,'2\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(i),ratio(i),LL(i));
%         elseif i == candi_best
%             fprintf(fid,'1\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(i),ratio(i),LL(i));
%         else
%             fprintf(fid,'0\t%f\t%f\t%f\t%f\t%f\n',p_total_del(i),beta_all(end,i),aCN(i),ratio(i),LL(i));
%         end
%     end
% %     fprintf(fid,'\n');
%     fprintf(fid,'----------------------------------------------------------\n');
%     fclose(fid);
    %-----------------------------------------------------------------------------------
    
    % Use Bayesian information criterion (BIC) to select the optimal model
    % With an increase in the number of subclonal populations (K to K+1), the number of parameters
    % that need to be estimated is increased by (2*K+1)*(S-1)^2+2*(S-1)+1
    
    alpha = 0.02;
    diff_th = 0;
    if cp_num == 1
        BIC_diff = -Inf;
        k = 0;
    else
        S = sum(depend_table(:,2)~=0);
        k = (2*(cp_num-1)+1)*(S-1)^2+2*(S-1)+1;
        BIC_diff = pre_LL-LL(best_indx)+alpha*0.5*k*log(num_loci(best_indx));
    end
%     fprintf(bic_fid,'%f\t%f\t%f\t%f\n',pre_LL,LL(best_indx),pre_LL-LL(best_indx),k*0.5*log(num_loci(best_indx)));
    
    if BIC_diff >= diff_th
        break;
    end
    
%     
    fprintf(1, '\n');
    disp('--------------- new optimal parameters -----------------');
    disp(['sp number:' num2str(cp_num)]);
    disp(['beta:' num2str(reshape(beta_all(:,best_indx),1,[]))]);
    disp(['p:' num2str(reshape(p_all(:,best_indx),1,[]))]);
    disp(['lambda:' num2str(lambda_all(best_indx)) ', LL:' num2str(LL(best_indx),'%5.1f')]);
    disp('--------------- new optimal parameters -----------------');
    fprintf(1, '\n');
    
    cp_flag = 1;
    pre_LL = LL(best_indx);
    pre_best_indx = best_indx;
    pre_CLImATHET_paras = CLImATHET_paras;
    init_CLImATHET_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]}];
    for i = 1:6
        pre_best_paras{i} = CLImATHET_paras{i}{best_indx};
    end
end

% fprintf(bic_fid,'----------------------------------------------------------\n\n');
% fclose(bic_fid);

% fid = fopen('normal_pro.txt', 'a+');
% fprintf(fid, '%s:', fn_nosuffix);
% for i = 1:max_cp_num
%     fprintf(fid, '%f\t', normal_pro(i));
% end
% fprintf(fid, '\n');
% fclose(fid);

% Now, use the selected parameters to call CNA/LOH
%-------------------------------------------------------------------
if cp_flag == 1
    init_CLImATHET_paras = CLImATHET_Init_paras(pre_CLImATHET_paras,depend_table,cp_num,pre_best_indx);
else
    init_CLImATHET_paras = CLImATHET_Init_paras(pre_CLImATHET_paras,depend_table,cp_num-1,pre_best_indx);
end
[temp,CLImATHET_paras] = CLImATHET_screening(1,init_CLImATHET_paras,depend_table,5*thres_EM,5,verbose1);
%-------------------------------------------------------------------

end


function CLImATHET_paras = CLImATHET_Init_paras(init_CLImATHET_paras,depend_table,cp_num,best_indx)
%This function is used to initialize/process parameters for CLImATHET training,
%best_indx is used to indicate which parameter configuration (ususally
%multiple generated in previous screening procedure) are selected. If
%best_indx equals 0, parameters will be initialized
global rc_median
global clamp_thres
global pre_best_paras

CLImATHET_paras = cell(1,8);
%parameter initialization
if best_indx == 0
    % copy neutral read counts
%     lambda_s = rc_median;
    lambda_s = [rc_median rc_median*2/3 rc_median/2];
    
    %---beta---, ###variable in grid searching###
    if isempty(init_CLImATHET_paras{3})
        %first generate initial parameters
        %   0<beta1<beta2<beta3<...<betaK<=1
        temp = linspace(0.15,0.99,5);
%         temp = [0.99];
        if cp_num == 1
            beta_0 = temp;% use 0 to indicate that no cells display aberrations in the loci
        else
%             beta_0 = [0.15 0.2 0.3;0.2 0.3 0.5;0.35 0.5 0.8];
%             beta_0 = [0.15;0.2;0.35];
            pre_beta = pre_best_paras{3};
            
% 			max_beta = pre_beta(end);
%             k = min(5,floor(max_beta*5)+2);
% 			temp = linspace(0.2,floor(max_beta*10)/10+0.09,k);
%             if length(temp) < cp_num
%                 beta_0 = linspace(0.2,0.99,cp_num)';
%             else
%                 beta_0 = combntns(temp,cp_num);
%                 beta_0 = beta_0';
%             end

            if cp_num == 2
                ext_value = 0.1:0.1:pre_beta-0.05;
                beta_0 = [ext_value;repmat(pre_beta,1,length(ext_value))];
%                 ext_value1 = 0.1:0.1:pre_beta-0.05;
%                 ext_value2 = pre_beta+0.05:0.1:0.99;
%                 beta_0 = [ext_value1 repmat(pre_beta,1,length(ext_value2));repmat(pre_beta,1,length(ext_value1)) ext_value2];
            else                                
                beta1 = 0.1:0.05:pre_beta(end)-0.05;
                temp = abs(repmat(beta1,length(pre_beta),1)-repmat(pre_beta,1,length(beta1)));
                x = min(temp,[],1);
                tv = x >= 0.05;
                beta1 = beta1(tv);
                beta_0 = zeros(cp_num,length(beta1));
                for i = 1:length(beta1)
                    tv1 = pre_beta > beta1(i);
                    tv2 = pre_beta < beta1(i);
                    beta_0(:,i) = [pre_beta(tv2);beta1(i);pre_beta(tv1)];
                end
            end           
        end
        if isempty(beta_0)
            beta_0 = linspace(0.2,0.99,cp_num)';
        end
        
        %---grid searching---
        if cp_num == 1
%             beta_0 = [0.99 0.8 0.6 0.4 0.2 0.15 0.1];
%             lambda_0 = [rc_median rc_median rc_median rc_median rc_median rc_median rc_median];           
            beta_0 = repmat(beta_0,1,length(lambda_s));
            lambda_0 = repmat(lambda_s,length(beta_0)/length(lambda_s),1);
            lambda_0 = lambda_0(:)';
            tv = ~((beta_0 <= 0.2) & (lambda_0 == rc_median*2/3 | lambda_0 == rc_median/2));
            beta_0 = beta_0(tv);
            lambda_0 = lambda_0(tv);
        end
        
    else %
        beta_0 = init_CLImATHET_paras{3};
    end
    
    
    N = size(beta_0,2);
    CLImATHET_paras{3} = mat2cell(beta_0,size(beta_0,1),ones(1,N));

    %---lambda---
    if isempty(init_CLImATHET_paras{5}) % lambda
        if cp_num == 1
            tv1 = lambda_0 == rc_median*2/3;
            tv2 = lambda_0 == rc_median/2;
            lambda_0(tv1) = rc_median*2./(2+beta_0(tv1));
            lambda_0(tv2) = rc_median./(1+beta_0(tv2));  
        else
            lambda_0 = repmat(pre_best_paras{5},1,N);
        end 
    else
        lambda_0 = init_CLImATHET_paras{5};
    end
    CLImATHET_paras{5} = mat2cell(lambda_0,size(lambda_0,1),ones(1,N));
    
    %---epsilon---
    if isempty(init_CLImATHET_paras{6}) % epsilon
          epsilon = 0;
    else
        epsilon = init_CLImATHET_paras{6};
    end
    CLImATHET_paras{6} = repmat({epsilon},1,N);

    %---pi---
    if isempty(init_CLImATHET_paras{1})
        S = sum(depend_table(:,2) ~= 0);
        K = size(beta_0,1);
        prior_0 = 1/((S-1)*K+1)*ones((S-1)*K+1,1);
    else
        prior_0 = init_CLImATHET_paras{1};
    end
    CLImATHET_paras{1} = repmat({prior_0},1,N);

    %---A---
    if isempty(init_CLImATHET_paras{2})
        S = sum(depend_table(:,2) ~= 0);
        K = size(beta_0,1);
        transmat_0 = norm_trans(ones((S-1)*K+1,(S-1)*K+1),clamp_thres);
    else
        transmat_0 = init_CLImATHET_paras{2};
    end
    CLImATHET_paras{2} = repmat({transmat_0},1,N); 
    
    %---p--- 
    if isempty(init_CLImATHET_paras{4})
%         if cp_num < 3
%             tv = depend_table(:,2) ~= 0;
%             cn_u = unique(depend_table(tv,3));
%             p_0 = ones(length(cn_u),1)*0.5;
%         else
%             p_0 = pre_best_paras{4};
%         end
        tv = depend_table(:,2) ~= 0;
        cn_u = unique(depend_table(tv,3));
        p_0 = ones(length(cn_u),1)*0.5;
    else
        p_0 = init_CLImATHET_paras{4};
    end
    CLImATHET_paras{4} = repmat({p_0},1,N);
    
    %---indicator vector---
    if isempty(init_CLImATHET_paras{7}) %indicator vector: '1' for update '0' fixed
        adj_all = ones(1,6);
    else
        adj_all = init_CLImATHET_paras{7};
    end
    CLImATHET_paras{7} = repmat({adj_all},1,N);

    %---other parameters---
    if isempty(init_CLImATHET_paras{8}) %parameters for observation function and EM algorithm
        other_paras = [1 1.9 1.9 1.9]; %normal prior(0,2,4,6 copies)
%         other_paras = [1 1.7 1.7 1.7];
    else
        other_paras = init_CLImATHET_paras{8};
    end
    CLImATHET_paras{8} = repmat({other_paras},1,N);

else %parse the results from previous screening
    %--pi--
    CLImATHET_paras{1} = init_CLImATHET_paras{1}(best_indx);
    %--A--
    CLImATHET_paras{2} = init_CLImATHET_paras{2}(best_indx);
    %--beta--
    CLImATHET_paras{3} = init_CLImATHET_paras{3}(best_indx);
    %--p--
    CLImATHET_paras{4} = init_CLImATHET_paras{4}(best_indx);
    %--lambda--
    CLImATHET_paras{5} = init_CLImATHET_paras{5}(best_indx);
    %--epsilon--
    CLImATHET_paras{6} = init_CLImATHET_paras{6}(best_indx);
    %--indicator vector--
    CLImATHET_paras{7} = init_CLImATHET_paras{7}(best_indx);
    %--parameters for observation function and EM algorithm--
    CLImATHET_paras{8} = {[1 1.9 1.9 1.9]};
%     CLImATHET_paras{8} = {[1 1.7 1.7 1.7]};
end

end
