function [LL, prior, transmat, beta, p, lambda, epsilon, nrIterations] = ...
    CLImATHET_EM_Newton(init_CLImATHET_paras,depend_table, thresh, max_iter,verbose)
% 10/12/2015 by Zhenhua
% This EM algorithm uses univariate method to update 
% parameters sperately, Newton-Raphson method is adopted

global clamp_thres

previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];

prior = init_CLImATHET_paras{1};
transmat = init_CLImATHET_paras{2};
beta = init_CLImATHET_paras{3};
p = init_CLImATHET_paras{4};
lambda = init_CLImATHET_paras{5};
epsilon = init_CLImATHET_paras{6};
normal_prior = init_CLImATHET_paras{8};

lambda_lim = [lambda-50 lambda+50];
if lambda_lim(1) < 50
    lambda_lim(1) = 50;
end

while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [loglik, exp_num_trans, exp_num_visits1, beta_u, p_u, lambda_u, epsilon_u] = ...
        CLImATHET_compute_ess(prior,transmat,beta,p,lambda,lambda_lim,epsilon,depend_table,normal_prior);
    
    converged = em_converged(loglik,previous_loglik,verbose,thresh);
    
    % update parameters
    if init_CLImATHET_paras{7}(1)
        prior = norm_trans(exp_num_visits1',0)';
    end
    if init_CLImATHET_paras{7}(2) && ~isempty(exp_num_trans)
        transmat = norm_trans(exp_num_trans,clamp_thres);
    end
    if init_CLImATHET_paras{7}(3) %update beta here
        beta = beta_u;
    end
    if init_CLImATHET_paras{7}(4) %update p here
        p = p_u;
    end
    if init_CLImATHET_paras{7}(5) %update lambda here
        lambda = lambda_u;
    end
    if init_CLImATHET_paras{7}(6) %update epsilon here
        epsilon = epsilon_u;
    end
    
    if verbose
        disp(['beta:' num2str(reshape(beta,1,[]))]);
        disp(['p:' num2str(reshape(p,1,[]))]);
        disp(['lambda:' num2str(lambda)]);
        fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik);
    end
    
    num_iter =  num_iter + 1;
    previous_loglik = loglik;
    LL = [LL loglik];
end
nrIterations = num_iter - 1;

end

%--------------------------------------------------------------------------
function [loglik, exp_num_trans, exp_num_visits1, beta_u, p_u, lambda_u, epsilon_u] = ...
    CLImATHET_compute_ess(prior,transmat,beta,p,lambda,lambda_lim,epsilon,depend_table,normal_prior)

global data_rc_ds_sep
global data_baf_ds_sep

global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rc_ds_sep); % each row is a sample
S_all = length(transmat); % number of all states 
exp_num_trans = zeros(S_all,S_all);
exp_num_visits1 = zeros(S_all,1);

%-----------------------E step-----------------------------
gamma_sep = cell(1,numex);
condi_probs_sep = cell(1,numex);
condi_probs_fluct_sep = cell(1,numex);
loglik = 0;

for ex = 1:numex

    % conditional probabilities
    [obslik, condi_probs, condi_probs_fluct] = CLImATHET_get_obslik(Chromosomes(ex),data_baf_ds_sep{ex},data_rc_ds_sep{ex},beta,p,lambda,depend_table,normal_prior);
    % Forward and Backward algorithm
    [temp1, gamma, current_ll, temp2, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    
    gamma_sep{ex} = gamma;
    clear gamma;
    condi_probs_sep{ex} = condi_probs;
    clear condi_probs;
    condi_probs_fluct_sep{ex} = condi_probs_fluct;
    clear condi_probs_nofluct;
end

%-----------------------M step-----------------------------
%update parameters

max_iter = 10;
% unpate beta
beta_tol = 0.01;
beta_u = CLImATHET_update_beta(beta,p,lambda,epsilon,depend_table,beta_tol,max_iter);
% beta_u = beta;

%update p
p_tol = 0.005;
p_u = CLImATHET_update_p(beta_u,p,lambda,epsilon,depend_table,p_tol,max_iter);
% p_u = p;

%update lambda
lambda_tol = 0.5;
lambda_u = CLImATHET_update_lambda(beta_u,p_u,lambda,epsilon,depend_table,lambda_tol,lambda_lim,max_iter);
lambda_u = round(lambda_u);
% lambda_u = lambda;

%update epsilon
% epsilon_tol = 0.5;
% epsilon_u = CLImATHET_update_epsilon(beta_u,p_u,lambda_u,epsilon,depend_table,epsilon_tol,max_iter);
% epsilon_u = round(epsilon_u);
epsilon_u = 0;

end

%--------------------------------------------------------------------------
function beta_u = CLImATHET_update_beta(beta,p,lambda,epsilon,depend_table,beta_tol,max_iter)

global data_baf_ds_sep
global data_rc_ds_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rc_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
mus = 0.5;%BAF of stromal cells
Num_US = 20; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries
Muc = depend_table(tv,4); %BAF values of different entries
cn_u = unique(Nc);% unique copy number
% normal_indx = find(Nc == 2 & Muc == 0.5);

p_all = zeros(length(Nc),1);
for i = 1:length(cn_u)
	tv = cn_u(i) == Nc;
	p_all(tv) = p(i);
end

iter = 0;
beta_u = zeros(size(beta));
% unConverged = ones(length(beta),1);
while 1
    %first order differential
    ELL_D_D_1 = zeros(length(beta),1);
    ELL_B_D_1 = zeros(length(beta),1);
    %second order differential
    ELL_D_D_2 = zeros(length(beta),1);
    ELL_B_D_2 = zeros(length(beta),1);
    
    temp1 = repmat(Nc,1,length(beta));
    temp2 = repmat(Muc,1,length(beta));
    temp3 = repmat(beta',length(Nc),1);
    Y = temp1.*temp3+ns*(1-temp3);
    Z = temp1.*temp2.*temp3+ns*mus*(1-temp3);
    clear temp1 temp2 temp3;
    v = Z./Y;
    lambda_c = lambda*Y/2+epsilon;
    v_1 = repmat(Nc*ns.*(Muc-mus),1,length(beta))./Y;
    v_2 = repmat(2*Nc*ns.*(Muc-mus).*(ns-Nc),1,length(beta))./Y.^3;
    
    part1 = lambda*(Nc-ns).*(1-p_all)./(2*p_all);
	part2 = lambda_c.*repmat((1-p_all)./p_all,1,length(beta));
    
    for ex = 1:numex
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        obs_baf = data_baf_ds_sep{ex};
        obs_rc = data_rc_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs_not_fluct = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        post_probs = gamma_sep{ex}(indx,:).*condi_probs_sep{ex}(indx,:); 
        for i = 2:size(Y,1) % 21/05/2015, now only calcualte heter entries in the dependent table                  
            %---beta---
            for j = 1:size(Y,2)
                ELL_B_D_1(j) = ELL_B_D_1(j)+post_probs((i-2)*K+j,:)*(obs_baf(1,:)*v_1(i,j)/v(i,j)-(obs_baf(2,:)-obs_baf(1,:))*v_1(i,j)/(1-v(i,j)))';
                ELL_B_D_2(j) = ELL_B_D_2(j)+post_probs((i-2)*K+j,:)*(obs_baf(1,:)*(v_2(i,j)*v(i,j)-v_1(i,j)^2)/v(i,j)^2-(obs_baf(2,:)-obs_baf(1,:))*(v_2(i,j)*(1-v(i,j))+v_1(i,j)^2)/(1-v(i,j))^2)';
                ELL_D_D_1(j) = ELL_D_D_1(j)+post_probs_not_fluct((i-2)*K+j,:)*(part1(i)*(psi(obs_rc+part2(i,j)+eps)+log(1-p_all(i))-psi(part2(i,j)+eps)))';
                ELL_D_D_2(j) = ELL_D_D_2(j)+post_probs_not_fluct((i-2)*K+j,:)*(part1(i)^2*(psi(1,obs_rc+part2(i,j)+eps)-psi(1,part2(i,j)+eps)))';
            end
            
        end
    end
    
    ELL_ALL_D_1 = ELL_D_D_1+ELL_B_D_1;    
    ELL_ALL_D_2 = ELL_D_D_2+ELL_B_D_2;
    
    beta_adj = -ELL_ALL_D_1./ELL_ALL_D_2;
    
    % now determine if the update violates the constrains
    % let beta_u = beta+c*beta_adj, c is the minimal coefficient from a list of
    % coefficients c_all that activiate but do not violate the constrains
    % 0.01<=beta(1)<beta(2)<...<beta(K)<=1
    c_p = 1;
    
    %-beta(1)<=-0.01 => beta(1)+c*beta_adj(1)>=0.01
    if -beta_adj(1) > 0 % not a feasible direction
        temp = (0.01-beta(1))/beta_adj(1);
        if temp < 1
%             disp(['Constrain: beta(' num2str(1) ')>' num2str(0) ') is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    for j = 2:K
        %-beta(j)+beta(j-1)<=-0.05 => beta(j)+c*beta_adj(j)>=beta(j-1)+c*beta_adj(j-1)+0.05
        if -beta_adj(j)+beta_adj(j-1) > 0 % not a feasible direction
            temp = (beta(j)-beta(j-1)-0.05)/(beta_adj(j-1)-beta_adj(j));
            if temp < 1
%                 disp(['Constrain: beta(' num2str(j-1) ')<beta(' num2str(j) ') is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
    end
    %beta(K)-1<=0 => beta(K)+c*beta_adj(K)<=1
    if beta_adj(K) > 0 % not a feasible direction
        temp = (0.99-beta(K))/beta_adj(K);
        if temp < 1
%             disp(['Constrain: beta(' num2str(K) ')<=1 is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    beta_u = beta+c_p*beta_adj;
    tv = isnan(beta_u);
    beta_u(tv) = beta(tv);
    
    iter = iter+1;
    
%     if sum(unConverged) == 0 || iter > max_iter
    if max(abs(beta_u-beta)) < beta_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        beta = beta_u;
    end    
    
end

end

%--------------------------------------------------------------------------
function lambda_u = CLImATHET_update_lambda(beta,p,lambda,epsilon,depend_table,lambda_tol,lambda_lim,max_iter)

global data_rc_ds_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rc_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries
Muc = depend_table(tv,4); %BAF values of different entries
cn_u = unique(Nc);% unique copy number
normal_indx = find(Nc == 2 & Muc == 0.5);

p_all = zeros(length(Nc),1);
for i = 1:length(cn_u)
	tv = cn_u(i) == Nc;
	p_all(tv) = p(i);
end

temp1 = repmat(Nc,1,length(beta));
temp3 = repmat(beta',length(Nc),1);
Y = temp1.*temp3+ns*(1-temp3);
clear temp1 temp3;

iter = 0;
while 1
    %first order differential
    ELL_D_D_1 = 0;
    %second order differential
    ELL_D_D_2 = 0;
    
    lambda_c = lambda*Y/2+epsilon;
    
	part2 = lambda_c.*repmat((1-p_all)./p_all,1,length(beta));
    part3 = Y.*repmat((1-p_all)*0.5./p_all,1,length(beta));
    
    for ex = 1:numex
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        obs_rc = data_rc_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs_not_fluct = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)
            % lambda    
            if i == normal_indx
                ELL_D_D_1 = ELL_D_D_1+post_probs_not_fluct(end,:)*(part3(i,1)*(psi(obs_rc+part2(i,1)+eps)+log(1-p_all(i))-psi(part2(i,1)+eps)))';
                ELL_D_D_2 = ELL_D_D_2+post_probs_not_fluct(end,:)*(part3(i,1)^2*(psi(1,obs_rc+part2(i,1)+eps)-psi(1,part2(i,1)+eps)))';
                continue;
            end  
            
            for j = 1:size(Y,2)
                ELL_D_D_1 = ELL_D_D_1+post_probs_not_fluct((i-2)*K+j,:)*(part3(i,j)*(psi(obs_rc+part2(i,j)+eps)+log(1-p_all(i))-psi(part2(i,j)+eps)))';
                ELL_D_D_2 = ELL_D_D_2+post_probs_not_fluct((i-2)*K+j,:)*(part3(i,j)^2*(psi(1,obs_rc+part2(i,j)+eps)-psi(1,part2(i,j)+eps)))';
            end
            
        end
    end
    lambda_adj = -ELL_D_D_1/ELL_D_D_2;
    c_p = 1;
    %-lambda<=-lambda_lim(1) => lambda+c*lambda_adj>=lambda_lim(1)
    if -lambda_adj > 0 % not a feasible direction
        temp = (lambda_lim(1)-lambda)/lambda_adj;
        if temp < 1
%             disp('Constrain: lambda>0 is active!!');
            if temp < c_p
                c_p = temp;
            end
        end
    end
    %lambda<=lambda_lim(2) => lambda+c*lambda_adj<=lambda_lim(2)
    if lambda_adj > 0 % not a feasible direction
        temp = (lambda_lim(2)-lambda)/lambda_adj;
        if temp < 1
%             disp('Constrain: lambda>0 is active!!');
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    lambda_u = lambda+c_p*lambda_adj;
    iter = iter+1;
	
    if abs(lambda_u-lambda) < lambda_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        lambda = lambda_u;
    end
    
end

end

%--------------------------------------------------------------------------
function epsilon_u = CLImATHET_update_epsilon(beta,p,lambda,epsilon,depend_table,epsilon_tol,max_iter)

global data_rc_ds_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rc_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries
Muc = depend_table(tv,4); %BAF values of different entries
cn_u = unique(Nc);% unique copy number
normal_indx = find(Nc == 2 & Muc == 0.5);

p_all = zeros(length(Nc),1);
for i = 1:length(cn_u)
	tv = cn_u(i) == Nc;
	p_all(tv) = p(i);
end

part1 = (1-p_all)./p_all;

temp1 = repmat(Nc,1,length(beta));
temp3 = repmat(beta',length(Nc),1);
Y = temp1.*temp3+ns*(1-temp3);
clear temp1 temp3;

iter = 0;
while 1
    %first order differential
    ELL_D_D_1 = 0;
    %second order differential
    ELL_D_D_2 = 0;
    
    lambda_c = lambda*Y/2+epsilon;
    
	part2 = lambda_c.*repmat(part1,1,length(beta));
    
    for ex = 1:numex
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        obs_rc = data_rc_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs_not_fluct = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)    
            % epsilon
            if i == normal_indx
                ELL_D_D_1 = ELL_D_D_1+post_probs_not_fluct(end,:)*(part1(i)*(psi(obs_rc+part2(i,1)+eps)+log(1-p_all(i))-psi(part2(i,1)+eps)))';
                ELL_D_D_2 = ELL_D_D_2+post_probs_not_fluct(end,:)*(part1(i)^2*(psi(1,obs_rc+part2(i,1)+eps)-psi(1,part2(i,1)+eps)))';
                continue;
            end
            
            for j = 1:size(Y,2)
                ELL_D_D_1 = ELL_D_D_1+post_probs_not_fluct((i-2)*K+j,:)*(part1(i)*(psi(obs_rc+part2(i,j)+eps)+log(1-p_all(i))-psi(part2(i,j)+eps)))';
                ELL_D_D_2 = ELL_D_D_2+post_probs_not_fluct((i-2)*K+j,:)*(part1(i)^2*(psi(1,obs_rc+part2(i,j)+eps)-psi(1,part2(i,j)+eps)))';
            end
            
        end
    end
    epsilon_adj = -ELL_D_D_1/ELL_D_D_2;
    c_p = 1;
    %-epsilon<=0 => epsilon+c*epsilon_adj>=0
    if -epsilon_adj > 0 % not a feasible direction
        temp = -epsilon/epsilon_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    %epsilon<lambda/10 => epsilon+c*epsilon_adj<lambda/10
    if epsilon_adj > 0 % not a feasible direction
        temp = (lambda/10-epsilon)/epsilon_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    epsilon_u = epsilon+c_p*epsilon_adj;
    iter = iter+1;
	
    if abs(epsilon_u-epsilon) < epsilon_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        epsilon = epsilon_u;
    end
    
end

end

%--------------------------------------------------------------------------
function p_u = CLImATHET_update_p(beta,p,lambda,epsilon,depend_table,p_tol,max_iter)


global data_rc_ds_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rc_ds_sep); % each row is a sample
K = length(beta);
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3); %copy numbers of different entries
Muc = depend_table(tv,4); %BAF values of different entries
cn_u = unique(Nc);% unique copy number
normal_indx = find(Nc == 2 & Muc == 0.5);

p_all = zeros(length(Nc),1);
for i = 1:length(cn_u)
	tv = cn_u(i) == Nc;
	p_all(tv) = p(i);
end

temp1 = repmat(Nc,1,length(beta));
temp3 = repmat(beta',length(Nc),1);
Y = temp1.*temp3+ns*(1-temp3);
lambda_c = lambda*Y/2+epsilon;
clear temp1 temp3;

p_u = zeros(size(p));
% unConverged = ones(length(p),1);

iter = 0;
while 1
    %first order differential
    ELL_D_D_1 = zeros(length(p),1);
    %second order differential
    ELL_D_D_2 = zeros(length(p),1);
    
    p_all = zeros(length(Nc),1);
    for i = 1:length(cn_u)
        tv = cn_u(i) == Nc;
        p_all(tv) = p(i);
    end
    
	part2 = lambda_c.*repmat((1-p_all)./p_all,1,length(beta));
    
    for ex = 1:numex
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        obs_rc = data_rc_ds_sep{ex};
        indx = 1:K*(size(Y,1)-1)+1;
        post_probs_not_fluct = gamma_sep{ex}(indx,:).*(1-condi_probs_fluct_sep{ex}(indx,:));
        for i = 1:size(Y,1)  
            % p
            k = find(Nc(i) == cn_u);
            if i == normal_indx
                ELL_D_D_1(k) = ELL_D_D_1(k)+post_probs_not_fluct(end,:)*(lambda_c(i,1)/p_all(i)^2*(psi(part2(i,1)+eps)-psi(obs_rc+part2(i,1)+eps)-log(1-p_all(i)))+(obs_rc-lambda_c(i,1))/p_all(i))';
                ELL_D_D_2(k) = ELL_D_D_2(k)+post_probs_not_fluct(end,:)*(1/p_all(i)^2*(2*lambda_c(i,1)/p_all(i)*(psi(obs_rc+part2(i,1)+eps)+log(1-p_all(i))-psi(part2(i,1)+eps))...
					+lambda_c(i,1)*(lambda_c(i,1)/p_all(i)^2*psi(1,obs_rc+part2(i,1)+eps)+1/(1-p_all(i))-lambda_c(i,1)/p_all(i)^2*psi(1,part2(i,1)+eps))-obs_rc+lambda_c(i,1)))';     
                continue;
            end  
                
            for j = 1:size(Y,2)
                ELL_D_D_1(k) = ELL_D_D_1(k)+post_probs_not_fluct((i-2)*K+j,:)*(lambda_c(i,j)/p_all(i)^2*(psi(part2(i,j)+eps)-psi(obs_rc+part2(i,j)+eps)-log(1-p_all(i)))+(obs_rc-lambda_c(i,j))/p_all(i))';
                ELL_D_D_2(k) = ELL_D_D_2(k)+post_probs_not_fluct((i-2)*K+j,:)*(1/p_all(i)^2*(2*lambda_c(i,j)/p_all(i)*(psi(obs_rc+part2(i,j)+eps)+log(1-p_all(i))-psi(part2(i,j)+eps))...
					+lambda_c(i,j)*(lambda_c(i,j)/p_all(i)^2*psi(1,obs_rc+part2(i,j)+eps)+1/(1-p_all(i))-lambda_c(i,j)/p_all(i)^2*psi(1,part2(i,j)+eps))-obs_rc+lambda_c(i,j)))';    
            end
   
        end
    end
    
    p_adj = -ELL_D_D_1./ELL_D_D_2;
    
    % now determine if the update violates the constrains
    % let p_u = p+c*p_adj, c is the minimal coefficient from a list of
    % coefficients c_all that activiate but do not violate the constrains
    % 0<p(i)<1
    M = length(p);

    for j = 1:M
        c_p = 1;
        %-p(j)<=-0.01 => p(j)+c*p_adj(j)>=0.01
        if -p_adj(j) > 0 % not a feasible direction
            temp = (0.01-p(j))/p_adj(j);
            if temp < 1
%                 disp(['Constrain: p(' num2str(j) ')>0 is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
        %p(j)<=0.99 => p(j)+c*p_adj(j)<=0.99
        if p_adj(j) > 0 % not a feasible direction
            temp = (0.99-p(j))/p_adj(j);
            if temp < 1
    %             disp(['Constrain: p(' num2str(j) ')<1 is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
        p_u(j) = p(j)+c_p*p_adj(j);
    end
    tv = isnan(p_u);
    p_u(tv) = p(tv);
    
    iter = iter+1;
    
%     if sum(unConverged) == 0 || iter > max_iter
    if max(abs(p_u-p)) < p_tol || iter > max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        p = p_u;
    end    
    
end

end

