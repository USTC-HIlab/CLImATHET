function [LL_all,CLImATHET_paras,p_states,num_loci,aCN] = CLImATHET_screening(stepsize1,init_CLImATHET_paras,depend_table,thres1,max_iter1,verbose1)
% 10/12/2015 by Zhenhua

global data_rc_sep
global data_baf_sep
global data_pos_sep

global data_rc_ds_sep
global data_baf_ds_sep
global data_pos_ds_sep

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_CLImATHET_paras{1};
transmat_all = init_CLImATHET_paras{2};
beta_all = init_CLImATHET_paras{3};
p_all = init_CLImATHET_paras{4};
lambda_all = init_CLImATHET_paras{5};
epsilon_all = init_CLImATHET_paras{6};
indivec_all = init_CLImATHET_paras{7};
otherparas_all = init_CLImATHET_paras{8};

numex = length(data_rc_sep); % each row is a sample
data_rc_ds_sep = cell(1,numex);
data_baf_ds_sep = cell(1,numex);
data_pos_ds_sep = cell(1,numex);

for ex = 1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_rc_sep{ex});
        data_baf_ds_sep{ex} = data_baf_sep{ex}(:,indx_ds);
        data_rc_ds_sep{ex} = data_rc_sep{ex}(indx_ds);
        data_pos_ds_sep{ex} = data_pos_sep{ex}(indx_ds);
    else %no ds
        data_baf_ds_sep{ex} = data_baf_sep{ex};
        data_rc_ds_sep{ex} = data_rc_sep{ex};
        data_pos_ds_sep{ex} = data_pos_sep{ex};
    end    
end

LL_all = zeros(1,length(beta_all));
CLImATHET_paras = cell(1,8);
if nargout > 2
    p_states = [];
    num_loci = zeros(1,length(beta_all));
    aCN = zeros(1,length(beta_all));
end

for i = 1:length(beta_all)
    %1x1 cell
    init_CLImATHET_paras(1) = prior_all(i);
    init_CLImATHET_paras(2) = transmat_all(i);
    init_CLImATHET_paras(3) = beta_all(i);
    init_CLImATHET_paras(4) = p_all(i);
    init_CLImATHET_paras(5) = lambda_all(i);
    init_CLImATHET_paras(6) = epsilon_all(i);
    init_CLImATHET_paras(7) = indivec_all(i);
    init_CLImATHET_paras(8) = otherparas_all(i);
    
    [LL, prior, transmat, beta, p, lambda, epsilon, iterations] = CLImATHET_EM_Newton(init_CLImATHET_paras,depend_table,thres1,max_iter1,verbose1);
        
    LL_all(i) = LL(end);
    CLImATHET_paras{1} = [CLImATHET_paras{1} {prior}];
    CLImATHET_paras{2} = [CLImATHET_paras{2} {transmat}];
    CLImATHET_paras{3} = [CLImATHET_paras{3} {beta}];
    CLImATHET_paras{4} = [CLImATHET_paras{4} {p}];
    CLImATHET_paras{5} = [CLImATHET_paras{5} {lambda}];
    CLImATHET_paras{6} = [CLImATHET_paras{6} {epsilon}];
    CLImATHET_paras{7} = [CLImATHET_paras{7} init_CLImATHET_paras(7)];
    CLImATHET_paras{8} = [CLImATHET_paras{8} init_CLImATHET_paras(8)];
    
    if nargout > 2
%         [temp,num_loci(i),aCN(i)] = CLImATHET_process_results(beta,depend_table);
        [temp,num_loci(i),aCN(i)] = CLImATHET_process_results(beta,depend_table);
		p_states = [p_states temp];
    end

    if verbose1
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, iterations:' num2str(iterations) ', sp_num:' num2str(length(beta))]);
        disp(['beta:' num2str(reshape(beta,1,[]))]);
        disp(['p:' num2str(reshape(p,1,[]))]);
        disp(['lambda:' num2str(lambda)]);
        disp(['LL:' num2str(LL(end),'%5.1f')]);
        disp('--------------- screening report -----------------')
    end
    
end

end