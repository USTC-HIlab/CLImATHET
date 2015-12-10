function [p_states,num_loci,aCN,alpha,segments] = CLImATHET_process_results(beta,depend_table)
% 21/05/2015 by Zhenhua
%-----------------------------------------------------
%------over-all information of the cancer sample------
%p_states: proportions of all hidden states
%num_loci: total number of loci investigated
%aCN: averaged copy numbers of subclonal populations

global data_pos_ds_sep
global gamma_sep

tv_S = depend_table(:,2)~=0;
CN_mapping = depend_table(tv_S,3)'; % copy number of different entries

K = length(beta);

%initialize output parameters
num_loci = 0;
if nargout > 2
    alpha = zeros(K+1,1);
    segments = [];
end

%initialize intermediate variables
exp_num_states = [];
pos_dist = [];

for i = 1:length(gamma_sep) %for the ith chromosome
    post_probs = gamma_sep{i};
    data_pos = data_pos_ds_sep{i};
    
    %---handle p_states and num_loci---
    if isempty(exp_num_states) %initialization
        exp_num_states = zeros(size(post_probs,1),1);
    end
    exp_num_states = exp_num_states+sum(post_probs,2);
    num_loci = num_loci+size(post_probs,2);

    %---handle MAP states---
    if nargout > 2 %output predicted MAP states and subclonal populations 
        [temp,I] = max(post_probs,[],1);
        MAP_state = floor((I-1)/K)+1;
        MAP_sp = rem(I-1,K)+1;
        tv = MAP_state == max(MAP_state);
        MAP_state = MAP_state+1;
        MAP_state(tv) = 1;
        MAP_sp(tv) = 0;
        
        results = CLImATHET_segment_results(MAP_state,MAP_sp);
        segments = [segments; ones(size(results,1),1)*i results];    
        pos_dist = [pos_dist data_pos(results(:,2))-data_pos(results(:,1))+1]; 
        
        clear results;
    end

end
pos_dist = pos_dist';

%---handle p_states---
% p_states = exp_num_states/num_loci;
% p_states = [p_states(end) sum(reshape(p_states(1:end-1),K,(length(p_states)-1)/K),1)]';

p_states = zeros(length(CN_mapping),1);
for i = 1:length(CN_mapping)
    tv = segments(:,4) == i;
    p_states(i) = sum(pos_dist(tv))/sum(pos_dist);
end

%---handle aCN---
aCN = CN_mapping*p_states;

%---handle alpha---
if nargout > 3   
    temp = sum(segments(:,3)-segments(:,2)+1);
    for i = 1:K+1
        tv = segments(:,end) == i-1;
        alpha(i) = sum(segments(tv,3)-segments(tv,2)+1)/temp;
    end
end

end
