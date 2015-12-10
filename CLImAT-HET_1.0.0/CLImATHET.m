function CLImATHET(Datasource, Outputsource, configfile)
% function CLImATHET(Datasource, Outputsource, configfile)
%--------------------------------------------------------------------%
%------------------>       version 1.0       <---------------------
%--------------------------------------------------------------------%
% 10/12/2015 by Zhenhua
% This is the first version of the CLImATHET method

%------Input and output------%
% Datasource: a directory containing data files to be processed
% Outputsource: directory of output files
% configfile: a file containing user-defined parameters

global current_version
global NoSolutionFlag
current_version = '1.0';
sourceflag = 1;

if nargin < 3
    error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m'] );
end

s = mkdir(Outputsource);
if ~s %failed to make a directory
    error(['Can not make directory: ' Outputsource]);
end

%parameters used in CLImATHET
%===============================================
%---for cn up to 7---
depend_table = [...
    %w1
    3 1 2 0.5;...
    1 1 0.001 0.5;...
    2 1 1 1.0;...
    4 1 2 1.0;...
    5 1 3 0.67;...
    6 1 3 1.0;...
    7 1 4 0.75;...
    8 1 4 0.5;...
    9 1 4 1.0;...
    10 1 5 0.8;...
    11 1 5 0.6;...
    12 1 5 1.0;...
    13 1 6 5/6;...
    14 1 6 4/6;...
    15 1 6 0.5;...
    16 1 6 1.0;...
    17 1 7 6/7;...
    18 1 7 5/7;...
    19 1 7 4/7;...
    20 1 7 1.0
    ];
%ABBBBB,AABBBB,AAABBB,BBBBBB
%ABBBBBB,AABBBBB,AAABBBB,BBBBBBB
%===============================================
thres_EM = 1e-4;
max_iter = 30;
verbose = 0;
init_CLImATHET_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]}]; % initial parameters will be assigned in the main function
                                                                   % parameters:pie,transmat,alpha,beta,p,lambda,epsilon

%initialization of global variable
global data_rc_sep
global data_baf_sep
global data_pos_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global tumor_range
global Chromosomes
global rc_median
global gender
global Het_prior
global max_copy

Het_prior = 0.7;
max_copy = max(depend_table(:,3));

% read configuration file
CLImATHET_read_config(configfile);

if sourceflag == 1 % reading a directory
    disp(['CLImAT-HET (version ' current_version ') is loading...'])
    datafilelist = dir(Datasource);
    if length(datafilelist) < 3
        error(['No data files in the directory ' Datasource]);
    else %now do batch annotation     
        filename = cell(1,(length(datafilelist)-2));
        for i = 3:length(datafilelist)
            filename{i-2} = datafilelist(i).name;
        end
        tumor_range_table = [0.1*ones(length(filename),1) ones(length(filename),1)];
        
        %--------------perform CLImAT-HET--------------------      
        %record all the time used for batch annotation
        t_all = 0;
        if length(filename) == 1
            disp('-----CLImAT-HET batch annotation starts now, ONE file is found-----');
        elseif length(filename) > 1
            disp(['-----CLImAT-HET batch annotation starts now, ' num2str(length(filename)) ' files are found-----']);
        end
        for fileindx = 1:length(filename)
            %clear global variables            
            gamma_sep = [];
            condi_probs_sep = [];
            condi_probs_fluct_sep = [];
            
            %record time cost
            tic
            tumor_range = tumor_range_table(fileindx,:);
            results = regexp(filename{fileindx}, '^(.+)\.+.+','tokens', 'once');
            if isempty(results)
                fn_nosuffix = filename{fileindx};
            else
                fn_nosuffix = results{1};
                if ~isempty(strfind(fn_nosuffix,'.'))
                    fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
                end
            end
            
            %--------open result files --------------
            o_fid = fopen([Outputsource '/' fn_nosuffix '.results'],'w');
            if o_fid == -1
                error(['Can not open result file for ' filename{fileindx}]);
            end
            
            %--------------load data--------------------
            [data_chr_all, data_pos_all, data_bd_all, data_td_all, data_rc_all] = CLImATHET_load_preprocessData([Datasource '/' filename{fileindx}]);

            rc_median = median(data_rc_all);
            Chromosomes = reshape(unique(data_chr_all),1,[]);
            
            nfilename = [Outputsource '/' fn_nosuffix '_normalized.mat'];
            
            eval(['save ' nfilename ' data_chr_all data_pos_all data_bd_all data_td_all data_rc_all']);
            
            if size(data_pos_all,1)>size(data_pos_all,2) %make sure it's 1byN
                data_pos_all = data_pos_all';
            end
            if size(data_bd_all,1)>size(data_bd_all,2) %make sure it's 1byN
                data_bd_all = data_bd_all';
            end
            if size(data_td_all,1)>size(data_td_all,2) %make sure it's 1byN
                data_td_all = data_td_all';
            end
            if size(data_rc_all,1)>size(data_rc_all,2) %make sure it's 1byN
                data_rc_all = data_rc_all';
            end
            
            %use at least 30000 data points for screening
            stepsize_ds = max(floor(length(data_chr_all)/30000),1);
            
            %-------divide into different chromosomes----------
            chr_num = length(Chromosomes);
            data_pos_sep = cell(chr_num,1);
            data_rc_sep = cell(chr_num,1);
            data_baf_sep = cell(chr_num,1);
            data_bafA_sep = cell(chr_num,1);
            data_het_sep = cell(chr_num,1);
                        
            for i = 1:chr_num
                tv = ismember(data_chr_all,Chromosomes(i));
                data_pos_sep{i} = data_pos_all(tv);
                data_rc_sep{i} = data_rc_all(tv);
                
                data_bd = data_bd_all(tv);
                data_td = data_td_all(tv);
                data_baf = data_bd./data_td;               
                tvA = data_baf < 0.5;           
                data_bafA_sep{i} = tvA; 
                data_bd(tvA) = data_td(tvA)-data_bd(tvA);      
                data_baf_sep{i} = [data_bd; data_td];
                data_het_sep{i} = zeros(1,length(data_rc_all(tv))) > 0;
                clear tv tvA data_baf data_bd data_td;
            end
            clear data_chr_all data_pos_all data_bd_all data_td_all data_rc_all;
            
            %------------------ call CLImATHET --------------------
            CLImATHET_paras = CLImATHET_main(init_CLImATHET_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose,fn_nosuffix);
            
            beta = CLImATHET_paras{3}{1};           
            p = CLImATHET_paras{4}{1};
            lambda = CLImATHET_paras{5}{1};
            
            %-------------- process results ----------------------          
%             [p_states,num_loci,aCN,alpha,segments] = CLImATHET_process_results(beta,depend_table);
            [p_states,num_loci,aCN,alpha,segments] = CLImATHET_process_results(beta,depend_table);
            
            tv_S = depend_table(:,2)~=0;
            cn = depend_table(tv_S,3)';
            cn(cn<1) = 0;
            Bcn = round(cn.*depend_table(tv_S,4)');
            AI_mapping = [3 1 1 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];%DEL(1) LOH(2) HET(3)
%             cn = [0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];
%             Bcn = [0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
%             AI_mapping = [1 1 3 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];%DEL(1) LOH(2) HET(3)
                       
            cn_segs_all = zeros(size(segments,1),1); % copy number
            mcn_segs_all = zeros(size(segments,1),1); % major copy number 
            AI_segs_all = zeros(size(segments,1),1); % allelic imblance
            sp_segs_all = zeros(size(segments,1),1); % subclonal population
            scores = zeros(size(segments,1),1); % reliability socre
            
            bp_len = 0;
            cn_w = 0;
            
            for i = 1:size(segments,1)
                chr_indx = segments(i,1);
                s_indx = segments(i,2);
                e_indx = segments(i,3);
                state_indx = segments(i,4);
                sp_indx = segments(i,5);
                St_pos = data_pos_sep{chr_indx}(s_indx);
                Ed_pos = data_pos_sep{chr_indx}(e_indx);
                
                data_rc = data_rc_sep{chr_indx}(s_indx:e_indx);
                data_baf = data_baf_sep{chr_indx}(:,s_indx:e_indx);

                % estimate copy number, minor copy number and subclonal population for highly amplified regions
%                 if cn(state_indx) == max_copy-1
                rc_median = median(data_rc);
                if sp_indx == 0
                    sp = 1;
                else
                    sp = sp_indx;
                end
                if gender == 1 && (Chromosomes(Chr_i) == 23 || Chromosomes(Chr_i) == 24)
                    CN = round((rc_median*2/lambda-(1-beta(sp)))/(beta(sp)+eps));
                else
                    CN = round(2*(rc_median/lambda-(1-beta(sp)))/(beta(sp)+eps));
                end
                if CN == cn(state_indx)
                    if gender == 1 && (Chromosomes(chr_indx) == 23 || Chromosomes(chr_indx) == 24)
                        cn_segs_all(i) = cn(state_indx);
                        mcn_segs_all(i) = cn(state_indx);
                        AI_segs_all(i) = -1;
                    else
                        cn_segs_all(i) = cn(state_indx);
                        mcn_segs_all(i) = Bcn(state_indx);
                        AI_segs_all(i) = AI_mapping(state_indx);
                    end
                    sp_segs_all(i) = sp_indx;
                else
                    if gender == 1 && (Chromosomes(Chr_i) == 23 || Chromosomes(Chr_i) == 24)
                        cn_segs_all(i) = cn(state_indx);
                        mcn_segs_all(i) = cn(state_indx);
                        AI_segs_all(i) = -1;
                        sp_segs_all(i) = sp_indx;
                    else
                        [CN, mCN, SP] = CLImATHET_segment_annotation(data_rc, data_baf, beta, p, lambda);
                        cn_segs_all(i) = CN;
                        mcn_segs_all(i) = mCN;
                        sp_segs_all(i) = SP;
                        if CN < 2
                            AI_segs_all(i) = 1;
                        elseif mCN == CN
                            AI_segs_all(i) = 2;
                        else
                            AI_segs_all(i) = 3;
                        end
                    end                  
                end
                
                if sp_segs_all(i) == 0
                    sp = 1;
                else
                    sp = sp_segs_all(i);
                end
                data_het_sep{chr_indx}(s_indx:e_indx) = CLImATHET_het_estimate(Chromosomes(chr_indx), data_baf, cn_segs_all(i), mcn_segs_all(i), beta(sp));
                scores(i) = CLImATHET_reliability_score(Chromosomes(chr_indx),data_rc,data_baf,data_het_sep{chr_indx}(s_indx:e_indx),beta,p,lambda,cn_segs_all(i),mcn_segs_all(i),sp_segs_all(i));
                
                bp_len = bp_len+(Ed_pos-St_pos+1); 
                cn_w = cn_w+(Ed_pos-St_pos+1)*cn_segs_all(i);
            end
            
            ACN = cn_w/bp_len;
            
            m_score = (segments(:,3)-segments(:,2)+1)'*scores;
            m_score = m_score/sum(segments(:,3)-segments(:,2)+1);
            
            scores = Scale_score(scores,m_score);
            
            listfile = 'LOG.txt';
            fp_list = fopen(listfile,'a+');
            if fp_list == -1
                warning('Can not open the file: "%s"!\n',listfile);
            else
                if fileindx == 1
                    fprintf(fp_list,'Version:\tDate\tTime\tSample\tTumor purity\tSubclonal populations\tAverage copy number\tCopy neutral RC\n');
%                     fprintf(fp_list,'Version:\tDate\tTime\tSample\tTumor purity\tTumor ploidy\tCopy neutral RD\tEpsilon\n');
                end
                fprintf(fp_list,'CLImATHET%s\t%s\t%s\t%s\t%f\t',current_version,datestr(clock,'mmm-dd-yyyy'),datestr(clock,'HH:MM:SS'),fn_nosuffix,beta(end));
                for i = 1:length(beta)
                    if i == length(beta)
                        fprintf(fp_list,'%f\t',beta(i));
                    else
                        fprintf(fp_list,'%f,',beta(i));
                    end
                end
                fprintf(fp_list,'%f\t%3.0f\n',ACN,lambda);
            end
            fclose(fp_list);
            
            %-------------- output summary of the results--------------
            fprintf(o_fid,'---------------------------------------------------------------\n');
            fprintf(o_fid,['             Summary of CLImAT-HET results (version ' current_version ')          \n']);
            if NoSolutionFlag
                fprintf(o_fid,'Warning: Prediction results may be inaccurate due to the failure\n');
                fprintf(o_fid,'in finding optimal initial global parameters!\n');
            end
            
            fprintf(o_fid,'General information of this cancer sample:                      \n');
            fprintf(o_fid,'   Proportion of abnormal cells in the sample: %6.4f\n',beta(end));
            fprintf(o_fid,'   Subclonal populations in the sample:');
            for i = 1:length(beta)
                fprintf(o_fid,' %6.4f',beta(i));
            end
            fprintf(o_fid,'\n');
            fprintf(o_fid,'   Frequencies of subclonal population in the sample:');
            for i = 1:length(beta)
                fprintf(o_fid,' %6.4f',alpha(i+1));
            end
            fprintf(o_fid,'\n');
            fprintf(o_fid,'   Average copy number: %1.2f\n',ACN);
            fprintf(o_fid,'   Copy neutral read counts: %3.0f\n',lambda);
            fprintf(o_fid,'   Parameter p of NB distributions: ');
            for i = 1:length(p)
                fprintf(o_fid,' %6.4f',p(i));
            end
            fprintf(o_fid,'\n');
            
            fprintf(o_fid,'---------------------------------------------------------------\n');
            fprintf(o_fid,'\n');          
            
            %--------------output state assignment in segments--------------
            fprintf(o_fid,'Chr\tStartPos\tEndPos\tCN\tmCN\tAI\tTumor fraction\tScore\n');
            for i = 1:size(segments,1)
                chr_indx = segments(i,1);
                s_indx = segments(i,2);
                e_indx = segments(i,3);
                
                s_pos = data_pos_sep{chr_indx}(s_indx);
                e_pos = data_pos_sep{chr_indx}(e_indx);
                if sp_segs_all(i) == 0
                    tumor_fraction = 0;
                else
                    tumor_fraction = beta(sp_segs_all(i));
                end
                fprintf(o_fid,'%d\t%d\t%d\t%d\t%d\t%d\t%f\t%3.1f\n',Chromosomes(chr_indx),...
                            s_pos,e_pos,cn_segs_all(i),mcn_segs_all(i),AI_segs_all(i),tumor_fraction,scores(i));
            end
            fclose(o_fid);
            
            %6.4finall display a brief report
            t = toc;
            t_all = t_all+t;
            disp ([num2str(fileindx) '. ' filename{fileindx} ' is done, time used: ' num2str(t/60) ' minites']);
            
        end %for fileindx = 3:length(filename)
    end %if length(datafilelist) < 3
    
    disp(['-----CLImAT-HET batch annotation is finished, totally ' num2str(t_all/60) ' minites were used-----']);
    clear all
end %if sourceflag  == 1

clear all
close all

end
