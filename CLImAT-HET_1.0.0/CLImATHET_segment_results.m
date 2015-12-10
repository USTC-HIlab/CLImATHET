function segments = CLImATHET_segment_results(state_seq,sp_seq)
% 04/02/2014 by Zhenhua
% This function is used to generate segments from a sequence of assigned states

pre_state = [];
pre_sp = [];
s_indx = [];
segments = [];
for i = 1:length(state_seq)
    if isempty(pre_state)
        pre_state = state_seq(i);
        pre_sp = sp_seq(i);
        s_indx = i;
    elseif state_seq(i) ~= pre_state || sp_seq(i) ~= pre_sp
        segments = [segments;s_indx i-1 pre_state pre_sp];
        pre_state = state_seq(i);
        pre_sp = sp_seq(i);
        s_indx = i;
    end
end

if s_indx <= length(state_seq)
    segments = [segments;s_indx length(state_seq) pre_state pre_sp];
end

end