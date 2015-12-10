function transmat = CLImATHET_init_transmat(data_pos, state_num)

L = 2000000;
N = length(data_pos);
transmat = zeros(state_num, state_num, N-1);
rho = 0.5*(1-exp(-(data_pos(2:end)-data_pos(1:end-1))/(2*L)));

for i = 1:N-1
    transmat(:,:,i) = (ones(state_num)-eye(state_num))*rho(i)/(state_num-1)+eye(state_num)*(1-rho(i));
end

end