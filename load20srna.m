% One time script: convert database to matlab format

%% R=loadcellfile('current_prokMSA_unaligned.txt');

%% R_names = R(1:2:end);
%% R_seq = R(2:2:end);

load('current_prokMSA_unaligned_seq.mat');

cur_len=-1; min_len = 9999999;
for i=1:length(R_seq)
    if(mod(i,500) == 0)
        i_is = i
    end
    cur_len = max(cur_len,length(R_seq{i}));
    min_len = min(min_len,length(R_seq{i}));
end

R_seq_num = round(4*rand(length(R_seq), ceil(cur_len/16)));

for i=1:length(R_seq)
    if(mod(i,500) == 0)
        i_is = i
    end
    R_seq_num{i} = convert_nuc_to_num(R_seq{i});
end

save 'current_prokMSA_unaligned_seq_num.mat' R_seq_num 

