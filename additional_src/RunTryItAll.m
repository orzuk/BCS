in_matlab_flag = 1; % 1 - run interactively 0 - run jobs 
reconstruct_random = 0; % 0: reconstruct a random mixture. 1: reconstruct the true mixture
sqrt_database_flag = 1;
if(sqrt_database_flag)
    tau_vec = round(sqrt([1000 5000 10000 15000 20000]));
    sqrt_str = '_sqrt_';
else
    tau_vec = [1000 5000 10000 15000 20000];
    sqrt_str = '';
end
for tau=tau_vec
    w_vec = tau .* [0:20]; % w should be on the same order of magnitude as tau but probably a bit bigger
    for w=w_vec
        eval_str = ['TryItAll(' num2str(w) ', ' ...
            num2str(tau) ', ' num2str(reconstruct_random) ', ' ...
            num2str(sqrt_database_flag) ')'];
        if(in_matlab_flag)
           eval(eval_str); 
        else
        SubmitMatlabJobToFarm(eval_str, ...
            ['out/new_try_metagenomics_gpsr_w_' num2str(w) '_tau_' ...
            num2str(tau) sqrt_str '.out'], 'broad');
        end
    end
end
