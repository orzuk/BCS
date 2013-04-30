% Try reconstruction of a random or true mixture
function TryItAll(w_vec, tau, reconstruct_random, sqrt_database_flag)

% Try the binned reconstruction
%                       seq DB  , exp mix ,start,stop,tau
AssignGeneralConstants;
[machine machine_delim html_outdir] = get_machine_type();% machine = UNIX; % UNIX;
new_data_flag = 2; % use new preprocessing - should be cleaner. 2 means new 2010 data
if(~exist('reconstruct_random', 'var') || isempty(reconstruct_random))
    reconstruct_random = 0; % try a random reconstruction
end
if(~exist('sqrt_database_flag', 'var') || isempty(sqrt_database_flag))
    sqrt_database_flag = 1;
end
if(sqrt_database_flag)
    sqrt_str = '_sqrt_';
else
    sqrt_str = '';
end

broad = 0;
switch machine
    case UNIX
        data_dir = '/seq/orzuk/zuk_weizmann/research/compressed_sensing/metagenomics/data';
        %         load('/seq/orzuk/zuk_weizmann/research/compressed_sensing/metagenomics/data/data.mat'); % load RedName, predbinSqr, mixnvhbin
        %         new_vec = load('/seq/orzuk/zuk_weizmann/research/compressed_sensing/metagenomics/data/data.mat'); % load RedName, predbinSqr, mixnvhbin
        %        load('~/public_html/matlab/libs/metagenomics/data.mat'); % load RedName, predbinSqr, mixnvhbin
    case PC
        if(broad)
            data_dir = 'y:\public_html\matlab\libs\metagenomics';
            %             load('y:\public_html\matlab\libs\metagenomics\data.mat'); % load RedName, predbinSqr, mixnvhbin
            %             new_vec = load('c:/research/metagenomics/data/new_data/data.mat');
        else
            data_dir = 'c:/research/metagenomics/data';
            %             load('c:/research/metagenomics/data/data.mat');
            %             new_vec = load('c:/research/metagenomics/data/new_data/data.mat');
        end
end
if(sqrt_database_flag) % for sqrt we also need to change tau!!!!
    new_vec = load(fullfile(data_dir, 'binchromnvhsqrt.mat'));
    load(fullfile(data_dir, 'PredBinSqrt.mat')); % load data. Var name is: npred
else
    load(fullfile(data_dir, 'data.mat')); % load data. Var. name is: predbinSqr
end
new_vec = load(fullfile(data_dir, 'new_data', 'data.mat')); % load new data
new_vec_2010 = load(fullfile(data_dir, 'Mixnvh-12-1-2010.mat'));

if(~exist('tau', 'var') || isempty(tau))
    tau = 10000; % GPSR parameter
end
if(~exist('w', 'var') || isempty(w))
    w = 0; % sparisty parameter for slack variables
end

if(new_data_flag)
    mixallbin = new_vec.mixallbin;
    mixn28bin = new_vec.mixn28bin;
end
switch new_data_flag
    case 1
        mixnvhbin = new_vec.mixnvhbin;
    case {2, 2010} % new data - for now works only on nvh
        mixnvhbin = new_vec_2010.binchromnvhNew;
end

pos_start_vec = [100 150 200 300 400];
pos_end_vec = [490 500 600 700 800]; %  900 1000]; % must get different values for different regions
for pos_start = vec2row(pos_start_vec)
    for pos_end = vec2row(pos_end_vec)
        try_start = pos_start
        try_end = pos_end
        L = pos_end - pos_start + 1
        if(reconstruct_random) % Generate random mixture of five species (equal frequencies)
            N = 1000;
            predbinSqr = predbinSqr(1:N,:,:); % take a smaller matrix - easier computation
            N = size(predbinSqr, 1); L = size(predbinSqr, 2);
            
            num_bad_equations = 10; % set number of bad equations
            bad_equations_inds = randperm(N);
            bad_equations_inds = bad_equations_inds(1:num_bad_equations);
            
            rand_inds = randperm(N); rand_inds = rand_inds(1:5);
            rand_weights = rand(5,1); rand_weights = rand_weights ./ sum(rand_weights);
            rand_mixture = reshape(sum(predbinSqr(rand_inds,:,1:4)), L, 4); % .* rand_weights;
            rand_mixture(bad_equations_inds,:) = rand(num_bad_equations,4);
            rand_mixture = rand_mixture ./ max(epsilon,repmat(sum(rand_mixture,2), 1, 4));
            
            bad_equations_inds_display = bad_equations_inds - pos_start+1;
            bad_equations_inds_display = [bad_equations_inds_display ...
                bad_equations_inds_display + L ...
                bad_equations_inds_display + 2*L ...
                bad_equations_inds_display + 3*L]
            
            for w=vec2row(w_vec)
                [rand_freqs slack_vars] = ...
                    TryGPSRBinAlign(predbinSqr,rand_mixture,pos_start,pos_end,tau,w); % try a random mixture
                rand_inds_reconstructed = find(rand_freqs)'
                rand_inds_original = rand_inds
                ignored_eqs = vec2row(find(slack_vars))
                save(['try_metagenomics_gpsr_rand_w_' num2str(w) '_tau_' ...
                    num2str(tau) '_pos_' ...
                    num2str(pos_start) '_' num2str(pos_end) '.mat'], ...
                    'rand_freqs', 'w', 'rand_inds_original', 'rand_inds_reconstructed', ...
                    'bad_equations_inds_display', 'slack_vars', 'ignored_eqs');
            end
        else % here do the real mixture
            %    freqs = [];
            %    w_vec = 2.^[-10:10]; % try different w's
            %    for i=1:length(w_vec)
            for w=vec2row(w_vec)
                if(exist('predbinSqr', 'var'))
                    freqs_nvh=TryGPSRBinAlign(predbinSqr,mixnvhbin,pos_start,pos_end,tau,w); % try the true mixture
                else
                    freqs_nvh=TryGPSRBinAlign(npred,mixnvhbin,pos_start,pos_end,tau,w); % try the true mixture
                end
                [sres_nvh,si_nvh]=sort(freqs_nvh); % sort according to frequency
                sres_nvh_top = sres_nvh(end-10:end) % view 10 most prominent species
                top_names_nvh = RedName{si_nvh(end-9:end)}; % and their names
                RedName{si_nvh(end-9:end)}
                
                save(['try_metagenomics_gpsr_w_' num2str(w) '_tau_' num2str(tau) '_pos_' ...
                    num2str(pos_start) '_' num2str(pos_end) sqrt_str '.mat'], ...
                    'freqs_nvh', 'w', 'sres_nvh', ...
                    'si_nvh', 'top_names_nvh');
                
                if(new_data_flag == 1) % run on all here
                    freqs_all = TryGPSRBinAlign(predbinSqr,mixallbin, ...
                        pos_start,pos_end,tau,w); % try the true mixture
                    [sres_all,si_all]=sort(freqs_all); % sort according to frequency
                    sres_all_top = sres_all(end-10:end) % view 10 most prominent species
                    top_names_all = RedName{si_all(end-9:end)}; % and their names
                    RedName{si_all(end-9:end)}
                    
                    save(['try_metagenomics_gpsr_w_' num2str(w) '_tau_' ...
                        num2str(tau) '_pos_' ...
                        num2str(pos_start) '_' num2str(pos_end) '.mat'], ...
                        'freqs_all', 'sres_all', 'top_names_all', '-append');
                end
                
            end
        end
    end
end

%end
% We expect to see:
% Photobacter leignatti
% Vibrio Fischeri
% Escherichia coli
% staphylococcus epidermidis
% Enterococcus faecalis


