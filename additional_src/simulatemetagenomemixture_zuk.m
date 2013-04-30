data_flag = 0; % 1 real data, 0 rand data
if(data_flag == 0)
    M = 500; % Number of creatures    
else
    load('current_prokMSA_unaligned_seq.mat');
    M = length(R_seq);
end

N = 180; % Length of sequence

I = 30; % Number of actuall creatures
sigma = 0.0; % The amount of measurement noise

% First generate the sequences
A = 1; C = 2; G = 3; T = 4; U = T;


% Notation: A=1, C=2, G=3, T=4

[alpha_vec CreaturesInds] = GetRandomMixtureVec(M,I); % Get the vector of random mixture

if(data_flag == 0)
    X = floor(rand(M,N)*4)+1;
    X = X(CreaturesInds,:);
else
    X = zeros(I,N);
    for i=1:length(CreaturesInds)
        Temp = convert_nuc_to_num(R_seq{CreaturesInds(i)});
        X(i,1:length(Temp)) = Temp; 
    end
end

Y = zeros(4,N);
Y(1,:) = alpha_vec * (X == A);
Y(2,:) = alpha_vec * (X == C);
Y(3,:) = alpha_vec * (X == G);
Y(4,:) = alpha_vec * (X == T);

%% Y = alpha_vec * X(CreaturesInds,:); % Y is what we observe



% Add some noise
%% r = randn(4,N); r=r-mean(r);
%% r = sigma * r ./ norm(r);
%% Y = Y + r;

save 'ExampleSample.mat' alpha_vec CreaturesInds Y

return;

alpha_vec_estimate = ReconstructMixture(Y, R_seq, I, 'greedy'); % Call our reconstruction algorithm (currently greedy)

% Get some measures of how well we did
ReconstructionError = norm(alpha_vec - alpha_vec_estimate)
figure; plot(alpha_vec, alpha_vec_estimate, '.'); title('True freqs. and our estimation');
xlabel('True Freqs'); ylabel('Our Estimation');















