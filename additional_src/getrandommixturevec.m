% Get the number of creatures and the number of actual creatures and
% generate an alpha vec according to some distribution
function [alpha_vec CreaturesInds] = GetRandomMixtureVec(M,I)

% Here randomize alphas according to some distribution .. Start with uniform:
%% alpha_vec = rand(1,I);


% Power law
beta = -5;
alpha_vec = [1:I] .^ (-beta); 
alpha_vec = alpha_vec ./ sum(alpha_vec); % Normalize

% Last step: Choose which creatures to include:
P = randperm(M);
CreaturesInds = P(1:I);


%% temp_alpha_vec = zeros(1,M);
%% temp_alpha_vec(P(1:I)) = alpha_vec;
%% alpha_vec = temp_alpha_vec;


