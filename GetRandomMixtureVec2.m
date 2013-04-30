% Get the number of creatures in the mixture and the total number creatures
% in database, % generate a frequency vector according to some distribution
%
% input:
% M - number of species in database
% I - desired number of species in mixture
% distFlag - frequency distribution - 0 is uniform, 1 is powerLaw
%
% output:
% alpha_vec - the frequency of each creature present in the mixture
% CreturesInds - their indices from the database
%
function [alpha_vec SpeciesInds] = GetRandomMixtureVec2(M,I,distFlag)

% Here randomize alphas according to some distribution ...

% Uniform
if (distFlag==0)
    alpha_vec = rand(1,I);
end

% Power law
if (distFlag==1)
    % note changed 8.3.2010 to powerlaw from beta=-5
    beta = 1;
    alpha_vec = [1:I] .^ (-beta);
end

alpha_vec = alpha_vec ./ sum(alpha_vec); % Normalize

% Last step: Choose which creatures to include:
P = randperm(M);
SpeciesInds = P(1:I);


%% temp_alpha_vec = zeros(1,M);
%% temp_alpha_vec(P(1:I)) = alpha_vec;
%% alpha_vec = temp_alpha_vec;


