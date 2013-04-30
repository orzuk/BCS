% Generate a random species mixture and calculate the PSSM
%
% input:
% R_seqInt - The sequence database
% sigma - the standard deviation of the normal noise added to the mixture
% I - the desired number of species in the mixture
% distflag - 0 is uniform, 1 is powerlaw
%
% output:
% orig_ind - Indices of the species present in the mixture
% orig_alpha - their frequencies
% Y - The mixture PSSM
function [orig_ind orig_alpha Y]=SimulateMetagenomeMixture(R_seqInt,sigma,I,distFlag)
N = size(R_seqInt,2); % Length of sequences

%I = 300; % Number of actuall creatures
%sigma = 0.0; % The amount of measurement noise

% First generate the sequences
A = 1; C = 2; G = 3; T = 4; U = T;

%data_flag = 1; % 1 real data, 0 rand data

% Notation: A=1, C=2, G=3, T=4


M=size(R_seqInt,1);
[alpha_vec CreaturesInds] = GetRandomMixtureVec2(M,I,distFlag); % Get the vector of random mixture
X = zeros(I,N);
for i=1:length(CreaturesInds)
    cs=R_seqInt(CreaturesInds(i),:);
    % find the unknown nucleotides and generate a random fill for them
    Unknowns=find(cs==5);
    for a=1:length(Unknowns)
        nn=randperm(4);
        cs(Unknowns(a))=nn(1);
    end
    X(i,:)=cs;
end

Y = zeros(4,N);
Y(1,:) = alpha_vec * (X == A);
Y(2,:) = alpha_vec * (X == C);
Y(3,:) = alpha_vec * (X == G);
Y(4,:) = alpha_vec * (X == T);

%% Y = alpha_vec * X(CreaturesInds,:); % Y is what we observe



% Add some noise
r = randn(4,N);
% r=r-mean(r);
r = sigma * r*sum(alpha_vec)/4;
%./ norm(r);
Y = Y + r;

orig_alpha=alpha_vec;
orig_ind=CreaturesInds;
