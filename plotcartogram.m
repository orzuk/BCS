% Plots a 4*L matrix as a 'cartogram' 
% by simply plotting the 4 curves in different colors
%
function Dummy = PlotCartogram(PSSM)

color_vec = 'brgk';

figure; hold on; xlabel('Nucleotide'); ylabel('Intensity');
for i=1:4
    plot(PSSM(i,:), [color_vec(i)]);
end
legend('A', 'C', 'G', 'T');
