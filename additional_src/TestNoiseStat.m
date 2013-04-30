function [tres]=TestNoiseStat(predbinchrom,predbinchromall,startpos,endpos,tau,noise)
tres=[];
for a=1:10
    disp(a);
    [res]=TestNoise(predbinchrom,predbinchromall,startpos,endpos,tau,noise);
    save(['TestNoise-' num2str(a)],'res');
    tres=[tres;res];
end