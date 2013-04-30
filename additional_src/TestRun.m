% find the closest bacteria to each result
% seqred2b - the sequence database used for the reconstruction (non-binned)
% runres - the (sparse) frequencies of the reconstruction
%
% Example:
% TestRun(Rev-Red-seq, freqs_nvh, start_pos,end_pos)
%
function TestRun(SeqRed2B,runres,startpos,endpos, save_flag, output_file)
compset=[7637,7431,5259,4808,4583,4773];
[sres,si]=sort(runres);
for a=length(si)-20:length(si)
    [minspec,minrespos,minres] = ...
        CompareToTruth(si(a),SeqRed2B,compset,startpos,endpos);
    disp(['pos ' num2str(length(si)-a) ' Freq ' num2str(sres(a)) ...
        ' MinDist ' num2str(minres) ' MinSpecies ' num2str(minspec) ...
        ' minspecpos ' num2str(minrespos)]);
    if(save_flag)
        save(output_file,   'minspec','minrespos','minres', 'si', 'sres', 'a');
    end
end

