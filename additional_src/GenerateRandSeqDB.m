function [randSeqDB]=GenerateRandSeqDB(trueSeqDB)
randSeqDB=rand(size(trueSeqDB))*4;
randSeqDB=floor(randSeqDB)+1;
randSeqDB=uint8(randSeqDB);