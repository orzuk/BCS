% Compute a binary basis spreading all sequences.
% Assume no sequence is zero
function B = GetBinaryBasis(BinSeqs)
TOL = 0.000000000000000001;
NumSeqs = size(BinSeqs, 1)
SeqsLen = size(BinSeqs, 2)

B = BinSeqs(1,:);
max_b_ind = max(find(B)); 

for i=2:NumSeqs
    V = BinSeqs(i,:);
    base_vec = size(B,1)+1
    while(max(V) > TOL)  
        max_v_ind = max(find(V));
        base_vec = find(B(1:base_vec-1,max_v_ind)); % Only one can be ..
        if(isempty(base_vec))
            break;
        end
        V = V - B(base_vec,:) .* (V(max_v_ind) / B(base_vec, max_v_ind)); 
    end
    
    if(max(V) > TOL) % Didn't find it ... 
        for j=?:-1:1
             max_v_ind = max(find(V));
             V = B(j,:) .* (V(max_v_ind) / B(base_vec, max_v_ind));
        V = V ./ V(max_v_ind);
        B = [B' V']'; B=sortrows(B(:,end:-1:1)); B = B(:,end:-1:1);
    end


end