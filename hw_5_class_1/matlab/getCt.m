function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    %
    %
    %
    %
    Ct = zeros(n_seg*8, 8+4*(n_seg-1));
    
    I = eye(4);
    
    Ct(1:4, 1:4) = I;
    for i = 1:n_seg-1
        Ct(4+8*(i-1)+1, 8+(i-1)+1) = 1;
        Ct(4+8*(i-1)+2: 4+8*(i-1)+4, 8+n_seg-1+3*(i-1)+1: 8+n_seg-1+3*(i-1)+3) = eye(3);
        
        Ct(4+8*(i-1)+5, 8+(i-1)+1) = 1;
        Ct(4+8*(i-1)+6: 4+8*(i-1)+8, 8+n_seg-1+3*(i-1)+1: 8+n_seg-1+3*(i-1)+3) = eye(3);
    end
    Ct(end-3:end, 5:8) = I;
    
end