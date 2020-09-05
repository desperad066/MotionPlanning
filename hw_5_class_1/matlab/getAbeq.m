function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    %
    %
    %
    %
    Aeq_start(1,1) = 1; %p
    Aeq_start(2,2) = 1; %v
    Aeq_start(3,3) = 1*2; %a
    Aeq_start(4,4) = 1*2*3; %j
    beq_start = start_cond';
    
    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    %
    %
    %
    T = ts(end);
    Aeq_end(1,end-(n_order+1)+1:end) = [1, T, T^2, T^3, T^4, T^5, T^6, T^7]; %p
    Aeq_end(2,end-(n_order+1)+1:end) = [0, 1, 2*T, 3*T^2, 4*T^3, 5*T^4, 6*T^5, 7*T^6]; %v
    Aeq_end(3,end-(n_order+1)+1:end) = [0, 0, 2, 6*T, 12*T^2, 20*T^3, 30*T^4, 42*T^5]; %a
    Aeq_end(4,end-(n_order+1)+1:end) = [0, 0, 0, 6, 24*T, 60*T^2, 120*T^3, 210*T^4]; %j
    beq_end = end_cond';
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    %
    %
    %
    
    if n_seg > 1
        for i = 1: n_seg-1
            T = ts(i);
            Aeq_wp(i,(i-1)*(n_order+1)+1 : i*(n_order+1)) = [1, T, T^2, T^3, T^4, T^5, T^6, T^7];
            beq_wp(i) = waypoints(1+i);
        end
    end
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    %
    %
    %
    if n_seg > 1
        for i = 1: n_seg-1
            T = ts(i);
            Aeq_con_p(i,(i-1)*(n_order+1)+1 : i*(n_order+1)) = [1, T, T^2, T^3, T^4, T^5, T^6, T^7];
            Aeq_con_p(i,(i)*(n_order+1)+1 : (i+1)*(n_order+1)) = -[1, 0, 0, 0, 0, 0, 0, 0];
        end
        beq_con_p = zeros(n_seg-1, 1);
    end
    
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %
    if n_seg > 1
        for i = 1: n_seg-1
            T = ts(i);
            Aeq_con_v(i,(i-1)*(n_order+1)+1 : i*(n_order+1)) = [0, 1, 2*T, 3*T^2, 4*T^3, 5*T^4, 6*T^5, 7*T^6];
            Aeq_con_v(i,(i)*(n_order+1)+1 : (i+1)*(n_order+1)) = -[0, 1, 0, 0, 0, 0, 0, 0];
        end
        beq_con_v = zeros(n_seg-1, 1);
    end
    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %
    if n_seg > 1
        for i = 1: n_seg-1
            T = ts(i);
            Aeq_con_a(i,(i-1)*(n_order+1)+1 : i*(n_order+1)) = [0, 0, 2, 6*T, 12*T^2, 20*T^3, 30*T^4, 42*T^5];
            Aeq_con_a(i,(i)*(n_order+1)+1 : (i+1)*(n_order+1)) = -[0, 0, 1*2, 0, 0, 0, 0, 0];
        end
        beq_con_a = zeros(n_seg-1, 1);
    end
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    if n_seg > 1
        for i = 1: n_seg-1
            T = ts(i);
            Aeq_con_j(i,(i-1)*(n_order+1)+1 : i*(n_order+1)) = [0, 0, 0, 6, 24*T, 60*T^2, 120*T^3, 210*T^4];
            Aeq_con_j(i,(i)*(n_order+1)+1 : (i+1)*(n_order+1)) = -[0, 0, 0, 1*2*3, 0, 0, 0, 0];
        end
        beq_con_j = zeros(n_seg-1, 1);
    end
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end