function allocated_profit=compute_dp_two_stage(DP_DA,Real_Matrix,Marginal_Matrix)
    [num_rows, num_cols] = size(Real_Matrix);
    num_mem = num_cols - 1;
    VRN = Real_Matrix(num_rows,num_cols);

    % idx = (1:num_mem).';
    % VRi = Real_Matrix(num_mem + 1 - idx, num_cols);
    % VMari = Marginal_Matrix(num_mem + 1 - idx, num_cols);
    % VRn_i = Real_Matrix(num_rows - idx, num_cols);
    % VMarn_i = Marginal_Matrix(num_rows - idx, num_cols);
    for i=1:num_mem
     VRi(num_mem+1-i)=Real_Matrix(i,num_cols);
     VRn_i(num_mem+1-i)=Real_Matrix(num_rows-i,num_cols);
     VMari(num_mem+1-i)=Marginal_Matrix(i,num_cols);
     VMarn_i(num_mem+1-i)=Marginal_Matrix(num_rows-i,num_cols);
    end
    Smar = VRN - VMarn_i - VMari;
    x_initial = (1/(1+DP_DA)) * (VRN-VRn_i) + (DP_DA/(1+DP_DA)) * VRi;
    None_Allocation = VRN - sum(x_initial);
    total_smar = sum(Smar);
    if total_smar == 0
        adjust_allocation = zeros(num_mem,1);
    else
        adjust_allocation = None_Allocation * (Smar / total_smar);
    end
    allocated_profit = adjust_allocation + x_initial;
    allocated_profit = reshape(allocated_profit,[],1);
end
