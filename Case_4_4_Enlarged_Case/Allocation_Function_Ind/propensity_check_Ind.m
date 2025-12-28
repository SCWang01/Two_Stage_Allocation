function [chart]=propensity_check(matrix,method)
    [num_rows, num_cols] = size(matrix);
    num_mem=num_cols-1;
    Vi=zeros(num_mem,1);
    Vn_i=zeros(num_mem,1);
    vn=matrix(num_rows,num_cols);
    %x_i 计算
    if num_rows==2^num_mem
        for i=1:num_cols-1
            xi(i)=method(i);
        end
        %v_i 计算  V_(n\i)计算
        for i=1:num_cols-1
            Vi(num_mem+1-i)=matrix(2^(i-1)+1,num_cols);
            Vn_i(num_mem+1-i)=matrix(2^num_mem-2^(i-1),num_cols);
        end
        %x_(n\i) 计算
        for i=1:num_cols-1
            x_ni(i)=sum(method,"all")-xi(i);
        end
        %DP 计算
        for i=1:num_mem
            DP(i)=(x_ni(i)-Vn_i(i))/(xi(i)-Vi(i));
        end
        chart=[xi;Vi';x_ni;Vn_i';DP];
        %生成矩阵
    end
    if num_rows==2*num_mem+1
        for i=1:num_cols-1
            xi(i)=method(i);
        end
        %v_i 计算  V_(n\i)计算
        for i=1:num_mem
            Vi(num_mem+1-i)=matrix(i,num_cols);
            Vn_i(num_mem+1-i)=matrix(num_rows-i,num_cols);
        end
        %x_(n\i) 计算
        for i=1:num_mem
            x_ni(i)=sum(method,"all")-xi(i);
        end
        %DP 计算
        for i=1:num_mem
            DP(i)=(x_ni(i)-Vn_i(i))/(xi(i)-Vi(i));
        end
        chart=[xi;Vi';x_ni;Vn_i';DP];
    end
end