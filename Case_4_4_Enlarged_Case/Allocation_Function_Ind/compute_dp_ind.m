function [dp_values]=compute_dp_ind(matrix)
 [num_rows, num_cols] = size(matrix);
 num_mem=num_cols-1;
 Vn=matrix(num_rows,num_cols);%the total profit
 Vi=zeros(num_mem,1);
 Vn_i=zeros(num_mem,1);
 xi=zeros(num_mem,1);
 for i=1:num_mem
     Vi(num_mem+1-i)=matrix(i,num_cols);
     Vn_i(num_mem+1-i)=matrix(num_rows-i,num_cols);
 end
upper=(num_mem-1)*Vn-sum(Vn_i,"all");
lower=Vn-sum(Vi,"all");
DP=upper/lower;%the calculation of propensity to disrupt
for i=1:num_mem
    xi(1+num_mem-i)=1/(1+DP)*(Vn-Vn_i(1+num_mem-i))+DP/(1+DP)*Vi(1+num_mem-i);
end
dp_values=xi';
end