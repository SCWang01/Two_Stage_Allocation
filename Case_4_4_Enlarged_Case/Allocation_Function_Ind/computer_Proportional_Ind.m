% From the paper:Cost Allocation for Inertia and  Frequency Response Ancillary Services 
function Proportional_Value=computer_Proportional_Ind(matrix)
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
     for i=1:num_mem
         xi(i)=Vi(i)*Vn/(sum(Vi));
     end
     Proportional_Value=xi';
end

