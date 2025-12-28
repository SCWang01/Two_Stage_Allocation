function  Matrix_Considered_Operator=Operator_2_Ind(Matrix)
    Length=size(Matrix,2);
    Profit=Matrix(:,end);
    origin_num_mem=Length-1;
    added_num=origin_num_mem+1;
    new_matrix=ind_listg1(added_num);
    new_matrix=[new_matrix,zeros(2*added_num+1,1)];
    new_matrix(1:origin_num_mem,end)=Profit(1:origin_num_mem);
    new_matrix(added_num,end)=0;
    new_matrix(added_num+1,end)=sum(Profit(1:origin_num_mem));
    new_matrix(added_num+2:end,end)=Profit(origin_num_mem+1:end);
    Matrix_Considered_Operator=new_matrix;
end