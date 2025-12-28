addpath("Allocation_Function\");
Base_matrix=bin_listg1(10);
DP_Shapley=zeros(11,18);
tic
for i=1:18
    f_vec = cellfun(@(s) s.f, Real_Time_result(:,i)); % real-time generation profit
    Profit_matrix=[Base_matrix,f_vec]; Added_matrix=Operator_2(Profit_matrix); % Generate the real-time profit
    Allocation_Solution=compute_shapley(Added_matrix);
    Check_Shapley=propensity_check(Added_matrix,Allocation_Solution);
    DP_Shapley(:,i)=Check_Shapley(end,:)';
end
Time_Shapley=toc;