addpath("Test_Solution\");
addpath("Allocation_Function_Ind\");
DP_Scenario=zeros(18,1);
Base_matrix=ind_listg1(50);
f_vec = cellfun(@(s) s.f, Day_ahead_result);
Profit_matrix_DA=[Base_matrix,f_vec];
Added_O=Operator_2_Ind(Profit_matrix_DA);
Expected_Allocation=compute_dp_ind(Added_O);
DP_Expected=propensity_check_Ind(Added_O,Expected_Allocation);
DP_ahead=DP_Expected(end,end);
DP_Two=zeros(51,18);
Profit_Two=zeros(51,18);
DP_One=zeros(51,18);
Ind_Two=zeros(51,18);
for i=1:18
    f_vec = cellfun(@(s) s.f, Real_Time_result(:,i)); % real-time generation profit
    Profit_matrix=[Base_matrix,f_vec]; Added_matrix=Operator_2_Ind(Profit_matrix); % Generate the real-time profit
    mar_vec = cellfun(@(s) s.f, Marginal_result(:,i)); % marginal generation profit
    Profit_mar_matrix=[Base_matrix,mar_vec]; Added_matrix_mar=Operator_2_Ind(Profit_mar_matrix); % Generate the marginaltime profit
    Test11=compute_dp_ind(Profit_matrix);
    DP_Solution=compute_dp_ind(Added_matrix);
    Check_DP=propensity_check_Ind(Added_matrix,DP_Solution);
    DP_One(:,i)=Check_DP(end,:)';
    Two_Solution=compute_dp_two_stage(DP_ahead,Added_matrix,Added_matrix_mar);
    Check_DP_Two=propensity_check_Ind(Added_matrix,Two_Solution);
    DP_Two(:,i)=Check_DP_Two(end,:)';
    Profit_Two(:,i)=Check_DP_Two(1,:)';
    Ind_Two(:,i)=Check_DP_Two(2,:)';
end




