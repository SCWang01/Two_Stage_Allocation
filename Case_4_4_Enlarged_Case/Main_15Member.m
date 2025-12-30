clear all
tAll = tic;
addpath("Uncertain_Scenario\");
addpath("Large_Case_Optimization_Function\");
addpath("Optimization_Function\");
addpath("Allocation_Function\");
%%
VPP_List=[4,4,1,2,4];
Parameter=data_definition_large(0,VPP_List);
fg=ind_listg1(15);
nScen = size(Parameter.Lmpda,2);


nI = size(fg,1);
Day_ahead_result = cell(nI,1);
Real_Time_result = cell(nI,nScen);
Marginal_result=cell(nI,nScen);
Time_record = zeros(nI,1);
Time_record_2= zeros(nI,1);

for i = 1:nI
    tStart = tic;
    fn = fg(i,:);
    fn=fn';
    da = Large_Da_Considered_Demand(fn, Parameter);
    Day_ahead_result{i,1} = da;
    Pda = da.Pda;
    rtRow = cell(1,nScen);
    for j = 1:nScen
        rtRow{1,j} = Large_Rt_Considered_Demand(fn, Parameter, j, Pda);
    end
    Real_Time_result(i,:) = rtRow;
    Time_record(i) = toc(tStart);
end

Findal_Da= Day_ahead_result{31,1}; 
w_number=Parameter.w_number;%number of WPPs 
p_number=Parameter.p_number;%number of PVs
tg_number=Parameter.tg_number;%number of TGs
ess_number=Parameter.ess_number;%number of ESSs
dr_number=Parameter.dr_number;%number of DRs
Pwppda=Findal_Da.Pwppda;Ppvpda=Findal_Da.Ppvpda;Ptgda=Findal_Da.Ptgda;Pdeda=Findal_Da.Pdeda;
Pceda=Findal_Da.Pceda;Cut_demand=Findal_Da.Cut_demand;
for i=1:nI
    tStart_2 = tic;
    fn=fg(i,:).';
    idx=0;
    fn_w  = fn(idx+1:idx+w_number);        idx = idx +w_number;
    fn_p  = fn(idx+1:idx+p_number);        idx = idx +p_number;
    fn_tg = fn(idx+1:idx+tg_number);       idx = idx +tg_number;
    fn_ess= fn(idx+1:idx+ess_number);      idx = idx +ess_number;
    fn_dr = fn(idx+1:idx+dr_number);
    Pda_marginal= Pwppda * fn_w + ...
    Ppvpda * fn_p + ...
    Ptgda * fn_tg + ...
    (Pdeda - Pceda) * fn_ess + ...
    Cut_demand * fn_dr;
    for j=1:nScen
            Marginal_result{i,j}=Large_Rt_Internal(fn,Parameter,j,Pda_marginal);
    end
    Time_record_2(i) = toc(tStart_2);
end
TotalTime = toc(tAll);
disp(TotalTime);
% Pda=Optimization_test_day_ahead.Pda;
% K=1;
% OPtimization_test_real_time=Large_Rt_Considered_Demand(fn,Parameter,1,Pda);
% Optimization_test_internal=Large_Rt_Internal(fn,Parameter,1,Pda);
    
