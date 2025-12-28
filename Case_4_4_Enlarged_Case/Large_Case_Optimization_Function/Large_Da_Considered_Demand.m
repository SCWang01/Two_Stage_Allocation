function [Result_dabe]=Large_Da_Considered_Demand(fn,Parameter)
    %% 变量设置
    %%%  here we consider the shape of the uncertain is 24* K (scenario number) * index
    %%% 1-K(scenario) Two-stage optimization
    ps=Parameter.P;
    w_number=Parameter.w_number;%number of WPPs 
    p_number=Parameter.p_number;%number of PVs
    tg_number=Parameter.tg_number;%number of TGs
    ess_number=Parameter.ess_number;%number of ESSs
    dr_number=Parameter.dr_number;%number of DRs

    %day-ahead quantities
    Pda=sdpvar(24,1);% day-ahead bid quantity_1
    %day-ahead wpp generation
    Pwppda=sdpvar(24,w_number);% wind power generation_1
    %day-ahead pv generation
    Ppvpda=sdpvar(24,p_number);%
    %day-ahead tg generation
    Ptgda=sdpvar(24,tg_number);%
    stg=binvar(24,tg_number);%
    %day-ahead ess charge/discharge
    Pceda=sdpvar(24,ess_number);%charging power
    Pdeda=sdpvar(24,ess_number);%discharging power
    soc=sdpvar(24,ess_number);%SOC variable_1
    cess=binvar(24,ess_number);%charging state change
    dess=binvar(24,ess_number);%discharging state change
    %day-ahead cutting load
    Cut_demand=sdpvar(24,dr_number);%curtailed load

    %cost of the VPP
    Cda=sdpvar(24,1);%day-ahead cost
    %cost of internal tgs
    Ctgda=sdpvar(24,tg_number);% tg day-ahead cost
    %cost of ESSs
    Cessda=sdpvar(24,ess_number);% ess cost
    %cost of DRs
    Cdrda=sdpvar(24,dr_number);%    curtailed demand cost
    
    n=24;
    
    lmpdas=Parameter.Lmpda; 
    
    lmprts=Parameter.Lmprt;
    
    K=size(Parameter.Wind_f(:,:,1),2);
    
    F=[];
    % index setting 
    idx = 0;
    fn_w  = fn(idx+1:idx+w_number);        idx = idx + w_number;
    fn_p  = fn(idx+1:idx+p_number);        idx = idx + p_number;
    fn_tg = fn(idx+1:idx+tg_number);       idx = idx + tg_number;
    fn_ess= fn(idx+1:idx+ess_number);      idx = idx + ess_number;
    fn_dr = fn(idx+1:idx+dr_number);
    has_w = any(fn_w ~= 0);
    has_p = any(fn_p ~= 0);
    has_tg = any(fn_tg ~= 0);
    has_ess = any(fn_ess ~= 0);
    has_dr = any(fn_dr ~= 0);
    has_w = any(fn_w ~= 0);
    has_p = any(fn_p ~= 0);
    has_tg = any(fn_tg ~= 0);
    has_ess = any(fn_ess ~= 0);
    has_dr = any(fn_dr ~= 0);
    %% Constraints
    %-----------------------------------------------------------------------------------------
    %day-ahead market bidding constraints
    F = [F, Pda(:,1) == ...
    Pwppda * fn_w + ...
    Ppvpda * fn_p + ...
    Ptgda * fn_tg + ...
    (Pdeda - Pceda) * fn_ess + ...
    Cut_demand * fn_dr ];


    %-----------------------------------------------------------------------------------------
    %WPP and PV output constraints
    if has_w
        Wind_min = squeeze(min(Parameter.Wind_f,[],2)); % 24 x w_number
        Wind_max = squeeze(max(Parameter.Wind_f,[],2));
        fn_w_ex = repmat(reshape(fn_w,1,[]), n, 1);
        F = [F, Wind_min .* fn_w_ex <= Pwppda <= Wind_max .* fn_w_ex];
    else
        F = [F, Pwppda == 0];
    end

    if has_p
        Solar_min = squeeze(min(Parameter.Solar_f,[],2)); % 24 x p_number
        Solar_max = squeeze(max(Parameter.Solar_f,[],2));
        fn_p_ex = repmat(reshape(fn_p,1,[]), n, 1);
        F = [F, Solar_min .* fn_p_ex <= Ppvpda <= Solar_max .* fn_p_ex];
    else
        F = [F, Ppvpda == 0];
    end
    

    %day-ahead constraints-tg
    %-----------------------------------------------------------------------------------------
    if has_tg
        Pgmin = reshape(Parameter.Pgmin,1,[]);
        Pgmax = reshape(Parameter.Pgmax,1,[]);
        Dg = reshape(Parameter.Dg,1,[]);
        Ug = reshape(Parameter.Ug,1,[]);
        Pgmin_ex = repmat(Pgmin, n, 1);
        Pgmax_ex = repmat(Pgmax, n, 1);
        fn_tg_ex = repmat(reshape(fn_tg,1,[]), n, 1);
        F = [F, Pgmin_ex .* fn_tg_ex .* stg <= Ptgda <= Pgmax_ex .* fn_tg_ex .* stg];
        F = [F, -repmat(Dg,n-1,1) <= diff(Ptgda) <= repmat(Ug,n-1,1)];
        F = [F, stg(1,:)==0];
    else
        F = [F, Ptgda == 0, stg == 0];
    end

    %day-ahead constraints-ESS
    %-----------------------------------------------------------------------------------------
    socstart = reshape(Parameter.socstart,1,[]);
    if has_ess
        p_cmax = reshape(Parameter.p_cmax,1,[]);
        p_dmax = reshape(Parameter.p_dmax,1,[]);
        socmin = reshape(Parameter.socmin,1,[]);
        socmax = reshape(Parameter.socmax,1,[]);
        yt_sin = reshape(Parameter.yt_sin,1,[]);
        yt_sout = reshape(Parameter.yt_sout,1,[]);
        Esmax = reshape(Parameter.Esmax,1,[]);
        p_cmax_ex = repmat(p_cmax, n, 1);
        p_dmax_ex = repmat(p_dmax, n, 1);
        fn_ess_ex = repmat(reshape(fn_ess,1,[]), n, 1);
        F = [F, 0 <= Pceda <= p_cmax_ex .* fn_ess_ex .* cess];
        F = [F, 0 <= Pdeda <= p_dmax_ex .* fn_ess_ex .* dess];
        F = [F, 0 <= cess + dess <= 1];
        socmin_ex = repmat(socmin, n, 1);
        socmax_ex = repmat(socmax, n, 1);
        yt_sin_ex = repmat(yt_sin, n, 1);
        yt_sout_ex = repmat(yt_sout, n, 1);
        Esmax_ex = repmat(Esmax, n, 1);
        sin_over_es = yt_sin_ex ./ Esmax_ex;
        inv_sout_es = 1 ./ (yt_sout_ex .* Esmax_ex);
        F = [F, socmin_ex <= soc <= socmax_ex];
        F = [F, soc(24,:) == socstart];
        F = [F, soc(1,:) == socstart + (Pceda(1,:).*sin_over_es(1,:) - Pdeda(1,:).*inv_sout_es(1,:))];
        F = [F, diff(soc) == (Pceda(2:end,:).*sin_over_es(2:end,:) - Pdeda(2:end,:).*inv_sout_es(2:end,:))];
    else
        F = [F, Pceda == 0, Pdeda == 0, cess == 0, dess == 0, soc == repmat(socstart,n,1)];
    end

    %curtailed constraints-dr
    if has_dr
        fn_dr_ex = repmat(reshape(fn_dr,1,[]), n, 1);
        F = [F, 0 <= Cut_demand <= 0.85 * Parameter.Load_base .* fn_dr_ex];
    else
        F = [F, Cut_demand == 0];
    end


    %-----------------------------------------------------------------------------------------

    %cost function 
    ag = reshape(Parameter.ag,1,[]);
    bg = reshape(Parameter.bg,1,[]);
    cg = reshape(Parameter.cg,1,[]);
    ag_ex = repmat(ag, n, 1);
    bg_ex = repmat(bg, n, 1);
    cg_ex = repmat(cg, n, 1);
    Ctgda = ag_ex .* Ptgda + cg_ex .* Ptgda.^2 + bg_ex .* stg;
    Cessda = repmat(reshape(Parameter.as,1,[]), n, 1) .* (Pceda + Pdeda);
    Cdrda = repmat(reshape(Parameter.Cutting,1,[]), n, 1) .* Cut_demand;

    Cda(:,1)=sum(Ctgda,2)+sum(Cessda,2)+sum(Cdrda,2);

    %%%day-head profit calculation 
    Profitda = sum(lmpdas .* repmat(Pda(:,1),1,K) - repmat(Cda(:,1),1,K), 1);

    

    %%%实时投标
    %%% real scenario market setting
    M=100000;
    mu_1=binvar(24,K);
    K=size(Parameter.Lmpda,2);
    Pre=sdpvar(24,K);%real-time generation bid quantitiy
    Preal=sdpvar(24,K);%real-time generation
    %-----------------------------------------------------------------------------------------
    %DER definition

    Pwppre=sdpvar(24,K,w_number);%wpp definition
    Ppvpre=sdpvar(24,K,p_number);%光伏出力-1

    Ptgre=sdpvar(24,K,tg_number);%TG出力-1
    stgre=binvar(24,K,tg_number);%TG real-time
    Ctgre=sdpvar(24,K,tg_number); %tg实时成本-1


    Pcere=sdpvar(24,K,ess_number);%充电电量-1
    Pdere=sdpvar(24,K,ess_number);%放电电量-1
    socre=sdpvar(24,K,ess_number);%实时SOC情况-1
    Re_cess=binvar(24,K,ess_number);%charging
    Re_dess=binvar(24,K,ess_number);%discharging

    Cutre_demand=sdpvar(24,K,dr_number);% DR-1
    Cdere=sdpvar(24,K,dr_number);%DR-1 cost
    
    Cre=sdpvar(24,K);%实时成̜�

    % Real_Time Generation and market setting
    idx = 0;
    fn_w  = fn(idx+1:idx+w_number);        idx = idx + w_number;
    fn_p  = fn(idx+1:idx+p_number);        idx = idx + p_number;
    fn_tg = fn(idx+1:idx+tg_number);       idx = idx + tg_number;
    fn_ess= fn(idx+1:idx+ess_number);      idx = idx + ess_number;
    fn_dr = fn(idx+1:idx+dr_number);

    % 24×K
    Pwpp_sum = zeros(n,K);
    Ppv_sum = zeros(n,K);
    Ptg_sum = zeros(n,K);
    Pess_sum = zeros(n,K);
    Pdr_sum = zeros(n,K);
    if has_w
        fn_w_ex3 = repmat(reshape(fn_w,1,1,[]), n, K, 1);
        Pwpp_sum  = squeeze(sum(Pwppre  .* fn_w_ex3, 3)); % 24×K
    else
        F = [F, Pwppre == 0];
    end
    if has_p
        fn_p_ex3 = repmat(reshape(fn_p,1,1,[]), n, K, 1);
        Ppv_sum   = squeeze(sum(Ppvpre  .* fn_p_ex3, 3)); % 24×K
    else
        F = [F, Ppvpre == 0];
    end
    if has_tg
        fn_tg_ex3 = repmat(reshape(fn_tg,1,1,[]), n, K, 1);
        Ptg_sum   = squeeze(sum(Ptgre   .* fn_tg_ex3, 3)); % 24×K
    else
        F = [F, Ptgre == 0, stgre == 0];
    end
    if has_ess
        fn_ess_ex3 = repmat(reshape(fn_ess,1,1,[]), n, K, 1);
        Pess_sum  = squeeze(sum((Pdere - Pcere) .* fn_ess_ex3, 3)); % 24×K
    else
        F = [F, Pcere == 0, Pdere == 0, Re_cess == 0, Re_dess == 0, socre == 0];
    end
    if has_dr
        fn_dr_ex3 = repmat(reshape(fn_dr,1,1,[]), n, K, 1);
        Pdr_sum   = squeeze(sum(Cutre_demand .* fn_dr_ex3, 3)); % 24×K
    else
        F = [F, Cutre_demand == 0];
    end

    F = [F, 0 <= Preal <= Pwpp_sum + Ppv_sum + Ptg_sum + Pess_sum + Pdr_sum];
    %
    F = [F Preal<=1.3*repmat(Pda,1,K)];
    F = [F Preal-0.8*repmat(Pda,1,K)>=-M*(ones(24,K)-mu_1)];
    F = [F Preal<=M*mu_1];
    F = [F Pre==Preal-repmat(Pda,1,K)];

    %Real-Time Generation Constraints for WPP
    if has_w
        Wind_f = Parameter.Wind_f;
        if ndims(Wind_f) == 2
            Wind_f = repmat(reshape(Wind_f,n,1,w_number),1,K,1);
        end
        F = [F, Pwppre == Wind_f .* fn_w_ex3];
    end
    % Real-Time Generation for PV
    if has_p
        Solar_f = Parameter.Solar_f;
        if ndims(Solar_f) == 2
            Solar_f = repmat(reshape(Solar_f,n,1,p_number),1,K,1);
        end
        F = [F, Ppvpre == Solar_f .* fn_p_ex3];
    end
    %-----------------------------------------------------------------------------------------
    
    %tg real-time constraints (matrix form)
    if has_tg
        Pgmin = reshape(Parameter.Pgmin,1,1,[]);%size:1*1*tg_number
        Pgmax = reshape(Parameter.Pgmax,1,1,[]);%size:1*1*tg_number
        Dg = reshape(Parameter.Dg,1,1,[]);%size:1*1*tg_number
        Ug = reshape(Parameter.Ug,1,1,[]);%size:1*1*tg_number
        Pgmin_ex3 = repmat(Pgmin, n, K, 1);%size:24*K*tg_number
        Pgmax_ex3 = repmat(Pgmax, n, K, 1);%size:24*K*tg_number
        F = [F, Pgmin_ex3 .* fn_tg_ex3 .* stgre <= Ptgre <= Pgmax_ex3 .* fn_tg_ex3 .* stgre];
        F = [F, stgre(1,:,:) == 0];
        F = [F, -repmat(Dg,n-1,K,1) <= diff(Ptgre) <= repmat(Ug,n-1,K,1)];
    end
    %-----------------------------------------------------------------------------------------

    %ESS real-time constraints (matrix form)
    if has_ess
        p_cmax = reshape(Parameter.p_cmax,1,1,[]); %size:1*1*ess_number
        p_dmax = reshape(Parameter.p_dmax,1,1,[]); %size:1*1*ess_number
        socmin = reshape(Parameter.socmin,1,1,[]); %socmin:1*1*ess_number
        socmax = reshape(Parameter.socmax,1,1,[]); %socmax:1*1*ess_number
        socstart = reshape(Parameter.socstart,1,1,[]); %socstart:1*1*ess_number
        yt_sin = reshape(Parameter.yt_sin,1,1,[]); %size:1*1*ess_number
        yt_sout = reshape(Parameter.yt_sout,1,1,[]); %size:1*1*ess_number
        Esmax = reshape(Parameter.Esmax,1,1,[]); %size:1*1*ess_number
        p_cmax_ex3 = repmat(p_cmax, n, K, 1); %size:24*K*ess_number
        p_dmax_ex3 = repmat(p_dmax, n, K, 1);  %size:24*K*ess_number
        F = [F, 0 <= Pcere <= p_cmax_ex3 .* fn_ess_ex3 .* Re_cess]; 
        F = [F, 0 <= Pdere <= p_dmax_ex3 .* fn_ess_ex3 .* Re_dess];
        F = [F, 0 <= Re_cess + Re_dess <= 1];
        socmin_ex3 = repmat(socmin, n, K, 1);
        socmax_ex3 = repmat(socmax, n, K, 1);
        socstart_ex3 = repmat(socstart, 1, K, 1);
        yt_sin_ex3 = repmat(yt_sin, n, K, 1);
        yt_sout_ex3 = repmat(yt_sout, n, K, 1);
        Esmax_ex3 = repmat(Esmax, n, K, 1);
        sin_over_es3 = yt_sin_ex3 ./ Esmax_ex3;
        inv_sout_es3 = 1 ./ (yt_sout_ex3 .* Esmax_ex3);
        F = [F, socmin_ex3 <= socre <= socmax_ex3];
        F = [F, socre(24,:,:) == socstart_ex3];
        F = [F, socre(1,:,:) == socstart_ex3 + (Pcere(1,:,:).*sin_over_es3(1,:,:) - Pdere(1,:,:).*inv_sout_es3(1,:,:))];
        F = [F, diff(socre) == (Pcere(2:end,:,:).*sin_over_es3(2:end,:,:) - Pdere(2:end,:,:).*inv_sout_es3(2:end,:,:))];
    end
    %-----------------------------------------------------------------------------------------
    %DR real-time constraints (matrix form)
    if has_dr
        Load_f = Parameter.Load_f;
        if ndims(Load_f) == 2
            Load_f = repmat(reshape(Load_f,n,1,dr_number),1,K,1);
        end
        F = [F, 0 <= Cutre_demand <= 0.85 * Load_f .* fn_dr_ex3];
    end

    %Real_time cost (matrix form)
    ag = reshape(Parameter.ag,1,1,[]);
    bg = reshape(Parameter.bg,1,1,[]);
    cg = reshape(Parameter.cg,1,1,[]);
    as = reshape(Parameter.as,1,1,[]);
    cutting = reshape(Parameter.Cutting,1,1,[]);

    ag_ex3 = repmat(ag, n, K, 1);
    bg_ex3 = repmat(bg, n, K, 1);
    cg_ex3 = repmat(cg, n, K, 1);
    as_ex3 = repmat(as, n, K, 1);
    cutting_ex3 = repmat(cutting, n, K, 1);
    Ctgre = ag_ex3 .* Ptgre + cg_ex3 .* Ptgre.*Ptgre + bg_ex3 .* stgre;
    Cessre = as_ex3 .* (Pcere + Pdere);
    Cdere = cutting_ex3 .* Cutre_demand;

    Cre = sum(Ctgre,3) + sum(Cessre,3) + sum(Cdere,3);
    Cre = squeeze(Cre);

    fre=0;
    %%%real-time profit calculation
    
    Profit = sum(lmpdas .* repmat(Pda,1,K) + lmprts .* Pre - Cre - Parameter.pane * abs(Pre), 1);
    
    %Profit
    f=sum(ps'.*Profit);

    opf=f;
    ops=sdpsettings('solver','gurobi');
    output=optimize(F,-opf,ops);
    
    % Data recording
    Result_dabe=struct;
    Result_dabe.Pda=double(Pda);
    Result_dabe.Pwppda=double(Pwppda);
    Result_dabe.Ppvpda=double(Ppvpda);
    Result_dabe.Ptgda=double(Ptgda);
    Result_dabe.stg=double(stg);
    Result_dabe.Pceda=double(Pceda);
    Result_dabe.Pdeda=double(Pdeda);
    Result_dabe.soc=double(soc);
    Result_dabe.cess=double(cess);
    Result_dabe.dess=double(dess);
    Result_dabe.Cut_demand=double(Cut_demand);

    Result_dabe.mu_1=double(mu_1);
    Result_dabe.Pre=double(Pre);
    Result_dabe.Preal=double(Preal);
    Result_dabe.Pwppre=double(Pwppre);
    Result_dabe.Ppvpre=double(Ppvpre);
    Result_dabe.Ptgre=double(Ptgre);
    Result_dabe.stgre=double(stgre);
    Result_dabe.Pcere=double(Pcere);
    Result_dabe.Pdere=double(Pdere);
    Result_dabe.socre=double(socre);
    Result_dabe.Re_cess=double(Re_cess);
    Result_dabe.Re_dess=double(Re_dess);
    Result_dabe.Cutre_demand=double(Cutre_demand);

    Result_dabe.f=double(opf);
    Result_dabe.Profitda=double(Profitda);
    Result_dabe.Cda=double(Cda);
    Result_dabe.Ctgda=double(Ctgda);
    Result_dabe.Cessda=double(Cessda);
    Result_dabe.Cdrda=double(Cdrda);
    Result_dabe.Cre=double(Cre);
    yalmip('clear');
end




