function Result_rt=Large_Rt_Internal(fn,Parameter,scen,Pda) 
    %% definition of members 
    ps=Parameter.P;
    Lmpda=Parameter.Lmpda(:,scen);
    Lmprt=Parameter.Lmprt(:,scen);
    w_number=Parameter.w_number;%number of WPPs 
    p_number=Parameter.p_number;%number of PVs
    tg_number=Parameter.tg_number;%number of TGs
    ess_number=Parameter.ess_number;%number of ESSs
    dr_number=Parameter.dr_number;%number of DRs

    pwpprt=sdpvar(24,w_number);%wpp output
    ppvprt=sdpvar(24,p_number);%pv output

    ptgprt=sdpvar(24,tg_number);%tg output
    stgrt=binvar(24,tg_number);

    pdert=sdpvar(24,ess_number);%ess discharge power
    pcert=sdpvar(24,ess_number);%ess charge power  
    socrt=sdpvar(24,ess_number);%ess soc
    cessrt=binvar(24,ess_number);%ess charge state
    dessrt=binvar(24,ess_number);%ess discharge state  

    cutre=sdpvar(24,dr_number);%dr load reduction

    ctgre=sdpvar(24,tg_number);%tg cost
    cessre=sdpvar(24,ess_number);%ess cost
    cdrre=sdpvar(24,dr_number);%dr cost
    crt=sdpvar(24,1);%real-time cost
    %% fn index
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

    %% real-time generation
    Prt=sdpvar(24,1);
    Preal=sdpvar(24,1);
    n = size(pwpprt,1);
    fn_w_ex = repmat(reshape(fn_w,1,[]), n, 1);
    fn_p_ex = repmat(reshape(fn_p,1,[]), n, 1);
    fn_tg_ex = repmat(reshape(fn_tg,1,[]), n, 1);
    fn_ess_ex = repmat(reshape(fn_ess,1,[]), n, 1);
    fn_dr_ex = repmat(reshape(fn_dr,1,[]), n, 1);
    %% constraints-deviation and constraints
    F=[];
    F=[F Prt==Preal-Pda];
    %% 实时出力约束
    F=[F 0<=Preal<=pwpprt*reshape(fn_w,[],1)+ppvprt*reshape(fn_p,[],1)+ptgprt*reshape(fn_tg,[],1)+(pdert-pcert)*reshape(fn_ess,[],1)+cutre*reshape(fn_dr,[],1)];%实时总出力约束


    % real-time renewable generation constraints (matrix form)
    if has_w
        Wind_f = Parameter.Wind_f;
        if ndims(Wind_f) == 3
            Wind_r = squeeze(Wind_f(:,scen,:));
        else
            Wind_r = Wind_f;
        end
        F = [F, 0 <= pwpprt <= Wind_r .* fn_w_ex];
    else
        F = [F, pwpprt == 0];
    end
    if has_p
        Solar_f = Parameter.Solar_f;
        if ndims(Solar_f) == 3
            Solar_r = squeeze(Solar_f(:,scen,:));
        else
            Solar_r = Solar_f;
        end
        F = [F, 0 <= ppvprt <= Solar_r .* fn_p_ex];
    else
        F = [F, ppvprt == 0];
    end

    % TG constraints (matrix form)
    if has_tg
        Pgmin = reshape(Parameter.Pgmin,1,[]);
        Pgmax = reshape(Parameter.Pgmax,1,[]);
        Dg = reshape(Parameter.Dg,1,[]);
        Ug = reshape(Parameter.Ug,1,[]);
        Pgmin_ex = repmat(Pgmin, n, 1);
        Pgmax_ex = repmat(Pgmax, n, 1);
        F = [F, Pgmin_ex .* fn_tg_ex .* stgrt <= ptgprt <= Pgmax_ex .* fn_tg_ex .* stgrt];
        F = [F, stgrt(1,:) == 0];
        F = [F, -repmat(Dg,n-1,1) <= diff(ptgprt) <= repmat(Ug,n-1,1)];
    else
        F = [F, ptgprt == 0, stgrt == 0];
    end

    % ESS constraints (matrix form)
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
        F = [F, 0 <= pcert <= p_cmax_ex .* fn_ess_ex .* cessrt];
        F = [F, 0 <= pdert <= p_dmax_ex .* fn_ess_ex .* dessrt];
        F = [F, 0 <= cessrt + dessrt <= 1];
        socmin_ex = repmat(socmin, n, 1);
        socmax_ex = repmat(socmax, n, 1);
        yt_sin_ex = repmat(yt_sin, n, 1);
        yt_sout_ex = repmat(yt_sout, n, 1);
        Esmax_ex = repmat(Esmax, n, 1);
        sin_over_es = yt_sin_ex ./ Esmax_ex;
        inv_sout_es = 1 ./ (yt_sout_ex .* Esmax_ex);
        F = [F, socmin_ex <= socrt <= socmax_ex];
        F = [F, socrt(1,:) == 0.5 + (pcert(1,:).*sin_over_es(1,:) - pdert(1,:).*inv_sout_es(1,:))];
        F = [F, socrt(24,:) == 0.5];
        F = [F, diff(socrt) == (pcert(2:end,:).*sin_over_es(2:end,:) - pdert(2:end,:).*inv_sout_es(2:end,:))];
    else
        F = [F, pcert == 0, pdert == 0, cessrt == 0, dessrt == 0, socrt == 0];
    end

    % DR constraints (matrix form)
    if has_dr
        Load_f = Parameter.Load_f;
        if ndims(Load_f) == 3
            Load_r = squeeze(Load_f(:,scen,:));
        else
            Load_r = Load_f;
        end
        F = [F, 0 <= cutre <= 0.85 * Load_r .* fn_dr_ex];
    else
        F = [F, cutre == 0];
    end
       
   
    ag = reshape(Parameter.ag,1,[]);
    bg = reshape(Parameter.bg,1,[]);
    cg = reshape(Parameter.cg,1,[]);
    ag_ex = repmat(ag, n, 1);
    bg_ex = repmat(bg, n, 1);
    cg_ex = repmat(cg, n, 1);
    Ctgrt = ag_ex .* ptgprt + cg_ex .* (ptgprt .* ptgprt) + bg_ex .* stgrt;
    Cessrt = repmat(reshape(Parameter.as,1,[]), n, 1) .* (pcert + pdert);
    Cdert = repmat(reshape(Parameter.Cutting,1,[]), n, 1) .* cutre;
    Crt = sum(Ctgrt,2) + sum(Cessrt,2) + sum(Cdert,2);

    f=sum(Lmpda.*Pda+Lmprt.*Prt-Parameter.pane*(abs(Prt))-Crt);
    opf=f;
    ops=sdpsettings('solver','gurobi');
    output=optimize(F,-opf,ops);
    rtfresult=struct;
    rtfresult.f=double(opf);
    rtfresult.Prt=double(Prt);
    rtfresult.Preal=double(Preal);
    rtfresult.Pwpprt=double(pwpprt);
    rtfresult.Ppvprt=double(ppvprt);
    rtfresult.Ptgrt=double(ptgprt);
    rtfresult.stgrt=double(stgrt);
    rtfresult.Pdert=double(pdert);
    rtfresult.Pcert=double(pcert);
    rtfresult.socrt=double(socrt);
    rtfresult.cessrt=double(cessrt);
    rtfresult.dessrt=double(dessrt);
    rtfresult.Cut_demand=double(cutre);
    rtfresult.Ctgrt=double(Ctgrt);
    rtfresult.Cessrt=double(Cessrt);
    rtfresult.Cdert=double(Cdert);
    rtfresult.Crt=double(Crt);
    Result_rt=rtfresult;
    yalmip('clear');
end
