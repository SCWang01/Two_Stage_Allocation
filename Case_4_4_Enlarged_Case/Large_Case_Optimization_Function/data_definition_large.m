%% 该文件为定义文件，用于定义内部所用的数据集
%definite k to get the real-time response
function Parameter = data_definition_large(k,Num_List)
    Parameter=Creating(k,Num_List);
    Parameter.w_number=Num_List(1);
    Parameter.p_number=Num_List(2);
    Parameter.tg_number=Num_List(3);
    Parameter.ess_number=Num_List(4);
    Parameter.dr_number=Num_List(5);
    
    % Uncertain Generation
    Parameter.Wind_f=Parameter.Wind_f(:,:,1:Parameter.w_number);
    Parameter.Solar_f=Parameter.Solar_f(:,:,1:Parameter.p_number);
    Parameter.Load_f=Parameter.Load_f(:,:,1:Parameter.dr_number);
    Parameter.Load_base=Parameter.Load_base(:,1:Parameter.dr_number);
    %% tg constraints
    Parameter.ag=[100,80,100];%%￥/(MW*15min)
    Parameter.bg=[40,20,80];%%￥/(15min)
    Parameter.cg=[80,70,90];%%￥/(MW*15min)^2
    Parameter.Pgmin=[0.1,0.1,0.1];%
    Parameter.Pgmax=[6,8,5];%出力上限
    Parameter.Ug=[1.5,2.0,2.0];
    Parameter.Dg=[1.5,2.0,2.0];
    Parameter.ag= Parameter.ag(1:Parameter.tg_number);
    Parameter.bg= Parameter.bg(1:Parameter.tg_number);
    Parameter.cg= Parameter.cg(1:Parameter.tg_number);
    Parameter.Pgmin= Parameter.Pgmin(1:Parameter.tg_number);
    Parameter.Pgmax= Parameter.Pgmax(1:Parameter.tg_number);
    Parameter.Ug= Parameter.Ug(1:Parameter.tg_number);
    Parameter.Dg= Parameter.Dg(1:Parameter.tg_number);


    %% 储能的情况 ESS1
    Parameter.as=zeros(Parameter.ess_number,1);%运行成本
    Parameter.yt_sout=[0.95,0.99,0.95,0.95,0.98];
    Parameter.yt_sin=[0.95,0.99,0.95,0.95,0.98];%能量消耗常数
    Parameter.p_cmax=[1.5,1.5,1.8,2.0,2.5];%充电功率
    Parameter.p_dmax=[1.5,1.5,1.8,2.0,2.5];%放电功率
    Parameter.socmin=[0.05,0.05,0.1,0.1,0.08];%最小荷电状态
    Parameter.socmax=[0.95,0.95,0.98,0.98,0.99];%最大荷电状态
    Parameter.Esmax=[8,8,10,10,12];%MW*15min
    Parameter.socstart=[0.5,0.5,0.5,0.5,0.5];
    Parameter.as= Parameter.as(1:Parameter.ess_number);
    Parameter.yt_sout= Parameter.yt_sout(1:Parameter.ess_number);
    Parameter.yt_sin= Parameter.yt_sin(1:Parameter.ess_number);
    Parameter.p_cmax= Parameter.p_cmax(1:Parameter.ess_number);
    Parameter.p_dmax= Parameter.p_dmax(1:Parameter.ess_number);
    Parameter.socmin= Parameter.socmin(1:Parameter.ess_number);
    Parameter.socmax= Parameter.socmax(1:Parameter.ess_number);
    Parameter.Esmax= Parameter.Esmax(1:Parameter.ess_number);
    Parameter.socstart= Parameter.socstart(1:Parameter.ess_number);

    %% cutting
    Parameter.Cutting=[685.00,651.45,687.20,700.56,689.86,633.36,685.50,613.38,689.13,695.31,695.09,...
        691.04,675.74,685.54,662.64];%￥/MW*15min
    Parameter.Cutting=Parameter.Cutting(1:Parameter.dr_number);
    
    %% 惩罚水平
    Parameter.pane=500;%惩罚水平￥/MW*15min
    Parameter.res=0;
end

