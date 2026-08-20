function [nucleolus, sigma_trace] = compute_Nucleolus(matrix)
    % 基础参数
    % It is based on the paper Risk-based profit allocation to DERs integrated with a virtual power plant using cooperative Game theory    [num_rows, num_cols] = size(matrix);
    n = num_cols - 1; % 玩家数量
    fg = matrix(2:end, 1:n); % 联盟的二进制表示（每一行为一个联盟）
    v = matrix(2:end, end);  % 每个联盟的收益
    % 初始化变量
    x = sdpvar(n, 1);     % 玩家分配变量
    sigma = sdpvar(1, 1); % 最大超额
    tol = 1e-6;
    ops = sdpsettings('solver','gurobi','verbose',0);
    % 初始化约束：效率 + non-negativity + 最大超额约束
    F_base = [sum(x) == v(end)];
    F = [F_base, fg(1:end-1,:) * x + sigma >= v(1:end-1)];
    %Some nucleolus methods require the allocation needs to be the
    %imputation
    % for i=1:n
    %     F=[F x(i)>=v(2^(n-i))];
    % end
    sigma_trace = [];
    tight_set = [];
  
    optimize(F, sigma, ops);
    sigma_val = value(sigma);
    x_val = value(x);
    sigma_trace(end+1) = sigma_val;
    excess = v(1:end-1) - fg(1:end-1,:) * x_val;
    tight_idx = find(abs(excess - sigma_val) < tol);
    while true
        % 添加当前 tight constraints 为等式约束
        eqs = [];
        for i = tight_idx(:)'
            eqs = [eqs, fg(i,:) * x == v(i) - sigma_val];
        end
        % 构建新模型
        sigma = sdpvar(1, 1); % 重定义 sigma
        F_new = [F_base, eqs, fg(1:end-1,:) * x + sigma >= v(1:end-1)];
        optimize(F_new, sigma, ops);
        sigma_val_new = value(sigma);
        x_val_new = value(x);
        sigma_trace(end+1) = sigma_val_new;
    
        % 找出新的 tight constraints
        excess = v(1:end-1) - fg(1:end-1,:) * x_val_new;
        tight_idx_new = find(abs(excess - sigma_val_new) < tol);
    
        % 判断是否终止
        if abs(sigma_val - sigma_val_new) < tol && ...
           length(tight_idx_new) == length(tight_idx) && ...
           all(ismember(tight_idx_new, tight_idx)) && all(ismember(tight_idx, tight_idx_new))
            break;
        end
        % 更新
        tight_idx = unique([tight_idx; tight_idx_new]);
        sigma_val = sigma_val_new;
        x_val = x_val_new;
    end

    nucleolus = x_val_new;
end
