function shapley_values = compute_shapley(matrix)
% The reference is Risk-based profit allocation to DERs integrated with a virtual power plant using cooperative Game theory
    [num_rows, num_cols] = size(matrix);
    n = num_cols - 1;
    shapley_values = zeros(1, n);
    
    for player = 1:n
        temp_value = 0;
        for k = 0:2^n - 1
            coalition = find(dec2bin(k, n) == '1');
            
            if ~ismember(player, coalition)
                binary_coalition = zeros(1, n);
                binary_coalition(coalition) = 1;
                
                matching_rows = ismember(matrix(:, 1:n), binary_coalition, 'rows');
                prev_value = 0;
                if any(matching_rows)
                    prev_value = matrix(matching_rows, end);
                end
                
                coalition_with_player = [coalition, player];
                binary_coalition_with_player = zeros(1, n);
                binary_coalition_with_player(coalition_with_player) = 1;
                
                matching_rows = ismember(matrix(:, 1:n), binary_coalition_with_player, 'rows');
                new_value = 0;
                if any(matching_rows)
                    new_value = matrix(matching_rows, end);
                end
                
                marginal_contribution = new_value - prev_value;
                temp_value = temp_value + (factorial(length(coalition)) * factorial(n - length(coalition) - 1) / factorial(n)) * marginal_contribution;
            end
        end
        shapley_values(player) = temp_value;
    end
end
