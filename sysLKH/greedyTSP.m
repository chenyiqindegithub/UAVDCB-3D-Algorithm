function [best_path, best_length] = greedyTSP(cost_mat)
    n = size(cost_mat, 1);
    M = inf;

    temp_cost_mat = cost_mat;
    temp_cost_mat(1:n+1:end) = M;

    best_path = 1;

    while length(best_path) < n
        current_city = best_path(end);
        costs_from_current = temp_cost_mat(current_city, :);
        costs_from_current(best_path) = M;
        [~, next_city] = min(costs_from_current);
        best_path = [best_path, next_city];
    end

    path_segments_cost = sum(cost_mat(sub2ind(size(cost_mat), best_path(1:end-1), best_path(2:end))));
    closing_cost = cost_mat(best_path(end), best_path(1));
    best_length = path_segments_cost + closing_cost;
end
