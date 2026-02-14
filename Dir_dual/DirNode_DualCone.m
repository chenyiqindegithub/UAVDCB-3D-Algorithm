function [schemes_nodes_cell, schemes_coverage_cell, run_time] = DirNode_DualCone(points_global, D1, D2, pos_0, cone_angle_phi1_deg, cone_angle_phi2_deg, angle_between_axes_varphi_deg, cover_mat_global, sid)


    % 初始化
    run_time = 0; tic;
    schemes_nodes_cell = {}; 
    schemes_coverage_cell = {};
    

    cover_node_global_indices = find(cover_mat_global(:, sid));
    if isempty(cover_node_global_indices), run_time = toc; return; end
    
    points_local_coords = points_global(cover_node_global_indices, :);
    vectors_0_to_points = points_local_coords - pos_0;
    num_covered_nodes = size(points_local_coords, 1);
    dims = size(points_global, 2);

    is_node_at_pos0 = false;
    spot_node_local_idx_list = find(vecnorm(vectors_0_to_points, 2, 2) <= 1.0e-6);
    if ~isempty(spot_node_local_idx_list)
        is_node_at_pos0 = true;
    end

    % 为每个节点生成候选姿态 
    candidate_orientations = [];
    
    for i = 1:num_covered_nodes
        v_i = vectors_0_to_points(i, :);
        if norm(v_i) < 1e-6, continue; end
        
        dir1_A = v_i;
        dir2_A = get_rotated_dir(dir1_A, angle_between_axes_varphi_deg, dims);
        candidate_orientations = [candidate_orientations; dir1_A, dir2_A];
        
        dir2_B = v_i;
        dir1_B = get_rotated_dir(dir2_B, -angle_between_axes_varphi_deg, dims);
        candidate_orientations = [candidate_orientations; dir1_B, dir2_B];
    end
    

    all_nodes_temp = {};
    all_coverage_temp = {};

    for i = 1:size(candidate_orientations, 1)
        orientation = candidate_orientations(i, :);
        dir1 = orientation(1:dims); dir2 = orientation(dims+1:end);
        
        set1 = is_point_in_cone_local(points_local_coords, pos_0, dir1, D1, cone_angle_phi1_deg);
        set2 = is_point_in_cone_local(points_local_coords, pos_0, dir2, D2, cone_angle_phi2_deg);
        
        total_set = union(set1, set2);
        if is_node_at_pos0, total_set = union(total_set, spot_node_local_idx_list); end
        
        if isempty(total_set), continue; end

        coverage = zeros(1, length(total_set));
        is_in1 = ismember(total_set, set1); is_in2 = ismember(total_set, set2);
        coverage(is_in1 & ~is_in2) = 1; coverage(~is_in1 & is_in2) = 2; coverage(is_in1 & is_in2) = 3;
        coverage(ismember(total_set, spot_node_local_idx_list)) = 3;
        
        all_nodes_temp{end+1} = cover_node_global_indices(total_set)';
        all_coverage_temp{end+1} = coverage;
    end
    

    if isempty(all_nodes_temp) && is_node_at_pos0
        final_nodes_global = cover_node_global_indices(spot_node_local_idx_list)';
        final_coverage_status = 3 * ones(1, length(final_nodes_global));
        all_nodes_temp{1} = final_nodes_global;
        all_coverage_temp{1} = final_coverage_status;
    end
    
    if isempty(all_nodes_temp), run_time = toc; return; end
    
    schemes_nodes_cell = all_nodes_temp;
    schemes_coverage_cell = all_coverage_temp;

    run_time = toc;
end

% 辅助函数

function rotated_vec = get_rotated_dir(vec_to_rotate, angle_deg, dims)
    if norm(vec_to_rotate) < 1e-9, rotated_vec = vec_to_rotate; return; end
    if dims == 2
        theta_rad = deg2rad(angle_deg);
        R = [cos(theta_rad), -sin(theta_rad); sin(theta_rad), cos(theta_rad)];
        rotated_vec = (vec_to_rotate(1:2) * R);
    else % 3D
        up_vector = [0 0 1]; 
        axis_of_rotation = cross(vec_to_rotate, up_vector);
        if norm(axis_of_rotation) < 1e-6
            axis_of_rotation = cross(vec_to_rotate, [0 1 0]);
        end
        rotated_vec = rodrigues_rot(vec_to_rotate, axis_of_rotation, angle_deg);
    end
end

function v_rot = rodrigues_rot(v, k, theta_deg)
    if norm(k) < 1e-9, v_rot = v; return; end
    theta_rad = deg2rad(theta_deg);
    k = k/norm(k);
    v_rot = v * cos(theta_rad) + cross(k, v) * sin(theta_rad) + k * dot(k, v) * (1 - cos(theta_rad));
end

function set_indices = is_point_in_cone_local(points, charge_pos, charge_direction, line_len, angle_deg)
    set_indices = [];
    if norm(charge_direction) < 1e-9, return; end
    charge_direction_norm = charge_direction / norm(charge_direction);
    cos_half_angle = cos(deg2rad(angle_deg) / 2);
    for i = 1:size(points, 1)
        vector = points(i, :) - charge_pos;
        dist = norm(vector);
        if dist > line_len + 1e-6, continue; end
        if dist < 1e-6, set_indices = [set_indices, i]; continue; end
        cos_angle_between = dot(vector, charge_direction_norm) / dist;
        if cos_angle_between >= cos_half_angle - 1e-6
            set_indices = [set_indices, i];
        end
    end
end