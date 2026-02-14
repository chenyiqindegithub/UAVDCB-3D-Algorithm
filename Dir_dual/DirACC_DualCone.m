function [schemes_nodes_cell, schemes_coverage_cell, run_time] = DirACC_DualCone(points_global, D1, D2, pos_0, cone_angle_phi1_deg, cone_angle_phi2_deg, angle_between_axes_varphi_deg, cover_mat_global, sid)

    run_time = 0; tic;
    schemes_nodes_cell = {}; schemes_coverage_cell = {};
    
    cover_node_global_indices = find(cover_mat_global(:, sid));
    if isempty(cover_node_global_indices), run_time = toc; return; end
    
    points_covered_coords = points_global(cover_node_global_indices, :);
    vectors_0_to_points = points_covered_coords - pos_0;
    
    is_at_origin_mask = vecnorm(vectors_0_to_points, 2, 2) < 1e-4;
    origin_nodes_global_indices = cover_node_global_indices(is_at_origin_mask);
    
    non_origin_nodes_mask = ~is_at_origin_mask;
    cover_node_global_indices_for_algo = cover_node_global_indices(non_origin_nodes_mask);
    points_covered_coords_for_algo = points_covered_coords(non_origin_nodes_mask, :);
    vectors_0_to_points_for_algo = vectors_0_to_points(non_origin_nodes_mask, :);
    num_covered_nodes_for_algo = size(points_covered_coords_for_algo, 1);
    
    if isempty(cover_node_global_indices_for_algo)
        if ~isempty(origin_nodes_global_indices)
            schemes_nodes_cell{1,1} = sort(origin_nodes_global_indices)';
            schemes_coverage_cell{1,1} = 3 * ones(1, length(origin_nodes_global_indices));
            run_time = toc; return;
        else
            run_time = toc; return;
        end
    end

    cone_half_angle_phi1_rad = deg2rad(cone_angle_phi1_deg) / 2;
    cone_half_angle_phi2_rad = deg2rad(cone_angle_phi2_deg) / 2;
    varphi_rad = deg2rad(angle_between_axes_varphi_deg);
    
    cone1_angle_threshold_deg = cone_angle_phi1_deg;
    cone2_angle_threshold_deg = cone_angle_phi2_deg;
    up_vector_for_e2_calc_default = [0, 0, 1];

    candidate_e1_directions = [];
    
    for i = 1:num_covered_nodes_for_algo
        vector_0i = vectors_0_to_points_for_algo(i, :);
        candidate_e1_directions = [candidate_e1_directions; normalize_vec(vector_0i)];

        for j = (i+1):num_covered_nodes_for_algo
            vector_0j = vectors_0_to_points_for_algo(j, :);
            angle_ij_deg = acosd(max(min(dot(vector_0i, vector_0j) / (norm(vector_0i) * norm(vector_0j)), 1), -1));
            
            Del_bisector_ij = normalize_vec(vector_0i + vector_0j);
            if norm(Del_bisector_ij) < 1e-6, continue; end

            if angle_ij_deg <= cone1_angle_threshold_deg + 1e-4
                candidate_e1_directions = [candidate_e1_directions; Del_bisector_ij];
            end
            
            if angle_ij_deg <= cone2_angle_threshold_deg + 1e-4
                up_vector_for_e1_calc_ij = determine_up_vector(Del_bisector_ij, up_vector_for_e2_calc_default);
                e1_from_bisector_as_e2 = calculate_e1_from_e2(Del_bisector_ij, varphi_rad, up_vector_for_e1_calc_ij);
                candidate_e1_directions = [candidate_e1_directions; e1_from_bisector_as_e2];
            end
        end
    end
    
    if ~isempty(candidate_e1_directions)
        candidate_e1_directions = unique(candidate_e1_directions, 'rows');
    end

    generated_schemes_data = struct('nodes',{}, 'e1_dir',{}, 'e2_dir',{});
    if ~isempty(candidate_e1_directions)
        for k_dir = 1:size(candidate_e1_directions, 1)
            current_e1_dir = candidate_e1_directions(k_dir, :);
            up_vector_for_e2_final = determine_up_vector(current_e1_dir, up_vector_for_e2_calc_default);
            current_e2_dir = calculate_e2_from_e1(current_e1_dir, varphi_rad, up_vector_for_e2_final);
            final_covered_local_indices = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0, current_e1_dir, current_e2_dir, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
            
            if ~isempty(final_covered_local_indices)
                current_set_global_indices = cover_node_global_indices_for_algo(final_covered_local_indices);
                idx_new_scheme = length(generated_schemes_data) + 1;
                generated_schemes_data(idx_new_scheme).nodes = sort(current_set_global_indices(:)');
                generated_schemes_data(idx_new_scheme).e1_dir = current_e1_dir;
                generated_schemes_data(idx_new_scheme).e2_dir = current_e2_dir;
            end
        end
    end
    
    final_schemes_data = generated_schemes_data;

    all_currently_covered_nodes = [];
    if ~isempty(final_schemes_data), all_currently_covered_nodes = unique([final_schemes_data.nodes]); end
    uncovered_nodes_global_indices = setdiff(cover_node_global_indices_for_algo', all_currently_covered_nodes);
    
    if ~isempty(uncovered_nodes_global_indices)
        newly_generated_schemes = struct('nodes',{}, 'e1_dir',{}, 'e2_dir',{});
        for i = 1:length(uncovered_nodes_global_indices)
            node_idx_global = uncovered_nodes_global_indices(i);
            vector_to_uncovered_node = points_global(node_idx_global, :) - pos_0;
            dir_to_node = normalize_vec(vector_to_uncovered_node);
            up_vec_hint = determine_up_vector(dir_to_node, up_vector_for_e2_calc_default);
            node_local_idx_for_algo = find(cover_node_global_indices_for_algo == node_idx_global, 1);
            
            e1_candidate_A = dir_to_node;
            e2_candidate_A = calculate_e2_from_e1(e1_candidate_A, varphi_rad, up_vec_hint);
            covered_indices_A = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0, e1_candidate_A, e2_candidate_A, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
            
            if ~isempty(node_local_idx_for_algo) && ismember(node_local_idx_for_algo, covered_indices_A)
                nodes_A_global = sort(cover_node_global_indices_for_algo(covered_indices_A)');
                new_scheme_idx = length(newly_generated_schemes) + 1;
                newly_generated_schemes(new_scheme_idx).nodes = nodes_A_global;
                newly_generated_schemes(new_scheme_idx).e1_dir = e1_candidate_A;
                newly_generated_schemes(new_scheme_idx).e2_dir = e2_candidate_A;
                continue;
            end
    
            e2_target_dir_B = dir_to_node;
            e1_candidate_B = calculate_e1_from_e2(e2_target_dir_B, varphi_rad, up_vec_hint);
            covered_indices_B = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0, e1_candidate_B, e2_target_dir_B, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
            
            if ~isempty(covered_indices_B)
                nodes_B_global = sort(cover_node_global_indices_for_algo(covered_indices_B)');
                new_scheme_idx = length(newly_generated_schemes) + 1;
                newly_generated_schemes(new_scheme_idx).nodes = nodes_B_global;
                newly_generated_schemes(new_scheme_idx).e1_dir = e1_candidate_B;
                newly_generated_schemes(new_scheme_idx).e2_dir = e2_target_dir_B;
            end
        end
        combined_schemes = [final_schemes_data, newly_generated_schemes];
        final_schemes_data = combined_schemes;
    end

    if ~isempty(origin_nodes_global_indices)
        if isempty(final_schemes_data)
            final_schemes_data(1).nodes = sort(origin_nodes_global_indices(:)');
            final_schemes_data(1).e1_dir = [0 0 1];
            final_schemes_data(1).e2_dir = calculate_e2_from_e1([0 0 1], varphi_rad, [0 1 0]);
        else
            for k = 1:length(final_schemes_data)
                final_schemes_data(k).nodes = union(final_schemes_data(k).nodes, origin_nodes_global_indices(:)');
            end
        end
    end

    if ~isempty(final_schemes_data)
        for k_scheme = 1:length(final_schemes_data)
            current_scheme = final_schemes_data(k_scheme);
            current_e1_dir = current_scheme.e1_dir;
            current_e2_dir = current_scheme.e2_dir;
            nodes_in_scheme_global = current_scheme.nodes;
            
            nodes_in_cone1_local = is_point_in_single_cone(points_covered_coords, pos_0, current_e1_dir, D1, cone_half_angle_phi1_rad);
            nodes_in_cone2_local = is_point_in_single_cone(points_covered_coords, pos_0, current_e2_dir, D2, cone_half_angle_phi2_rad);
            
            coverage_status_for_this_scheme = zeros(1, length(nodes_in_scheme_global));
            
            for k_node = 1:length(nodes_in_scheme_global)
                node_global_idx = nodes_in_scheme_global(k_node);
                node_local_idx = find(cover_node_global_indices == node_global_idx, 1);
                
                if ~isempty(node_local_idx)
                    is_in_c1 = ismember(node_local_idx, nodes_in_cone1_local);
                    is_in_c2 = ismember(node_local_idx, nodes_in_cone2_local);
                    
                    if is_in_c1 && is_in_c2
                        coverage_status_for_this_scheme(k_node) = 3;
                    elseif is_in_c1
                        coverage_status_for_this_scheme(k_node) = 1;
                    elseif is_in_c2
                        coverage_status_for_this_scheme(k_node) = 2;
                    end
                end
            end
            schemes_nodes_cell{k_scheme, 1} = nodes_in_scheme_global;
            schemes_coverage_cell{k_scheme, 1} = coverage_status_for_this_scheme;
        end
    end
    run_time = toc;
end

function normalized_v = normalize_vec(v)
    norm_v = norm(v);
    if norm_v < 1e-9, normalized_v = zeros(size(v)); else, normalized_v = v / norm_v; end
end

function up_vec = determine_up_vector(e_unit_vec, default_up_vec)
    up_vec = default_up_vec;
    if norm(cross(e_unit_vec, up_vec)) < 1e-6
        if all(abs(e_unit_vec - [0,0,1]) < 1e-6) || all(abs(e_unit_vec - [0,0,-1]) < 1e-6)
            up_vec = [0, 1, 0];
        else
            up_vec = [0, 0, 1];
        end
        if norm(cross(e_unit_vec, up_vec)) < 1e-6
            if all(abs(e_unit_vec - [0,1,0]) < 1e-6) || all(abs(e_unit_vec - [0,-1,0]) < 1e-6)
                 up_vec = [1, 0, 0];
            else
                 up_vec = [0, 1, 0];
                 if norm(cross(e_unit_vec, up_vec)) < 1e-6
                     up_vec = [1,0,0];
                 end
            end
        end
    end
end

function e2_unit_vec = calculate_e2_from_e1(e1_unit_vec, varphi_rad, up_vector_hint)
    e1_unit_vec_n = normalize_vec(e1_unit_vec(:)');
    k = cross(e1_unit_vec_n, up_vector_hint);
    if norm(k) < 1e-60
        if abs(e1_unit_vec_n(1)) < (1-1e-6) && abs(e1_unit_vec_n(2)) < (1-1e-6)
            k_temp_axis = [0,0,1];
        elseif abs(e1_unit_vec_n(1)) < (1-1e-6) && abs(e1_unit_vec_n(3)) < (1-1e-6)
            k_temp_axis = [0,1,0];
        else
            k_temp_axis = [1,0,0];
        end
        k = cross(e1_unit_vec_n, k_temp_axis);
    end
    k = normalize_vec(k);
    v_to_rotate = e1_unit_vec_n;
    theta = varphi_rad;
    e2_unit_vec = v_to_rotate*cos(theta) + cross(k,v_to_rotate)*sin(theta) + k*dot(k,v_to_rotate)*(1-cos(theta));
    e2_unit_vec = normalize_vec(e2_unit_vec);
end

function e1_unit_vec = calculate_e1_from_e2(e2_unit_vec, varphi_rad, up_vector_hint)
    e2_unit_vec_n = normalize_vec(e2_unit_vec(:)');
    k = cross(e2_unit_vec_n, up_vector_hint);
    if norm(k) < 1e-60
        if abs(e2_unit_vec_n(1)) < (1-1e-6) && abs(e2_unit_vec_n(2)) < (1-1e-6)
            k_temp_axis = [0,0,1];
        elseif abs(e2_unit_vec_n(1)) < (1-1e-6) && abs(e2_unit_vec_n(3)) < (1-1e-6)
            k_temp_axis = [0,1,0];
        else
            k_temp_axis = [1,0,0];
        end
        k = cross(e2_unit_vec_n, k_temp_axis);
    end
    k = normalize_vec(k);
    v_to_rotate = e2_unit_vec_n;
    theta = -varphi_rad;
    e1_unit_vec = v_to_rotate*cos(theta) + cross(k,v_to_rotate)*sin(theta) + k*dot(k,v_to_rotate)*(1-cos(theta));
    e1_unit_vec = normalize_vec(e1_unit_vec);
end

function set_indices = is_point_in_single_cone(points_coords, cone_vertex_pos, cone_axis_unit_vector, cone_max_dist_D, cone_half_angle_rad)
    set_indices = [];
    num_points = size(points_coords, 1);
    if num_points == 0, return; end
    
    cone_axis_unit_vector_n = normalize_vec(cone_axis_unit_vector);
    if norm(cone_axis_unit_vector_n) < 1e-9, return; end
    
    for i = 1:num_points
        vector_to_point = points_coords(i, :) - cone_vertex_pos;
        dist_to_vertex = norm(vector_to_point);
        
        if dist_to_vertex > cone_max_dist_D + 1e-5, continue; end
        if dist_to_vertex < 1e-5, set_indices = [set_indices, i]; continue; end
        
        unit_vector_to_point = vector_to_point / dist_to_vertex;
        cos_angle_val = dot(unit_vector_to_point, cone_axis_unit_vector_n);
        
        if cos_angle_val >= cos(cone_half_angle_rad + 1e-5)
            set_indices = [set_indices, i];
        end
    end
end

function covered_indices_local = get_nodes_in_dual_cone_explicit(points_local_coords, pos_0, e1_unit_vec, e2_unit_vec, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad)
    nodes_in_cone1 = is_point_in_single_cone(points_local_coords, pos_0, e1_unit_vec, D1, cone_half_angle_phi1_rad);
    nodes_in_cone2 = is_point_in_single_cone(points_local_coords, pos_0, e2_unit_vec, D2, cone_half_angle_phi2_rad);
    covered_indices_local = union(nodes_in_cone1, nodes_in_cone2);
end
