function [schemes_nodes_cell, schemes_coverage_cell, run_time] = DirGCC_DualCone_Sampled(points_global, D1, D2, pos_0, cone_angle_phi1_deg, cone_angle_phi2_deg, angle_between_axes_varphi_deg, cover_mat_global, sid)
% DirGCC_DualCone_NoDedup: 通过启发式与采样确定双锥充电方向（无去重）

%% 0. 内部参数与初始化
sample_angle_step_e2_deg = 30;
run_time = 0;
tic;

schemes_nodes_cell = {};
schemes_coverage_cell = {};

cone_half_angle_phi1_rad = deg2rad(cone_angle_phi1_deg) / 2;
cone_half_angle_phi2_rad = deg2rad(cone_angle_phi2_deg) / 2;
varphi_rad = deg2rad(angle_between_axes_varphi_deg);
up_vector_for_e2_calc_default = [0, 0, 1]; 

dim_points = size(points_global, 2);
if dim_points == 2
    points_global_3d = [points_global, zeros(size(points_global, 1), 1)];
    pos_0_3d = [pos_0, 0];
elseif dim_points == 3
    points_global_3d = points_global;
    pos_0_3d = pos_0;
else
    error('输入点的维度必须是2或3。');
end

%% 1. 阶段一：确定相关节点
cover_node_global_indices = find(cover_mat_global(:, sid));
if isempty(cover_node_global_indices)
    run_time = toc;
    return;
end

points_covered_coords = points_global_3d(cover_node_global_indices, :);
vectors_0_to_points = points_covered_coords - pos_0_3d;

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

%% 2. 阶段二：启发式生成候选双锥体朝向
candidate_e1_axes = [];
paired_node_local_ids = [];

for idx_A = 1:num_covered_nodes_for_algo
    for idx_B = idx_A:num_covered_nodes_for_algo
        if idx_A == idx_B, continue; end
        
        vec_OA = vectors_0_to_points_for_algo(idx_A, :);
        vec_OB = vectors_0_to_points_for_algo(idx_B, :);
        
        norm_OA = norm(vec_OA); 
        norm_OB = norm(vec_OB);
        
        cos_angle_AOB = max(min(dot(vec_OA, vec_OB) / (norm_OA * norm_OB), 1), -1);
        angle_AOB_rad = acos(cos_angle_AOB);
        
        if rad2deg(angle_AOB_rad) <= cone_angle_phi1_deg + 1e-4
            node_A_coord = points_covered_coords_for_algo(idx_A, :);
            node_B_coord = points_covered_coords_for_algo(idx_B, :);

            p_dirs_vecs = findCircleCenter(node_A_coord, node_B_coord, pos_0_3d, cone_angle_phi1_deg);
            if ~isempty(p_dirs_vecs), candidate_e1_axes = [candidate_e1_axes; p_dirs_vecs]; end
            
            gcc_dirs_vecs = gcc_direction(node_A_coord, node_B_coord, pos_0_3d, cone_angle_phi1_deg);
            if ~isempty(gcc_dirs_vecs), candidate_e1_axes = [candidate_e1_axes; gcc_dirs_vecs]; end

            paired_node_local_ids = [paired_node_local_ids, idx_A, idx_B];
        end
    end
end

lone_nodes_local_ids = setdiff(1:num_covered_nodes_for_algo, unique(paired_node_local_ids));

for idx_lone = lone_nodes_local_ids
    vec_to_lone = vectors_0_to_points_for_algo(idx_lone, :);
    candidate_e1_axes = [candidate_e1_axes; vec_to_lone];
end

if ~isempty(candidate_e1_axes)
    valid_rows = vecnorm(candidate_e1_axes, 2, 2) > 1e-6;
    candidate_e1_axes = candidate_e1_axes(valid_rows, :);
    if ~isempty(candidate_e1_axes)
        candidate_e1_axes = candidate_e1_axes ./ vecnorm(candidate_e1_axes, 2, 2);
    end
end

candidate_e1_e2_pairs = {};
sample_step_rad = deg2rad(sample_angle_step_e2_deg);

for i = 1:size(candidate_e1_axes, 1)
    e1 = candidate_e1_axes(i, :);
    
    aux_vec = [0,0,1];
    if abs(dot(e1, aux_vec)) > 0.995, aux_vec = [1,0,0]; end
    
    v_perp1 = cross(e1, aux_vec);
    if norm(v_perp1) < 1e-6, continue; end
    v_perp1_unit = v_perp1 / norm(v_perp1);
    v_perp2_unit = cross(e1, v_perp1_unit);
    
    for theta = 0 : sample_step_rad : (2*pi - sample_step_rad/2)
        e2 = e1 * cos(varphi_rad) + (v_perp1_unit * cos(theta) + v_perp2_unit * sin(theta)) * sin(varphi_rad);
        if norm(e2) > 1e-6, candidate_e1_e2_pairs{end+1, 1} = [e1; e2 / norm(e2)]; end
    end
end


generated_schemes_data = struct('nodes',{}, 'e1_dir',{}, 'e2_dir',{});
if ~isempty(candidate_e1_e2_pairs)
    for k = 1:length(candidate_e1_e2_pairs)
        axes_pair = candidate_e1_e2_pairs{k};
        current_e1_dir = axes_pair(1, :);
        current_e2_dir = axes_pair(2, :);
        
        final_covered_local_indices = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0_3d, current_e1_dir, current_e2_dir, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
        
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
        vector_to_uncovered_node = points_global_3d(node_idx_global, :) - pos_0_3d;
        dir_to_node = normalize_vec(vector_to_uncovered_node);
        up_vec_hint = determine_up_vector(dir_to_node, up_vector_for_e2_calc_default);
        node_local_idx_for_algo = find(cover_node_global_indices_for_algo == node_idx_global, 1);
        
       
        e1_candidate_A = dir_to_node;
        e2_candidate_A = calculate_e2_from_e1(e1_candidate_A, varphi_rad, up_vec_hint);
        covered_indices_A = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0_3d, e1_candidate_A, e2_candidate_A, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
        
        if ~isempty(node_local_idx_for_algo) && ismember(node_local_idx_for_algo, covered_indices_A)
            nodes_A_global = sort(cover_node_global_indices_for_algo(covered_indices_A)');
            new_scheme_idx = length(newly_generated_schemes) + 1;
            newly_generated_schemes(new_scheme_idx).nodes = nodes_A_global;
            newly_generated_schemes(new_scheme_idx).e1_dir = e1_candidate_A;
            newly_generated_schemes(new_scheme_idx).e2_dir = e2_candidate_A;
            continue; 
        end

      
        e2_target_dir_B = dir_to_node;
        e1_candidate_B = calculate_e2_from_e1(e2_target_dir_B, varphi_rad, up_vec_hint);
        covered_indices_B = get_nodes_in_dual_cone_explicit(points_covered_coords_for_algo, pos_0_3d, e1_candidate_B, e2_target_dir_B, D1, cone_half_angle_phi1_rad, D2, cone_half_angle_phi2_rad);
        
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

%% 最终覆盖计算
if ~isempty(final_schemes_data)
    for k_scheme = 1:length(final_schemes_data)
        current_scheme = final_schemes_data(k_scheme);
        current_e1_dir = current_scheme.e1_dir;
        current_e2_dir = current_scheme.e2_dir;
        nodes_in_scheme_global = current_scheme.nodes;
        
        nodes_in_cone1_local = is_point_in_single_cone(points_covered_coords, pos_0_3d, current_e1_dir, D1, cone_half_angle_phi1_rad);
        nodes_in_cone2_local = is_point_in_single_cone(points_covered_coords, pos_0_3d, current_e2_dir, D2, cone_half_angle_phi2_rad);
        
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


%% 辅助函数

function normalized_v = normalize_vec(v)
    norm_v = norm(v);
    if norm_v < 1e-9, normalized_v = zeros(size(v)); else, normalized_v = v / norm_v; end
end

function up_vec = determine_up_vector(e1_unit_vec, default_up_vec)
    up_vec = default_up_vec;
    if norm(cross(e1_unit_vec, up_vec)) < 1e-6
        if all(abs(up_vec - [0,0,1]) < 1e-6) || all(abs(up_vec - [0,0,-1]) < 1e-6)
            up_vec = [0, 1, 0];
        else
            up_vec = [0, 0, 1];
        end
        if norm(cross(e1_unit_vec, up_vec)) < 1e-6
            if all(abs(up_vec - [0,1,0]) < 1e-6) || all(abs(up_vec - [0,1,-0]) < 1e-6)
                 up_vec = [1, 0, 0];
            else
                 up_vec = [0, 1, 0];
                 if norm(cross(e1_unit_vec, up_vec)) < 1e-6
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
        if abs(e1_unit_vec_n(1)) < (1-1e-6) && abs(e1_unit_vec_n(2)) < (1-1e-6), k_temp_axis = [0,0,1];
        elseif abs(e1_unit_vec_n(1)) < (1-1e-6) && abs(e1_unit_vec_n(3)) < (1-1e-6), k_temp_axis = [0,1,0];
        else k_temp_axis = [1,0,0]; end
        k = cross(e1_unit_vec_n, k_temp_axis);
    end
    k = normalize_vec(k);
    v_to_rotate = e1_unit_vec_n; theta = varphi_rad;
    e2_unit_vec = v_to_rotate*cos(theta) + cross(k,v_to_rotate)*sin(theta) + k*dot(k,v_to_rotate)*(1-cos(theta));
    e2_unit_vec = normalize_vec(e2_unit_vec);
end

function set_indices = is_point_in_single_cone(points_coords, cone_vertex_pos, cone_axis_unit_vector, cone_max_dist_D, cone_half_angle_rad)
    set_indices = []; num_points = size(points_coords, 1);
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

function P_dir_vecs = findCircleCenter(A_coord, B_coord, O_pos, cone_angle_deg)
    P_dir_vecs = [];
    OA_vec = A_coord - O_pos; 
    OB_vec = B_coord - O_pos;
    if norm(OA_vec) < 1e-9 || norm(OB_vec) < 1e-9, return; end 
    
    OA_unit = OA_vec / norm(OA_vec); 
    OB_unit = OB_vec / norm(OB_vec);
    if abs(dot(OA_unit, OB_unit)) > 1 - 1e-9, return; end
    
    angle_AOB_rad = acos(dot(OA_unit, OB_unit));
    if angle_AOB_rad > deg2rad(cone_angle_deg) + 1e-4, return; end
    
    vec_OD_mid = OA_unit + OB_unit;
    d_od_mid = norm(vec_OD_mid);
    if d_od_mid < 1e-9, return; end
    bisector_dir = vec_OD_mid / d_od_mid;

    angle_P_O_D_mid_rad = acos( cos(deg2rad(cone_angle_deg/2)) / cos(angle_AOB_rad/2) );
    if isnan(angle_P_O_D_mid_rad)
        P_dir_vecs = bisector_dir; 
        return; 
    end
    
    d_oc_len = cos(angle_P_O_D_mid_rad);
    C_point_on_plane = O_pos + bisector_dir * d_oc_len;
    plane_normal_vec = cross(OA_unit, OB_unit);
    if norm(plane_normal_vec) < 1e-9
        P_dir_vecs = bisector_dir;
        return; 
    end
    
    d_cp_len = sin(angle_P_O_D_mid_rad);
    plane_normal_unit_vec = plane_normal_vec / norm(plane_normal_vec);
    
    P1_abs = C_point_on_plane + plane_normal_unit_vec * d_cp_len;
    P2_abs = C_point_on_plane - plane_normal_unit_vec * d_cp_len;
    
    P1_dir = P1_abs - O_pos;
    P2_dir = P2_abs - O_pos;
    
    P_dir_vecs = [P1_dir / norm(P1_dir); P2_dir / norm(P2_dir)];
end

function directions_vec = gcc_direction(A_coord, B_coord, O_pos, cone_angle_deg)
    directions_vec = [];
    OA_vec = A_coord - O_pos; 
    OB_vec = B_coord - O_pos; 
    AB_vec = B_coord - A_coord;
    if norm(OA_vec) < 1e-9 || norm(OB_vec) < 1e-9 || norm(AB_vec) < 1e-9, return; end
    
    e_AB_unit = AB_vec / norm(AB_vec); 
    OA_unit = OA_vec / norm(OA_vec); 
    OB_unit = OB_vec / norm(OB_vec);
    
    sin_half_angle = sin(deg2rad(cone_angle_deg/2));
    dir1 = OA_unit + sin_half_angle * e_AB_unit;
    dir2 = OB_unit - sin_half_angle * e_AB_unit;

    if norm(dir1) > 1e-6, directions_vec = [directions_vec; dir1 / norm(dir1)]; end
    if norm(dir2) > 1e-6, directions_vec = [directions_vec; dir2 / norm(dir2)]; end
end
