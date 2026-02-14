function [info, tag, duration] = RuncMFEDS(N_nodes, X_chrg_pos, phi1, phi2, D1, D2, varphi, cover_mat, sid, parameter)
    TOL = 1e-9;
    phi1 = deg2rad(phi1 / 2);
    phi2 = deg2rad(phi2 / 2);
    varphi = deg2rad(varphi);

    tStart = tic;
    cover_node_global_indices = find(cover_mat(:,sid));
    if isempty(cover_node_global_indices)
        duration = toc(tStart);
        info = {};
        tag = {};
        return;
    end
    
    info = {};
    tag = {};

    for p_idx = 1:size(X_chrg_pos, 1)
        O = X_chrg_pos(p_idx, :);
        dists_to_O_all = vecnorm(N_nodes - O, 2, 2);
        collocated_node_indices = find(dists_to_O_all < TOL);
        
        non_collocated_nodes_in_sphere_idx = setdiff(cover_node_global_indices, collocated_node_indices);
        
        if isempty(non_collocated_nodes_in_sphere_idx) && ~isempty(collocated_node_indices)
            tmp_info = collocated_node_indices';
            tmp_tag = repmat(3, 1, numel(tmp_info));
            info{end + 1} = tmp_info;
            tag{end + 1} = tmp_tag;
            continue;
        end

        N_sphere = N_nodes(non_collocated_nodes_in_sphere_idx, :);
        
        if isempty(N_sphere)
            continue;
        end

        S_Ant_O_raw = struct('psi_Ant', {}, 'case', {});

        % 孤立节点
        S_Case1_Ant = generate_case1_directions(O, N_sphere, D1, D2, varphi, parameter);
        for i = 1:length(S_Case1_Ant)
            S_Ant_O_raw(end + 1) = struct('psi_Ant', S_Case1_Ant(i), 'case', 'Case1');
        end

        % 一个充电锥的边界接触到 >= 2 个节点
        [S_Ant_AB, S_Ant_C1_no_C2] = generate_case2_directions(O, N_sphere, D1, D2, phi1, phi2, varphi);
        for i = 1:length(S_Ant_AB)
            S_Ant_O_raw(end + 1) = struct('psi_Ant', S_Ant_AB(i), 'case', 'Case2_AB');
        end
        for i = 1:length(S_Ant_C1_no_C2)
            ea = S_Ant_C1_no_C2(i).ea;
            eb = rotate_vector(ea, varphi);
            S_Ant_O_raw(end + 1) = struct('psi_Ant', struct('ea', ea, 'eb', eb), 'case', 'Case2_AB_NoC2');
        end

        [S_Ant_BA_swapped, S_Ant_C2_no_C1] = generate_case2_directions(O, N_sphere, D2, D1, phi2, phi1, varphi);
        for i = 1:length(S_Ant_BA_swapped)
            psi.ea = S_Ant_BA_swapped(i).eb; psi.eb = S_Ant_BA_swapped(i).ea;
            S_Ant_O_raw(end + 1) = struct('psi_Ant', psi, 'case', 'Case2_BA');
        end
        for i = 1:length(S_Ant_C2_no_C1)
            eb = S_Ant_C2_no_C1(i).ea;
            ea = rotate_vector(eb, varphi);
            S_Ant_O_raw(end + 1) = struct('psi_Ant', struct('ea', ea, 'eb', eb), 'case', 'Case2_BA_NoC1');
        end


        S_Case3_Ant = generate_case3_directions(O, N_sphere, D1, D2, phi1, phi2, varphi);
        for i = 1:length(S_Case3_Ant)
            S_Ant_O_raw(end + 1) = struct('psi_Ant', S_Case3_Ant(i), 'case', 'Case3');
        end
        
        if isempty(S_Ant_O_raw) && ~isempty(collocated_node_indices)
            continue;
        elseif isempty(S_Ant_O_raw)
            continue;
        end

        candidate_entries = {};
        for i = 1:length(S_Ant_O_raw)
            psi = S_Ant_O_raw(i).psi_Ant;
            covered_idx_c1 = find_covered_node_indices(O, N_nodes, psi.ea, D1, phi1);
            covered_idx_c2 = find_covered_node_indices(O, N_nodes, psi.eb, D2, phi2);
            covered_union = union(covered_idx_c1, covered_idx_c2);

            if isempty(covered_union), continue; end

            nodes_in_c1_only = setdiff(covered_idx_c1, covered_idx_c2);
            nodes_in_c2_only = setdiff(covered_idx_c2, covered_idx_c1);
            nodes_in_both = intersect(covered_idx_c1, covered_idx_c2);
            
            current_efficiency = calculate_total_efficiency(dists_to_O_all, ...
                                                            nodes_in_c1_only, ...
                                                            nodes_in_c2_only, ...
                                                            nodes_in_both, ...
                                                            parameter);
            
            candidate_entries(end+1, :) = {covered_union, psi, current_efficiency};
        end

        if isempty(candidate_entries)
            continue;
        end
        
        final_selection_cell = filter_sets_code1_style(candidate_entries);
        
        num_final = size(final_selection_cell, 1);

        for i = 1:num_final
            final_union_raw = final_selection_cell{i, 1};
            final_psi_ant   = final_selection_cell{i, 2};

            if ~isempty(collocated_node_indices)
                final_union_with_collocated = union(final_union_raw, collocated_node_indices);
            else
                final_union_with_collocated = final_union_raw;
            end
            
            tmp_info = [];
            tmp_tag = [];

            for j = 1:numel(final_union_with_collocated)
                node_id = final_union_with_collocated(j);
                tmp_info(end + 1) = node_id;
                xyz = N_nodes(node_id, :);
                
                if ismember(node_id, collocated_node_indices)
                    tmp_tag(end + 1) = 3;
                else
                    in_c1 = is_point_in_cone(O, xyz, final_psi_ant.ea, D1, phi1);
                    in_c2 = is_point_in_cone(O, xyz, final_psi_ant.eb, D2, phi2);
                    if in_c1 && in_c2, tmp_tag(end + 1) = 3;
                    elseif in_c1, tmp_tag(end + 1) = 1;
                    elseif in_c2, tmp_tag(end + 1) = 2;
                    else, tmp_tag(end + 1) = 0;
                    end
                end
            end
            
            info{end + 1} = tmp_info;
            tag{end + 1} = tmp_tag;
        end
    end

    info = convertCellToColumn(info);
    tag = convertCellToColumn(tag);

    duration = toc(tStart);
end

% 情况1 辅助函数
function S_Case1_Ant = generate_case1_directions(O, N_nodes, D1, D2, varphi, parameter)
    S_Case1_Ant = struct('ea', {}, 'eb', {});

    for i = 1:size(N_nodes, 1)
        A = N_nodes(i, :);

        if A == O
            fprintf('case1出错：【出现与原点重合的节点】\n');
            continue;
        end

        e_OA = normalize_local(A - O);
        dist_OA = norm(A - O);

        if dist_OA <= D1 && dist_OA > D2
            S_Case1_Ant(end + 1).ea = e_OA;
            eb = rotate_vector(e_OA, varphi);
            S_Case1_Ant(end).eb = eb;
        elseif dist_OA > D1 && dist_OA <= D2
            eb = rotate_vector(e_OA, varphi);
            S_Case1_Ant(end + 1).ea = eb;
            S_Case1_Ant(end).eb = e_OA;
        elseif dist_OA <= D1 && dist_OA <= D2
            if (parameter.alpha1 / (dist_OA + parameter.beta1)^parameter.gamma1) > (parameter.alpha2 / (dist_OA + parameter.beta2)^parameter.gamma2)
                S_Case1_Ant(end + 1).ea = e_OA;
                eb = rotate_vector(e_OA, varphi);
                S_Case1_Ant(end).eb = eb;
            else
                eb = rotate_vector(e_OA, varphi);
                S_Case1_Ant(end + 1).ea = eb;
                S_Case1_Ant(end).eb = e_OA;
            end
        end
    end
end

% 情况2 辅助函数
function [S_Ant, S_Dir_C1_no_C2] = generate_case2_directions(O, N_sphere, D1, D2, phi1, phi2, varphi)
    TOL = 1e-9;
    S_Ant = struct('ea', {}, 'eb', {});
    S_Dir_C1_no_C2 = struct('ea', {});

    num_points = size(N_sphere, 1);

    if num_points < 2
        return;
    end

    for idx_A = 1:num_points
        for idx_B = 1:num_points
            if idx_A == idx_B
                continue;
            end

            A = N_sphere(idx_A, :);
            B = N_sphere(idx_B, :);

            if isequal(A, O) || isequal(B, O)
                fprintf('case2出错：【出现与原点重合的节点】\n');
                continue;
            end

            dist_OA = norm(A - O);
            dist_OB = norm(B - O);

            if dist_OA > D1 + TOL || dist_OB > D1 + TOL
                continue;
            end

            e_OA = normalize_local(A - O);
            e_OB = normalize_local(B - O);
            varphi_AOB = angle_between_vectors(e_OA, e_OB);

            if varphi_AOB > 2 * phi1 + TOL
                continue;
            end

            [e1, e2] = compute_v_AB(e_OA, e_OB, phi1);

            if isempty(e1) || isempty(e2)
                continue;
            end

            ea = e1;

            angle_set = AngleSetManager();
            container = AngleRangeContainer();
            e_REF = random_perpendicular_unit_vector(ea);

            indices_to_exclude = [idx_A, idx_B];
            mask = true(num_points, 1);
            mask(indices_to_exclude) = false;
            other_indices = find(mask);

            for j = 1:length(other_indices)
                idx_C = other_indices(j);
                C = N_sphere(idx_C, :);

                if isequal(C, O)
                    fprintf('case2出错：【出现与原点重合的节点】\n');
                    continue;
                end

                dist_OC = norm(C - O);

                if dist_OC > D2 + TOL || dist_OC < TOL
                    continue;
                end

                e_OC = normalize_local(C - O);

                [e_c2_dir1, e_c2_dir2] = compute_e_c2(ea, e_OC, angle_between_vectors(ea, e_OC), varphi, phi2);

                if isempty(e_c2_dir1) || isempty(e_c2_dir2)
                    continue;
                end

                angle_1 = compute_angle_with_e_REF(ea, e_REF, varphi, O, e_c2_dir1);
                angle_2 = compute_angle_with_e_REF(ea, e_REF, varphi, O, e_c2_dir2);

                angleMin = min(angle_1, angle_2) + pi;
                angleMax = max(angle_1, angle_2) + pi;

                tmpAngle = (angleMin + angleMax) / 2;
                tmpV = compute_e2_with_angle(ea, e_REF, varphi, O, tmpAngle - pi);

                if is_point_in_cone(O, C, tmpV, D1, phi1)
                    angle_set.add_range(rad2deg(angleMin), rad2deg(angleMax));
                    container.add_range(idx_C, rad2deg(angleMin), rad2deg(angleMax));
                else
                    angle_set.add_range(rad2deg(angleMax), 360);
                    container.add_range(idx_C, rad2deg(angleMax), 360);
                    angle_set.add_range(0, rad2deg(angleMin));
                    container.add_range(idx_C, 0, rad2deg(angleMin));
                end

            end

            representatives1 = angle_set.get_representative_angles();

            for ang = 1:length(representatives1)
                result = container.query_angle(representatives1(ang));

                if result.size == 0
                    continue
                end

                eb = compute_e2_with_angle(ea, e_REF, varphi, O, representatives1(ang) - pi);
                S_Ant(end + 1).ea = ea;
                S_Ant(end).eb = eb;
            end

            S_Dir_C1_no_C2(end + 1).ea = ea;
        end
    end
end


function S_Case3_Ant = generate_case3_directions(O, N_sphere, D1, D2, phi1, phi2, varphi)
    TOL = 1e-9;
    S_Case3_Ant = struct('ea', {}, 'eb', {});

    if size(N_sphere, 1) < 2
        return;
    end

    num_nodes = size(N_sphere, 1);
    temp_results = {};

    for i = 1:num_nodes

        for j = 1:num_nodes
            if i == j, continue; end

            A = N_sphere(i, :);
            B = N_sphere(j, :);

            if isequal(A, O) || isequal(B, O)
                fprintf('case3出错：【出现与原点重合的节点】\n');
                continue;
            end

            dist_OA = norm(A - O);
            dist_OB = norm(B - O);

            if dist_OA > max(D1, D2) + TOL || dist_OB > max(D1, D2) + TOL
                fprintf('case3出错：【出现覆盖不到的节点】\n');
                continue;
            end

            e_OA = normalize_local(A - O);
            e_OB = normalize_local(B - O);

            e_OA = e_OA / norm(e_OA);
            e_OB = e_OB / norm(e_OB);

            angle_AOB = angle_between_vectors(e_OA, e_OB);

            if angle_AOB + phi1 > varphi + phi2 || angle_AOB + phi2 > varphi + phi1
                [e1_dir1, ~] = compute_e_c2(e_OB, e_OA, angle_AOB, varphi + phi2, phi1);
                if ~isempty(e1_dir1)
                    e1_T = e1_dir1;
                    e2_T = rotate_vector_in_plane(e1_T, e_OB, varphi);
                    if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                        temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                    end
                end

                [e2_dir1, ~] = compute_e_c2(e_OA, e_OB, angle_AOB, varphi + phi1, phi2);
                if ~isempty(e2_dir1)
                    e2_T = e2_dir1;
                    e1_T = rotate_vector_in_plane(e2_T, e_OA, varphi);
                    if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                        temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                    end
                end
            end

            if angle_AOB + phi1 <= varphi + phi2
                e1_T_candidate = rotate_vector_in_plane(e_OB, e_OA, angle_AOB + phi1);
                if ~isempty(e1_T_candidate)
                    [e2_dir1, ~] = compute_e_c2(e1_T_candidate, e_OB, angle_between_vectors(e1_T_candidate, e_OB), varphi, phi2);
                    if ~isempty(e2_dir1)
                        e1_T = e1_T_candidate;
                        e2_T = e2_dir1;
                        if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                            temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                        end
                    end
                end
            end

            if angle_AOB + phi2 <= varphi + phi1
                e2_T_candidate = rotate_vector_in_plane(e_OA, e_OB, angle_AOB + phi2);
                if ~isempty(e2_T_candidate)
                    [e1_dir1, ~] = compute_e_c2(e2_T_candidate, e_OA, angle_between_vectors(e2_T_candidate, e_OA), varphi, phi1);
                    if ~isempty(e1_dir1)
                        e2_T = e2_T_candidate;
                        e1_T = e1_dir1;
                        if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                            temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                        end
                    end
                end
            end

            if abs(angle_AOB - phi1) >= abs(varphi - phi2)
                e1_T_candidate = rotate_vector_in_plane(e_OA, e_OB, phi1);
                if ~isempty(e1_T_candidate)
                    [e2_dir1, ~] = compute_e_c2(e1_T_candidate, e_OB, angle_between_vectors(e1_T_candidate, e_OB), varphi, phi2);
                    if ~isempty(e2_dir1)
                        e1_T = e1_T_candidate;
                        e2_T = e2_dir1;
                        if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                            temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                        end
                    end

                    [e2_dir1_alt, ~] = compute_e_c2(e1_T_candidate, e_OB, phi1 - angle_AOB, varphi, phi2);
                    if ~isempty(e2_dir1_alt)
                        e1_T = e1_T_candidate;
                        e2_T = e2_dir1_alt;
                        if ~isempty(e1_T) && ~isempty(e2_T) && numel(e1_T) == 3 && numel(e2_T) == 3
                            temp_results{end + 1} = [e1_T(:)', e2_T(:)'];
                        end
                    end
                end
            end
        end
    end

    if ~isempty(temp_results)
        psi_Ant_mat = cell2mat(temp_results');
        [~, unique_idx] = unique(round(psi_Ant_mat / TOL) * TOL, 'rows');
        unique_mat = psi_Ant_mat(unique_idx, :);

        for k = 1:size(unique_mat, 1)
            S_Case3_Ant(k).ea = unique_mat(k, 1:3);
            S_Case3_Ant(k).eb = unique_mat(k, 4:6);
        end
    end
end

% 通用辅助函数
function v_norm = normalize_local(v)
    tol = 1e-9;
    v = v(:);
    n = norm(v);
    if n < tol
        v_norm = [];
    else
        v_norm = v / n;
    end
end

function angle = angle_between_vectors(v1, v2)
    v1 = v1(:);
    v2 = v2(:);
    v1_norm = normalize_local(v1);
    v2_norm = normalize_local(v2);

    if isempty(v1_norm) || isempty(v2_norm)
        error('输入向量为零向量，无法计算夹角');
    end

    if length(v1_norm) ~= length(v2_norm)
        error('v1 和 v2 的长度必须相同');
    end

    dot_product = dot(v1_norm, v2_norm);
    dot_product = max(-1, min(1, dot_product));
    angle = acos(dot_product);
end

function v_norm = normalize_vec(v)
    TOL = 1e-9;
    v = reshape(v, [], 3);
    norm_v = vecnorm(v, 2, 2);
    valid_rows = norm_v > TOL;
    
    v_norm = v;
    if any(valid_rows)
        v_norm(valid_rows, :) = v(valid_rows, :) ./ norm_v(valid_rows);
    end
end

function covered_indices = find_covered_node_indices(O, all_nodes, e_dir, D, phi)
    TOL = 1e-9;
    vecs_from_O = all_nodes - O;
    dists = vecnorm(vecs_from_O, 2, 2);
    nodes_in_dist_mask = (dists <= D + TOL);

    if ~any(nodes_in_dist_mask)
        covered_indices = [];
        return;
    end

    nodes_in_dist_idx = find(nodes_in_dist_mask);
    vecs_to_nodes_in_dist = normalize_vec(vecs_from_O(nodes_in_dist_idx, :));
    e_dir_norm = normalize_vec(e_dir);

    dot_products = vecs_to_nodes_in_dist * e_dir_norm';
    dot_products_clamped = min(1, max(-1, dot_products));
    angles = acos(dot_products_clamped);
    nodes_in_angle_mask = (angles <= phi + TOL);
    covered_indices = nodes_in_dist_idx(nodes_in_angle_mask);
end

function [e1, e2] = compute_v_AB(ea, eb, phi1)
    TOL = 1e-9;
    ea = ea(:)';
    eb = eb(:)';
    varphi_AOB = acos(min(1, max(-1, dot(ea, eb))));
    if varphi_AOB > 2 * phi1 + TOL
        e1 = [];
        e2 = [];
        return;
    end

    tan_phi_sq = tan(phi1) ^ 2;
    tan_varphi_half_sq = tan(varphi_AOB / 2) ^ 2;
    if tan_phi_sq < tan_varphi_half_sq - TOL
        e1 = [];
        e2 = [];
        return;
    end

    cos_varphi_AOB = cos(varphi_AOB);
    term1 = (ea + eb) / (cos_varphi_AOB + 1);
    cross_vec = cross(eb, ea);
    sqrt_term = sqrt(max(0, tan_phi_sq - tan_varphi_half_sq));
    sin_varphi_AOB = sin(varphi_AOB);

    if abs(sin_varphi_AOB) < TOL
        e1 = [];
        e2 = [];
        return;
    end

    common_factor = 1 / sqrt(1 + tan_phi_sq);
    e1_unnorm = term1 + (sqrt_term / sin_varphi_AOB) * cross_vec;
    e2_unnorm = term1 - (sqrt_term / sin_varphi_AOB) * cross_vec;
    e1 = common_factor * normalize_local(e1_unnorm);
    e2 = common_factor * normalize_local(e2_unnorm);
end

function v_perp_unit = random_perpendicular_unit_vector(e_a)
    if ~(isvector(e_a) && numel(e_a) == 3)
        error('输入必须是一个3维向量 (1x3 或 3x1)。');
    end

    is_row_vector = isrow(e_a);
    if is_row_vector
        e_a = e_a';
    end

    e_a_unit = e_a / norm(e_a);

    while true
        v_rand = randn(3, 1);
        v_perp = v_rand - dot(v_rand, e_a_unit) * e_a_unit;
        norm_v_perp = norm(v_perp);

        if norm_v_perp > 1e-8
            v_perp_unit = v_perp / norm_v_perp;
            break;
        end
    end

    if is_row_vector
        v_perp_unit = v_perp_unit';
    end
end

function [e_c2_dir1, e_c2_dir2] = compute_e_c2(e_a, e_b, phi_AOB, psi, phi2)
    TOL = 1e-9;
    e_c2_dir1 = [];
    e_c2_dir2 = [];

    if ~isvector(e_a) || numel(e_a) ~= 3 || ~isvector(e_b) || numel(e_b) ~= 3
        error('输入向量 e_a 和 e_b 必须是3维向量。');
    end

    e_a = e_a(:);
    e_b = e_b(:);
    e_a = normalize_local(e_a);
    e_b = normalize_local(e_b);

    if isempty(e_a) || isempty(e_b)
        return;
    end

    sin_phi_AOB = sin(phi_AOB);
    sin_psi = sin(psi);

    if abs(sin_phi_AOB) < TOL || abs(sin_psi) < TOL
        return;
    end

    cos_term_num = cos(phi2) - cos(phi_AOB) * cos(psi);
    cos_term_den = sin_phi_AOB * sin_psi;
    cos_term_val = cos_term_num / cos_term_den;

    if cos_term_val < -1.0 - TOL || cos_term_val > 1.0 + TOL
        return;
    end

    cos_term_val = max(-1.0, min(1.0, cos_term_val));
    phi_ab = acos(cos_term_val);

    cos_phi_AOB = cos(phi_AOB);
    ctan_phi_AOB = cos_phi_AOB / sin_phi_AOB;
    csc_phi_AOB = 1.0 / sin_phi_AOB;

    e_perp = cross(e_a, e_b);
    norm_e_perp = norm(e_perp);

    if norm_e_perp < TOL
        return;
    end

    term1 = cos(psi) - sin(psi) * cos(phi_ab) * ctan_phi_AOB;
    term2 = sin(psi) * cos(phi_ab) * csc_phi_AOB;
    term3 = sin(psi) * sin(phi_ab);

    v_c2_1 = term1 * e_a + term2 * e_b + term3 * csc_phi_AOB * e_perp;
    v_c2_2 = term1 * e_a + term2 * e_b - term3 * csc_phi_AOB * e_perp;

    e_c2_dir1 = normalize_local(v_c2_1);
    e_c2_dir2 = normalize_local(v_c2_2);
end

function angle = compute_angle_with_e_REF(e1, e_REF, phi, O, e_2)
    TOL = 1e-8;

    if nargin ~= 5
        error('需要 5 个输入参数: e1, e_REF, phi, O, e_2');
    end

    e1 = e1(:);
    e_REF = e_REF(:);
    e_2 = e_2(:);
    O = O(:);

    if numel(e1) ~= 3 || numel(e_REF) ~= 3 || numel(e_2) ~= 3 || numel(O) ~= 3
        error('所有输入向量 (e1, e_REF, e_2, O) 的形状必须是 (3,1)');
    end

    e1 = e1 / norm(e1);
    e_REF = e_REF / norm(e_REF);
    e_2 = e_2 / norm(e_2);

    if abs(dot(e1, e_REF)) > TOL
        error('e1 和 e_REF 必须正交。点积为: %e', dot(e1, e_REF));
    end

    verAng = angle_between_vectors(e1, e_2);

    if (verAng - phi) > TOL
        warning('MATLAB:AngleMismatch', ...
            'e1 和 e_2 之间的角度 (%.6f rad) 与提供的 phi (%.6f rad) 不匹配。', ...
            verAng, phi);
    end

    if abs(phi) < TOL || abs(phi - pi) < TOL
        angle = 0.0;
        return;
    end

    e_perp_basis = cross(e1, e_REF);
    x_coord = dot(e_2, e_REF);
    y_coord = dot(e_2, e_perp_basis);
    angle = atan2(y_coord, x_coord);
end

function e2 = compute_e2_with_angle(e1, e_REF, phi, O, angle)
    TOL = 1e-8;

    if nargin ~= 5
        error('需要 5 个输入参数: e1, e_REF, phi, O, angle');
    end

    e1 = e1(:);
    e_REF = e_REF(:);
    O = O(:);

    if numel(e1) ~= 3 || numel(e_REF) ~= 3 || numel(O) ~= 3
        error('所有输入向量 (e1, e_REF, O) 的形状必须是 (3,1)');
    end

    e1 = e1 / norm(e1);
    e_REF = e_REF / norm(e_REF);

    if abs(dot(e1, e_REF)) > TOL
        error('e1 和 e_REF 必须正交。点积为: %e', dot(e1, e_REF));
    end

    cos_phi = cos(phi);
    sin_phi = sin(phi);
    cos_angle = cos(angle);
    sin_angle = sin(angle);

    e2_parallel = cos_phi * e1;
    e_perp_basis = cross(e1, e_REF);
    e2_perp = sin_phi * (cos_angle * e_REF + sin_angle * e_perp_basis);
    e2 = e2_parallel + e2_perp;
    e2 = e2 / norm(e2);
end

function isIn = is_point_in_cone(pos, point, direction, D, phi)
    TOL = 1e-9;

    pos = pos(:);
    point = point(:);
    direction = direction(:);

    to_point = point - pos;
    dist = norm(to_point);

    if dist < TOL
        isIn = true;
        return;
    end

    direction_norm = direction / norm(direction);
    to_point_norm = to_point / dist;

    dot_prod = dot(to_point_norm, direction_norm);
    dot_prod_clamped = max(-1.0, min(1.0, dot_prod));
    angle = acos(dot_prod_clamped);

    if angle > phi + TOL
        isIn = false;
        return;
    end

    if dist > D + TOL
        isIn = false;
        return;
    end

    isIn = true;
end

function e_oc = rotate_vector_in_plane(e_oa, e_ob, phi_aoc)
    TOL = 1e-9;

    e_oa = e_oa(:);
    e_ob = e_ob(:);

    if norm(e_oa) < TOL
        error('输入向量 e_oa 不能是零向量以进行旋转。');
    end

    if norm(e_ob) < TOL
        error('输入向量 e_ob 不能是零向量以定义旋转平面。');
    end

    e_oa_unit = e_oa / norm(e_oa);
    e_ob_unit = e_ob / norm(e_ob);

    dot_product = dot(e_oa_unit, e_ob_unit);
    dot_product = max(-1.0, min(1.0, dot_product));

    v_perp_unnormalized = e_ob_unit - dot_product * e_oa_unit;
    norm_perp = norm(v_perp_unnormalized);

    if norm_perp < TOL
        [~, min_idx] = min(abs(e_oa_unit));
        dummy_vec = zeros(3, 1);
        dummy_vec(min_idx) = 1;
        e_perp = cross(e_oa_unit, dummy_vec);
        e_perp = e_perp / norm(e_perp);
    else
        e_perp = v_perp_unnormalized / norm_perp;
    end

    cos_phi = cos(phi_aoc);
    sin_phi = sin(phi_aoc);

    v_oc_unnormalized = cos_phi * e_oa_unit + sin_phi * e_perp;
    e_oc_norm = norm(v_oc_unnormalized);

    if e_oc_norm < TOL
        e_oc = zeros(3, 1);
    else
        e_oc = v_oc_unnormalized / e_oc_norm;
    end
end

function [selected_entries] = greedySetCover(all_entries)
    num_entries = length(all_entries);

    if num_entries == 0
        selected_entries = [];
        return;
    end

    all_unions_cell = cell(1, num_entries);

    for i = 1:num_entries
        union_vector = all_entries(i).Covered_Nodes.Union;
        all_unions_cell{i} = union_vector(:).';
    end

    all_nodes_vector = [all_unions_cell{:}];

    if isempty(all_nodes_vector)
        selected_entries = [];
        return;
    end

    nodes_to_be_covered = unique(all_nodes_vector);

    candidate_indices = 1:num_entries;
    selected_indices = [];

    while ~isempty(nodes_to_be_covered)
        best_candidate_local_idx = -1;
        max_coverage_count = -1;

        for i = 1:length(candidate_indices)
            current_entry_idx = candidate_indices(i);
            current_union_set = all_entries(current_entry_idx).Covered_Nodes.Union;

            newly_covered = intersect(current_union_set, nodes_to_be_covered);
            coverage_count = length(newly_covered);

            if coverage_count > max_coverage_count
                max_coverage_count = coverage_count;
                best_candidate_local_idx = i;
            end
        end

        if best_candidate_local_idx == -1
            warning('存在无法被覆盖的节点，算法提前终止。');
            break;
        end

        best_entry_original_idx = candidate_indices(best_candidate_local_idx);
        selected_indices = [selected_indices, best_entry_original_idx];
        best_union_set = all_entries(best_entry_original_idx).Covered_Nodes.Union;
        nodes_to_be_covered = setdiff(nodes_to_be_covered, best_union_set);
        candidate_indices(best_candidate_local_idx) = [];
    end

    selected_entries = all_entries(selected_indices);
end

function v_rotated = rotate_vector(v, theta, axis)
    is_row = isrow(v);
    v = v(:);

    if nargin < 3 || isempty(axis)
        axis = [0, 0, 1]';
    else
        axis = axis(:);
    end

    axis = axis / norm(axis);

    k = axis;
    cos_theta = cos(theta);
    sin_theta = sin(theta);

    term1 = v * cos_theta;
    term2 = cross(k, v) * sin_theta;
    term3 = k * (k' * v) * (1 - cos_theta);

    v_rotated = term1 + term2 + term3;

    if is_row
        v_rotated = v_rotated';
    end
end

function columnCell = convertCellToColumn(rowCell)
    if ~iscell(rowCell)
        error('输入必须是一个 cell 数组。');
    end
    columnCell = rowCell(:);
end

function total_eff = calculate_total_efficiency(all_dists, nodes_c1_only, nodes_c2_only, nodes_both, parameter)
    total_eff = 0;

    alpha1 = parameter.alpha1; beta1 = parameter.beta1; gamma1 = parameter.gamma1;
    alpha2 = parameter.alpha2; beta2 = parameter.beta2; gamma2 = parameter.gamma2;

    if ~isempty(nodes_c1_only)
        dists_c1 = all_dists(nodes_c1_only);
        effs_c1 = alpha1 ./ (dists_c1 + beta1).^gamma1;
        total_eff = total_eff + sum(effs_c1);
    end

    if ~isempty(nodes_c2_only)
        dists_c2 = all_dists(nodes_c2_only);
        effs_c2 = alpha2 ./ (dists_c2 + beta2).^gamma2;
        total_eff = total_eff + sum(effs_c2);
    end

    if ~isempty(nodes_both)
        dists_both = all_dists(nodes_both);
        effs1_for_both = alpha1 ./ (dists_both + beta1).^gamma1;
        effs2_for_both = alpha2 ./ (dists_both + beta2).^gamma2;
        max_effs = max(effs1_for_both, effs2_for_both);
        total_eff = total_eff + sum(max_effs);
    end
end

function final_entries = filter_sets_code1_style(candidate_entries)
    final_entries = {};
    set_size = size(candidate_entries, 1);

    if set_size == 0
        return;
    end

    set_len = zeros(set_size, 1);
    for num = 1:set_size
        set_len(num) = length(candidate_entries{num, 1});
    end
    [~, sorted_indices] = sort(set_len, 'descend');
    sorted_entries = candidate_entries(sorted_indices, :);

    processed_rows = false(set_size, 1);

    for i = 1:set_size
        if processed_rows(i)
            continue;
        end

        set_i = sort(sorted_entries{i, 1});
        best_efficiency_for_set_i = sorted_entries{i, 3};
        best_entry_index = i;

        for j = i + 1:set_size
            if processed_rows(j)
                continue;
            end

            set_j = sort(sorted_entries{j, 1});

            if isequal(set_i, set_j)
                processed_rows(j) = true;
                if sorted_entries{j, 3} > best_efficiency_for_set_i
                    best_efficiency_for_set_i = sorted_entries{j, 3};
                    best_entry_index = j;
                end
            elseif all(ismember(set_j, set_i))
                 processed_rows(j) = true;
            end
        end
        
        final_entries(end + 1, :) = sorted_entries(best_entry_index, :);
    end
end
