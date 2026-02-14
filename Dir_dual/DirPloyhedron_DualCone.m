function [node_sets, status_sets, run_time] = DirPloyhedron_DualCone(points, D1, D2, pos, phi1, phi2, varphi_deg, cover_mat, sid)


    run_time = 0;
    tic;

    node_sets = {};
    status_sets = {};

    cover_node_indices = find(cover_mat(:, sid));
    if isempty(cover_node_indices)
        run_time = toc;
        return;
    end
    candidate_points = points(cover_node_indices, :);
    node_num = length(cover_node_indices);
    
    lable = false;
    spot_lable_local_idx = 0;
    pos_to_points = candidate_points - pos;
    in_zone_idx = find(vecnorm(pos_to_points, 2, 2) <= 1.0e-4, 1);
    if ~isempty(in_zone_idx)
        spot_lable_local_idx = in_zone_idx;
        lable = true;
    end

    if lable && node_num == 1
        node_sets{1, 1} = cover_node_indices;
        status_sets{1, 1} = 3;
        run_time = toc;
        return;
    end


    if max(phi1, phi2) > 75
        base_directions = calc_12polytope(pos, 1); 
    else
        base_directions = calc_32polytope(pos, 1);
    end
    

    all_nodes_temp = {};
    all_statuses_temp = {};
    varphi_rad = deg2rad(varphi_deg);
    up_vector = [0, 0, 1];

    for i = 1:size(base_directions, 1)
        e1 = base_directions(i, :);
        e2 = calculate_e2_from_e1_robust(e1, varphi_rad, up_vector);
        if any(isnan(e2)), continue; end
        
        [nodes, statuses] = evaluate_axes_pair([e1; e2], candidate_points, pos, D1, D2, phi1, phi2, cover_node_indices, lable, spot_lable_local_idx);
        if ~isempty(nodes)
            all_nodes_temp{end+1, 1} = nodes;
            all_statuses_temp{end+1, 1} = statuses;
        end
    end
    

    all_currently_covered_nodes = [];
    if ~isempty(all_nodes_temp)
        all_currently_covered_nodes = unique([all_nodes_temp{:}]);
    end

    uncovered_nodes_global_ids = setdiff(cover_node_indices', all_currently_covered_nodes);

    if ~isempty(uncovered_nodes_global_ids)
        for i = 1:length(uncovered_nodes_global_ids)
            node_id_global = uncovered_nodes_global_ids(i);
            local_idx = find(cover_node_indices == node_id_global, 1);
            if isempty(local_idx), continue; end
            
            target_direction_vec = pos_to_points(local_idx, :);
            if norm(target_direction_vec) < 1e-6, continue; end
            target_unit_vec = target_direction_vec / norm(target_direction_vec);


            e1_A = target_unit_vec;
            e2_A = calculate_e2_from_e1_robust(e1_A, varphi_rad, up_vector);
            if ~any(isnan(e2_A))
                [nodes_A, statuses_A] = evaluate_axes_pair([e1_A; e2_A], candidate_points, pos, D1, D2, phi1, phi2, cover_node_indices, lable, spot_lable_local_idx);
                if ~isempty(nodes_A)
                    all_nodes_temp{end+1, 1} = nodes_A;
                    all_statuses_temp{end+1, 1} = statuses_A;
                end
            end
            

            e2_B = target_unit_vec;
            e1_B = calculate_e2_from_e1_robust(e2_B, -varphi_rad, up_vector); % 反向计算 e1
            if ~any(isnan(e1_B))
                [nodes_B, statuses_B] = evaluate_axes_pair([e1_B; e2_B], candidate_points, pos, D1, D2, phi1, phi2, cover_node_indices, lable, spot_lable_local_idx);
                if ~isempty(nodes_B)
                    all_nodes_temp{end+1, 1} = nodes_B;
                    all_statuses_temp{end+1, 1} = statuses_B;
                end
            end
        end
    end


    if ~isempty(all_nodes_temp)
        node_sets = all_nodes_temp;
        status_sets = all_statuses_temp;
    end
    
    run_time = toc;
end


% 辅助函数

function [sorted_nodes, sorted_statuses] = evaluate_axes_pair(axes_pair, candidate_points, pos, D1, D2, phi1, phi2, cover_node_indices, lable, spot_lable_local_idx)
    sorted_nodes = [];
    sorted_statuses = [];
    
    e1 = axes_pair(1,:);
    e2 = axes_pair(2,:);
    
    [set_1, ~] = is_point_in_cone(candidate_points, pos, e1, D1, phi1);
    [set_2, ~] = is_point_in_cone(candidate_points, pos, e2, D2, phi2);
    
    intersect_set = intersect(set_1, set_2);
    only_set_1 = setdiff(set_1, set_2);
    only_set_2 = setdiff(set_2, set_1);
    
    total_set_local_indices = [intersect_set, only_set_1, only_set_2];
    status_vector = [3*ones(1, length(intersect_set)), 1*ones(1, length(only_set_1)), 2*ones(1, length(only_set_2))];
    
    if isempty(total_set_local_indices) && ~lable
        return;
    end
    
    [total_set_local_indices, sort_idx] = sort(total_set_local_indices);
    status_vector = status_vector(sort_idx);
    
    covered_nodes = cover_node_indices(total_set_local_indices)';
    
    if lable
        is_center_in_set = any(total_set_local_indices == spot_lable_local_idx);
        if ~is_center_in_set
            center_node_global_idx = cover_node_indices(spot_lable_local_idx);
            covered_nodes = [covered_nodes, center_node_global_idx];
            status_vector = [status_vector, 3];
            [covered_nodes, sort_idx] = sort(covered_nodes);
            status_vector = status_vector(sort_idx);
        else
            center_node_idx_in_sorted = find(covered_nodes == cover_node_indices(spot_lable_local_idx), 1);
            if ~isempty(center_node_idx_in_sorted)
                status_vector(center_node_idx_in_sorted) = 3;
            end
        end
    end
    
    if isempty(covered_nodes)
        return;
    end
    
    sorted_nodes = covered_nodes;
    sorted_statuses = status_vector;
end

function e2_unit_vec = calculate_e2_from_e1_robust(e1_unit_vec, varphi_rad, up_vector_hint)
    norm_e1 = norm(e1_unit_vec);
    if norm_e1 < 1e-9, e2_unit_vec = [NaN, NaN, NaN]; return; end
    e1_n = e1_unit_vec(:)' / norm_e1;
    k = cross(e1_n, up_vector_hint(:)'); 
    norm_k = norm(k);
    if norm_k < 1e-6
        if abs(e1_n(3)) > 1-1e-6, temp_axis = [1,0,0]; else, temp_axis = [0,0,1]; end
        k = cross(e1_n, temp_axis); norm_k = norm(k);
        if norm_k < 1e-6, if abs(varphi_rad) < 1e-9, e2_unit_vec = e1_n; else, e2_unit_vec = [NaN, NaN, NaN]; end; return; end
    end
    k = k / norm_k; v = e1_n; theta = varphi_rad;
    e2_unit_vec = v*cos(theta) + cross(k,v)*sin(theta) + k*dot(k,v)*(1-cos(theta));
    norm_e2 = norm(e2_unit_vec);
    if norm_e2 < 1e-9, e2_unit_vec = [NaN, NaN, NaN]; else, e2_unit_vec = e2_unit_vec / norm_e2; end
end

function [set,average_angle]=is_point_in_cone(points,charge_pos,charge_direction,line_len,angle_deg)
    set=[]; angleList=[]; average_angle=NaN;
    if isempty(points) || norm(charge_direction) < 1e-9, return; end
    charge_direction = charge_direction / norm(charge_direction);
    for i=1:size(points,1)
        point=points(i,:); vector=point-charge_pos; magnitude_vector = norm(vector);
        if magnitude_vector < 1e-9
             if line_len >= -1e-4, set=[set,i]; angleList=[angleList,0]; end
             continue;
        end
        dotProduct = dot(vector, charge_direction);
        cos_val = dotProduct / magnitude_vector;
        cos_val = max(min(cos_val, 1), -1);
        angleRad = acos(cos_val);
        angleDegrees = rad2deg(angleRad);
        if angleDegrees <= angle_deg/2 + 1e-4 && magnitude_vector <= line_len + 1e-4
            set=[set,i]; angleList=[angleList,angleDegrees];
        end
    end
    if ~isempty(angleList), average_angle=mean(angleList); end
end

function directions=calc_32polytope(pos,radius)
   theta=atan(1/2); p_upper=[]; p_under=[]; A=[pos(1),pos(2),pos(3)+radius]; phi=(1+sqrt(5))/2; L=2*radius/sqrt(phi*sqrt(5));
   for i=1:5
       x=pos(1)+radius*cos(theta)*cos(deg2rad(72)*i); y=pos(2)+radius*cos(theta)*sin(deg2rad(72)*i); z=pos(3)+radius*sin(theta);
       point_upper=[x,y,z]; p_upper=[p_upper;point_upper]; point_under=2*pos-point_upper; p_under=[p_under;point_under];
   end
   distances=pdist2(p_upper,p_under); neighbor_condition_upper=distances-L<1.0e-2; directions=[];
   for i=1:5
      j=mod(i,5)+1; B=p_upper(i,:); C=p_upper(j,:); O=calc_center(A,B,C); direction=O-pos; directions=[directions;direction];
      B=p_under(i,:); C=p_under(j,:); O=calc_center(2*pos-A,B,C); direction=O-pos; directions=[directions;direction];
      neighbor_point=p_under(neighbor_condition_upper(i,:),:); O=calc_center(p_upper(i,:),neighbor_point(1,:),neighbor_point(2,:)); direction=O-pos; directions=[directions;direction];
   end
   directions=[directions;p_under-pos;p_upper-pos;A-pos;pos-A]; directions = directions ./ vecnorm(directions,2,2);
end

function directions=calc_12polytope(pos,radius)
     phi=(1+sqrt(5))/2; d=2*radius/sqrt(3); p_orgen=[d/2,d/2,d/2; d/2,-d/2,d/2; -d/2,d/2,d/2; -d/2,-d/2,d/2]; p_orgen=[p_orgen;-p_orgen] + pos;
     p_green=[0,phi,1/phi; 0,phi,-1/phi; 0,-phi,1/phi; 0,-phi,-1/phi]; p_green = p_green ./ vecnorm(p_green,2,2) * radius + pos;
     p_blue=[1/phi,0,phi; -1/phi,0,phi; 1/phi,0,-phi; -1/phi,0,-phi]; p_blue = p_blue ./ vecnorm(p_blue,2,2) * radius + pos;
     p_pink=[phi,1/phi,0; phi,-1/phi,0; -phi,1/phi,0; -phi,-1/phi,0]; p_pink = p_pink ./ vecnorm(p_pink,2,2) * radius + pos;
     p_points=[p_orgen;p_green;p_blue;p_pink]; L=(4*radius)/(sqrt(3)*(1+sqrt(5))); distances=pdist2(p_points,p_points); neighbor_condition=abs(distances-L)<1.0e-2;
     directions=[];
     for i=1:20
         p=p_points(i,:); neighbor_points=p_points(neighbor_condition(i,:),:);
         if size(neighbor_points, 1) < 2, continue; end
         O=calc_flat_center(p,neighbor_points(1,:),neighbor_points(2,:),L); direction=O-pos; directions=[directions;direction];
     end
     directions = unique(round(directions./vecnorm(directions,2,2), 6), 'rows');
end

function O=calc_flat_center(A,B,C,L)
    AB=B-A; AC=C-A; AD=(AB+AC)/2; theta=deg2rad(54); if norm(AD) < 1e-9, O=A; return; end
    O=A+AD/norm(AD)*L/2/cos(theta);    
end

function O=calc_center(A,B,C)
    D=(B+C)/2; AD=D-A; O=A+2/3*AD;
end