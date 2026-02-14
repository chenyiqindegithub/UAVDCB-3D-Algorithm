% AngleSetManager.m

classdef AngleSetManager


    properties
        angle_set = {}
    end

    methods
        function obj = AngleSetManager()
        end

        function obj = add_range(obj, start_angle, end_angle)
            if ~isnumeric(start_angle) || ~isnumeric(end_angle)
                error('起始和结束角度必须是数字。');
            end
            obj.angle_set{end+1} = [double(start_angle), double(end_angle)];
        end

        function merged = get_merged_ranges(obj)
            if isempty(obj.angle_set)
                merged = {};
                return;
            end
            
            intervals = obj.split_wrap_ranges();
            inters = {};
            n = numel(intervals);
            for i = 1:n
                for j = i+1:n
                    a1 = intervals{i}(1); b1 = intervals{i}(2);
                    a2 = intervals{j}(1); b2 = intervals{j}(2);
                    s = max(a1, a2);
                    e = min(b1, b2);
                    if e >= s || abs(e-s) < 1e-9
                        inters{end+1} = [s, e];
                    end
                end
            end

            if isempty(inters)
                if ~isempty(intervals)
                    mat_intervals = cell2mat(intervals');
                    sorted_mat = sortrows(mat_intervals, [1, 2]);
                    merged = mat2cell(sorted_mat, ones(size(sorted_mat, 1), 1), 2)';
                else
                    merged = {};
                end
                return;
            end

            mat_inters = cell2mat(inters');
            sorted_inters = sortrows(mat_inters, [1, 2]);
            
            merged = {};
            if isempty(sorted_inters)
                return;
            end
            
            curr_s = sorted_inters(1, 1);
            curr_e = sorted_inters(1, 2);
            tol = 1e-9;
            
            for i = 2:size(sorted_inters, 1)
                s = sorted_inters(i, 1);
                e = sorted_inters(i, 2);
                if s <= curr_e || abs(s - curr_e) < tol
                    curr_e = max(curr_e, e);
                else
                    merged{end+1} = [curr_s, curr_e];
                    curr_s = s;
                    curr_e = e;
                end
            end
            merged{end+1} = [curr_s, curr_e];
        end

        function angles = get_representative_angles(obj)
            merged = obj.get_merged_ranges();
            if isempty(merged)
                angles = [];
                return;
            end
            
            angles = zeros(1, numel(merged));
            for i = 1:numel(merged)
                range = merged{i};
                angles(i) = obj.calculate_midpoint(range(1), range(2));
            end
        end
    end
    
    methods (Access = private)
        function split_ranges = split_wrap_ranges(obj)
            split_ranges = {};
            tol = 1e-9; % 浮点数比较容差
            for i = 1:numel(obj.angle_set)
                range = obj.angle_set{i};
                start_angle = range(1);
                end_angle = range(2);
                
                s = mod(start_angle, 360);
                e = mod(end_angle, 360);
                
                if start_angle <= end_angle
                    eff_end = e;
                    if abs(e - 0) < tol && end_angle > start_angle
                        eff_end = 360.0;
                    end
                    if abs(end_angle - 360.0) < tol && abs(start_angle - 360.0) > tol
                        eff_end = 360.0;
                    end
                    split_ranges{end+1} = [s, eff_end];
                else
                    split_ranges{end+1} = [s, 360.0];
                    if abs(e - 0) > tol
                        split_ranges{end+1} = [0.0, e];
                    end
                end
            end
        end

        function mid = calculate_midpoint(~, start_angle, end_angle)
            mid_val = (start_angle + end_angle) / 2.0;
            mid = mod(mid_val, 360.0);
        end
    end
end