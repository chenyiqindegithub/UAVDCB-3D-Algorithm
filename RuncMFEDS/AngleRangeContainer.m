classdef AngleRangeContainer
    properties
        % 使用结构体数组来存储带标识符的范围
        angle_ranges = struct('identifier', {}, 'start_angle', {}, 'end_angle', {});
    end
    
    methods
        function obj = AngleRangeContainer()
            % 构造函数
        end
        
        function obj = add_range(obj, identifier, start_angle, end_angle)
            % 添加一个新的角度范围
            new_range.identifier = identifier;
            new_range.start_angle = start_angle;
            new_range.end_angle = end_angle;
            
            obj.angle_ranges(end+1) = new_range;
        end
        
        function matching_identifiers = query_angle(obj, angle)
            
            matching_identifiers = {};
            
            for i = 1:numel(obj.angle_ranges)
                range = obj.angle_ranges(i);
                if angle >= range.start_angle && angle <= range.end_angle
                    matching_identifiers{end+1} = range.identifier;
                end
            end
        end
    end
end