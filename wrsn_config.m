classdef wrsn_config
    properties
        % node_num : 节点个数
        node_num
        % x_coords:空间范围参数 X
        x_coords
        % y_coords:空间范围参数 Y
        y_coords
        % z_coords: 空间范围参数 Z
        z_coords
        % radius:充电半径
        radius
        radius2
        angle
        angle2
        angle_between
        % alpha: 能量传输效率参数
        alpha
        alpha2
        % beta:能量传输效率参数
        beta
        beta2
        % gamma :能量传输效率参数
        gamma
        gamma2
        % E_orgmin :节点起始能量最小值
        E_orgmin
        % E_orgmax :节点起始能量最大值
        E_orgmax
        % E_expmin :节点期望能量最小值
        E_expmin
        % E_expmax :节点期望能量最大值
        E_expmax
        % n_need_ratio :需要捕获能量的节点比例
        n_need_ratio
        % egy_increase :能量补充值
        egy_increase
        % p :充电功率
        p
        % E_orgV :充电器的满电能量
        E_orgV
        % E_orgL : 充电器的最低能量
        E_orgL
        %充电器完成充电后能量
        Energy  
        % move_cost_rate: 移动能耗
        move_cost_rate
        % hover_cost_rate: 悬停能耗
        hover_cost_rate
        % move_v :飞行速度
        move_v  
        % pos: 无人机基站
        pos
        %mode_p: 生成充电位置的方式
        mode_p
        %mode_d : 充电阈形式
        mode_d
        %mode_tsp:求解TSP问题的方式
        mode_tsp
        %points :节点的坐标
        points
        % spots :充电位置坐标
        spots
        %E_org :节点的起始能量
        E_org
        %E_exp :节点的期望能量
        E_exp
        % list_need_egy :期望能量的节点位置
        list_need_egy
        % wrsn的初始总能量
        init_total_energy
        c_max
        center_num
        isDense
        cover_mat
        rate_mat
        rate_mat3D
        %path_loss 路径能耗
        path_loss
        %hover_loss 悬停能耗
        hover_loss
        %charge_loss 充电能耗
        charge_loss
        dist_len
        dist_mat
        move_t
        tsp_time
        charge_t
        cost_t
        cost_mat
        cone_number
        cdir_time
        lp_time
        entire_time
    end
    methods
        function wrsn = wrsn_config(node_num,radius,radius2,angle,angle2,angle_between)
            wrsn.mode_p=1;
            wrsn.mode_d=1;
            wrsn.mode_tsp=1;
            wrsn.isDense=0;
            wrsn.node_num = node_num;
            wrsn.x_coords=100;
            wrsn.y_coords=100;
            wrsn.z_coords=20;
            wrsn.alpha=12;
            wrsn.alpha2=12;
            wrsn.beta=4;
            wrsn.beta2=4;
            wrsn.gamma=2;
            wrsn.gamma2=2;
            wrsn.radius=radius;
            wrsn.radius2=radius2;
            wrsn.angle=angle;
            wrsn.angle2=angle2;
            wrsn.angle_between=angle_between;
            wrsn.E_orgmin=20; %单位 J
            wrsn.E_orgmax=90;%单位 J
            wrsn.E_expmin=20; %单位 J
            wrsn.E_expmax=90;% 单位 J
            wrsn.egy_increase=10; %单位 J
            wrsn.p=1; %单位 J/s
            wrsn.move_cost_rate=4; %单位 J/s
            wrsn.hover_cost_rate=4; %单位 J/s
            wrsn.move_v=1;%单位 m/s
            wrsn.n_need_ratio=1;
            wrsn.E_orgV=2500000;
            wrsn.Energy=wrsn.E_orgV;
            wrsn.E_orgL=100;
            wrsn.pos=[0,0,0];
            wrsn.c_max=0.9;
            wrsn.center_num=100;
            wrsn.entire_time=0;
            wrsn.lp_time=0;
            wrsn.cdir_time=0;
            points=generateNodes(wrsn);
            wrsn.points=points;
            wrsn.spots=[];
            wrsn.cover_mat=[];
            wrsn.rate_mat=[];
            wrsn.rate_mat3D=[];
            [E_org,E_exp,list_need_egy]=generateNodeInitEnergy(wrsn);
            wrsn.E_org=E_org;
            wrsn.E_exp=E_exp;
            wrsn.list_need_egy=list_need_egy;
            wrsn.init_total_energy=sum(E_org)+wrsn.E_orgV;
        end
        function points=generateNodes(wrsn)
            if wrsn.isDense==1
                points=generateNodesByDense(wrsn);
            elseif wrsn.isDense == 0
                points=generateNodeByRandom(wrsn);
            else
                points = generateNodesBySparse(wrsn);
            end
        end
        function points=generateNodeByRandom(wrsn)
            xmin = 0; xmax = wrsn.x_coords;
            ymin = 0; ymax = wrsn.y_coords;
            zmin = 0; zmax = wrsn.z_coords;
            n = wrsn.node_num;
            points=zeros(n,3);
            i = 1;
                while i <= n
                    rng('shuffle');
                    x=randi([xmin,xmax]);
                    y=randi([ymin,ymax]);
                    z=randi([zmin,zmax]);
                    newNode = [x, y, z];
                    if ~ismember(newNode, points, 'rows')
                        points(i, :) = newNode;
                        i = i + 1;
                    end
                 end
        end
        %生成密集节点位置
        function points=generateNodesByDense(wrsn)
            X=wrsn.radius:wrsn.x_coords-wrsn.radius;
            Y=wrsn.radius:wrsn.y_coords-wrsn.radius;
            Z=wrsn.radius:wrsn.z_coords-wrsn.radius;
            [x,y,z]=meshgrid(X,Y,Z);
            w_points=[reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
            num=size(w_points,1);
            points=[];
            centers=w_points(randperm(num,wrsn.center_num),:);
            for i=1:wrsn.center_num
                lable=true;
                spos=centers(i,:);
                g_num=randi([20,30]);
                while lable
                    spoints=defineNodes(spos,g_num,wrsn.radius);
                    angle_mat=calcAngleMatrix(spoints-spos,g_num);
                    angle_mat(angle_mat<=wrsn.angle)=1;
                    angle_mat(angle_mat>wrsn.angle)=0;
                    lable=~all(sum(angle_mat,2)>1);
                end
                points=[points;spos;spoints];
                if size(points,1)>=wrsn.node_num
                    points=points(1:wrsn.node_num,:);
                    break;
                end
            end
        end
        function points = generateNodesBySparse(wrsn)
            xmin = 0; xmax = wrsn.x_coords;
            ymin = 0; ymax = wrsn.y_coords;
            zmin = 0; zmax = wrsn.z_coords;
            n = wrsn.node_num;
            min_dist = 10;
            points = zeros(n, 3);
            max_attempts_per_node = 1000;
            rng('shuffle');
            points(1, :) = [randi([xmin, xmax]), randi([ymin, ymax]), randi([zmin, zmax])];
            i = 2;
            while i <= n
                attempts = 0;
                while true
                    candidate_node = [randi([xmin, xmax]), randi([ymin, ymax]), randi([zmin, zmax])];
                    if ismember(candidate_node, points(1:i-1, :), 'rows')
                        continue;
                    end
                    distances = pdist2(candidate_node, points(1:i-1, :));
                    if all(distances >= min_dist)
                        points(i, :) = candidate_node;
                        i = i + 1;
                        break;
                    end
                    attempts = attempts + 1;
                    if attempts > max_attempts_per_node
                        error('无法在指定的尝试次数内找到一个满足最小距离 (%f) 的新节点。', min_dist);
                    end
                end
            end
        end
        %生成节点的能量
        function [E_org,E_exp,list_need_egy]=generateNodeInitEnergy(wrsn)
            E_org=wrsn.E_orgmin+round(wrsn.E_orgmax-wrsn.E_orgmin)*rand(wrsn.node_num,1);
            n_need_egy=round(wrsn.node_num*wrsn.n_need_ratio);
            list_need_egy=randperm(wrsn.node_num,n_need_egy);
            E_exp=E_org;
            E_exp(list_need_egy)=min(E_org(list_need_egy)+wrsn.egy_increase,wrsn.E_expmax);
        end
        function [t,exitflag,wrsn]=solveLP(wrsn,cMat,cone_lable)
                tic
                n_zone = length(cone_lable);
                if n_zone == 0
                    t = [];
                    exitflag = -1;
                    wrsn.lp_time = toc;
                    wrsn.charge_loss = 0;
                    wrsn.hover_loss = 0;
                    wrsn.Energy = -1;
                    return;
                end
                f_hover_rate = wrsn.hover_cost_rate * ones(n_zone, 1);
                f_loss_rate = wrsn.p * ones(n_zone, 1) - sum(cMat, 1)';
                f = f_hover_rate + f_loss_rate;
                A1 = -cMat;
                b1 = wrsn.E_org - wrsn.E_exp;
                A2 = ones(1, n_zone) * wrsn.p;
                b2 = wrsn.E_orgV - wrsn.E_orgL;
                A = [A1; A2];
                b = [b1; b2];
                lb = zeros(n_zone, 1);
                ub = [];
                options = optimset('display','off');
                [t,~,exitflag] = linprog(f,A,b,[],[],lb,ub,options);
                run_time=toc;
                wrsn.lp_time = run_time;
                if exitflag == 1
                    t(t < 1e-6) = 0;
                    wrsn.hover_loss = wrsn.hover_cost_rate * sum(t);
                    wrsn.charge_loss = f_loss_rate' * t;
                    wrsn.Energy = wrsn.init_total_energy - wrsn.hover_loss - wrsn.charge_loss;
                else
                    wrsn.charge_loss = -1;
                    wrsn.hover_loss = -1;
                    wrsn.Energy = -1;
                end
        end
        function [res,wrsn,exitflag]=directionCompareMethod(wrsn,iter)
            if isempty(wrsn.spots)
                wrsn=buildChargePosition(wrsn);
                wrsn=calcTransRateMatrix(wrsn);
                wrsn=buildDistanceMatrix(wrsn,wrsn.spots);
                [best_path,wrsn]=solveTSP(wrsn);
                distance_len=sum(wrsn.dist_mat(sub2ind(size(wrsn.dist_mat),best_path(1:end-1),best_path(2:end))))+wrsn.dist_mat(best_path(end),1);
                move_time=distance_len/wrsn.move_v;
                wrsn.move_t=move_time;
                wrsn.dist_len=distance_len;
            end
            [wrsn,cMat,cone_lable,~,~]=calcTransMatrix(wrsn);
            [t,exitflag,wrsn]=solveLP(wrsn,cMat,cone_lable);
            if exitflag~=1
                return
            end
            wrsn.charge_t=sum(t(t>0));
            wrsn.cost_t=wrsn.charge_t+wrsn.move_t;
            res=buildResult(wrsn);
        end
        function [res,wrsn]=atspCompareMethod(wrsn,iter)
             if isempty(wrsn.spots)
                wrsn=buildChargePosition(wrsn);
                wrsn=calcTransRateMatrix(wrsn);
               [wrsn,cMat,cone_lable,~,~]=calcTransMatrix(wrsn);
                [t,exitflag,wrsn]=solveLP(wrsn,cMat,cone_lable);
                if exitflag~=1
                    return
                end
                wrsn=buildDistanceMatrix(wrsn,wrsn.spots);
                wrsn.charge_t=sum(t(t>0));
             end
            [best_path,wrsn]=solveTSP(wrsn);
            wrsn.dist_len=sum(wrsn.dist_mat(sub2ind(size(wrsn.dist_mat),best_path(1:end-1),best_path(2:end))))+wrsn.dist_mat(best_path(end),1);
            wrsn.move_t=wrsn.dist_len/wrsn.move_v;
            res=buildResult(wrsn);
        end
        function wrsn=createModel(wrsn)
            buildTimes=1;
            pos_mode=["greedy","k-means","grid","node"];
            cone_mode=["cMRFDS","node","acc","gcc","fix"];
            tsp_mode=["lkh","greedy","aco"];
            while true
                wrsn=buildChargePosition(wrsn);
                wrsn=calcTransRateMatrix(wrsn);
                [wrsn,cMat,cone_lable,~,~]=calcTransMatrix(wrsn);
                [t,exitfalg,wrsn]=solveLP(wrsn,cMat,cone_lable);
                if exitfalg==1
                    spot_lable=t>0;
                    posCone=unique(cone_lable(spot_lable));
                    spots_pos=wrsn.spots(posCone,:);
                    wrsn.charge_t=sum(t(t>0));
                    wrsn=buildDistanceMatrix(wrsn,spots_pos);
                    break ;
                else
                    [org,exp,list_need]=generateNodeInitEnergy(wrsn);
                    wrsn.E_org=org;
                    wrsn.E_exp=exp;
                    wrsn.list_need_egy=list_need;
                    wrsn.points=generateNodes(wrsn);
                    buildTimes=buildTimes+1;
                    if buildTimes>100
                        break;
                    end
                end
            end
        end
        function [wrsn,cMat,cone_lable,spotConeInfo,coneInfo]=calcTransMatrix(wrsn)
            if wrsn.mode_d==0
                coneInfo={};
                [cMat,cone_lable,coneNumber,spotConeInfo]=calc_globe_Cmat(wrsn);
            else
                [cMat,cone_lable,coneNumber,spotConeInfo,coneInfo,time]=calcConeTransMatrix(wrsn);
            end
            wrsn.cone_number=coneNumber;
            wrsn.cdir_time=time;
        end
        function [cMat,zone_lable,cone_number,spot_cone,cone_info,total_time]=calcConeTransMatrix(wrsn)
            cMat=[];
            zone_lable=[];
            spot_cone={};
            cone_number=0;
            cone_info={};
            spots_n=size(wrsn.spots,1);
            total_time=0;
            for i=1:spots_n
                temp=[];
                info={};
                tag={};
                switch wrsn.mode_d
                    case 1
                        parameter.alpha1 = wrsn.alpha;
                        parameter.alpha2 = wrsn.alpha2;
                        parameter.beta1 = wrsn.beta;
                        parameter.beta2 = wrsn.beta2;
                        parameter.gamma1 = wrsn.gamma;
                        parameter.gamma2 = wrsn.gamma2;
                        [info,tag,time]=RuncMFEDS(wrsn.points, wrsn.spots(i,:), wrsn.angle, wrsn.angle2, wrsn.radius, wrsn.radius2, wrsn.angle_between,wrsn.cover_mat,i,parameter);
                    case 2
                        [info,tag,time]=DirNode_DualCone(wrsn.points,wrsn.radius,wrsn.radius2,wrsn.spots(i,:),wrsn.angle,wrsn.angle2,wrsn.angle_between,wrsn.cover_mat,i);
                    case 3
                        [info,tag,time]=DirACC_DualCone(wrsn.points,wrsn.radius,wrsn.radius2,wrsn.spots(i,:),wrsn.angle,wrsn.angle2,wrsn.angle_between,wrsn.cover_mat,i);
                    case 4
                        [info,tag,time]=DirGCC_DualCone_Sampled(wrsn.points,wrsn.radius,wrsn.radius2,wrsn.spots(i,:),wrsn.angle,wrsn.angle2,wrsn.angle_between,wrsn.cover_mat,i);
                    case 5
                        [info, tag, time] = DirPloyhedron_DualCone(wrsn.points, wrsn.radius, wrsn.radius2, wrsn.spots(i,:), wrsn.angle, wrsn.angle2, wrsn.angle_between, wrsn.cover_mat, i);
                end
                total_time=total_time+time;
                zone_num=numel(info);
                for j=1:zone_num
                    cone_info(end+1,1)=info(j);
                    temp=[temp,cone_number+j];
                    cZone=zeros(wrsn.node_num,1);
                    point_indices = cell2mat(info(j));
                    status_vector = cell2mat(tag(j));
                    if ~isempty(point_indices)
                        for k = 1:length(point_indices)
                            node_idx = point_indices(k);
                            node_status = status_vector(k);
                            cZone(node_idx) = wrsn.rate_mat(node_idx, i, node_status);
                        end
                        relevant_rates = cZone(point_indices);
                        sum_of_rates = sum(relevant_rates);
                        if sum_of_rates > 1
                            cZone(point_indices) = (relevant_rates ./ sum_of_rates) * wrsn.p * wrsn.c_max;
                        else
                            cZone(point_indices) = relevant_rates * wrsn.p;
                        end
                        cMat=[cMat,cZone];
                        zone_lable(end+1)=i;
                    end
                end
                cone_number=cone_number+zone_num;
                spot_cone{i,1}=temp;
            end
        end
        function [cMat,zone_lable,cone_number,spot_cone]=calc_globe_Cmat(wrsn)
            spots_n=size(wrsn.spots,1);
            for i=1:spots_n
                cMat(:,i)=wrsn.rate_mat(:,i)*wrsn.p/length(find(wrsn.rate_mat(:,i)~=0));
            end
            zone_lable=1:spots_n;
            cone_number=0;
            spot_cone=[];
        end
        function wrsn=calcTransRateMatrix(wrsn)
            num = size(wrsn.points, 1);
            spot_num = size(wrsn.spots, 1);
            rateMat3D = zeros(num, spot_num, 3);
            all_points = [wrsn.points; wrsn.spots];
            distance_mat = squareform(pdist(all_points));
            p2s = distance_mat(1:num, num+1:end);
            r = max(wrsn.radius, wrsn.radius2);
            cover_mask = p2s <= r;
            for i = 1:num
                for j = 1:spot_num
                    if cover_mask(i, j)
                        d = p2s(i, j);
                        rateMat3D(i, j, 1) = wrsn.alpha / (d + wrsn.beta)^wrsn.gamma;
                        rateMat3D(i, j, 2) = wrsn.alpha2 / (d + wrsn.beta2)^wrsn.gamma2;
                        rateMat3D(i, j, 3) = max(rateMat3D(i, j, 1), rateMat3D(i, j, 2));
                    end
                end
            end
            wrsn.rate_mat = rateMat3D;
        end
        function wrsn=buildDistanceMatrix(wrsn,spots_pos)
            all_points=[wrsn.pos;spots_pos];
            n=size(all_points,1);
            distance_mat=squareform(pdist(all_points));
            move_cost_mat=distance_mat*wrsn.move_cost_rate;
            wrsn.cost_mat=move_cost_mat;
            wrsn.dist_mat=distance_mat;
        end
        function res=calcModelLoss(wrsn,iter)
            tic
            switch wrsn.mode_p
                case {1,2,3,4,5}
                    wrsn=createModel(wrsn);
                    [best_path, wrsn] = solveTSP(wrsn);
                    if isempty(best_path) || ~isfinite(wrsn.path_loss)
                        res = struct();
                        return;
                    end
                    if numel(best_path) > 1
                        if isrow(best_path), best_path = best_path'; end
                        path_indices = sub2ind(size(wrsn.dist_mat), best_path, [best_path(2:end); best_path(1)]);
                        distance_len = sum(wrsn.dist_mat(path_indices));
                    else
                        distance_len = 0;
                    end
                    wrsn.dist_len = distance_len;
                    wrsn.move_t = wrsn.dist_len / wrsn.move_v;
                    wrsn.cost_t = wrsn.charge_t + wrsn.move_t;
                    wrsn.entire_time = toc;
                    res = buildResult(wrsn);
            end
        end

        %生成充电位置
        function wrsn=buildChargePosition(wrsn)
            switch wrsn.mode_p
                case 1
                    wrsn=buildPositionByGreedy(wrsn);
                case 2
                    wrsn=buildPositionByKmeans(wrsn);
                case 3
                    wrsn=buildPositionByGrid(wrsn);
                case 4
                    wrsn=buildPositionByNodes(wrsn);
                case 5
                    wrsn=buildPositionConDir(wrsn);
                case 6
                    [wrsn.spots, wrsn.cover_mat,~,wrsn.points]=definePosition(wrsn,wrsn.node_num,wrsn.radius,20);
            end
        end

        %TSP问题求解
        function [bestPath, wrsn] = solveTSP(wrsn)
            bestPath = [];
            pathLoss = inf;
            tic
            switch wrsn.mode_tsp
            case 1
            [bestPath, pathLoss] = lkhTSP(wrsn.cost_mat);
            case 2
            [bestPath, pathLoss] = greedyTSP(wrsn.cost_mat);
            case 3
            [bestPath, pathLoss] = acoTSP(wrsn.cost_mat);
            end
            run_time = toc;
            wrsn.path_loss = pathLoss;
            wrsn.tsp_time = run_time;
        end


        %返回结果构造
        function [res]=buildResult(wrsn)
                res.hover_loss=wrsn.hover_loss;
                res.charge_loss=wrsn.charge_loss;
                res.path_loss=wrsn.path_loss;
                res.spots=size(wrsn.spots,1);
                res.cones=wrsn.cone_number;
                res.span=wrsn.cost_t;
                res.move_t=wrsn.move_t;
                res.charge_t=wrsn.charge_t;
                res.dist_len=wrsn.dist_len;
                res.atsp_time=wrsn.tsp_time;
                res.dir_time=wrsn.cdir_time;
                res.lp_time=wrsn.lp_time;
                res.entrie_time=wrsn.entire_time;
        end
        % 贪心的方式生成充电位置
        function  wrsn=buildPositionByGreedy(wrsn)
                nodeNum=wrsn.node_num;
                r = max(wrsn.radius, wrsn.radius2);
                nodes=wrsn.points;
                pdist_mat=squareform(pdist(nodes));
                neighborMat=pdist_mat<=r;
                [~,index]=sort(sum(neighborMat,1),"descend");
                firstChargePosIndex=index(1);
                chargePos=nodes(firstChargePosIndex,:);
                row=1;
                while true
                    all_pos=[nodes;chargePos];
                    distanceMat=squareform(pdist(all_pos));
                    p2s=distanceMat(1:nodeNum,nodeNum+1:end);
                    coverMat=(p2s<=r);
                    lable=all(sum(coverMat,2)>=1);
                    if lable
                        break
                    else
                        not_cover=find(~sum(coverMat,2)>=1);
                        [~,bs_index]=sort(sum(neighborMat(not_cover,:),2),"descend");
                        c_index=not_cover(bs_index(1));
                        chargePos=[chargePos;nodes(c_index,:)];
                    end
                end
                wrsn.spots=chargePos;
                wrsn.cover_mat=coverMat;
        end
function  wrsn=buildPositionConDir(wrsn)
        nodeNum=wrsn.node_num;
        r = max(wrsn.radius, wrsn.radius2);
        nodes=wrsn.points;
        pdist_mat=squareform(pdist(nodes));
        neighborMat=pdist_mat<=r;
        [~,index]=sort(sum(neighborMat,1),"descend");
        firstChargePosIndex=index(1);
        chargePos=nodes(firstChargePosIndex,:);
        row=1;
        while true
            all_pos=[nodes;chargePos];
            distanceMat=squareform(pdist(all_pos));
            p2s=distanceMat(1:nodeNum,nodeNum+1:end);
            coverMat=(p2s<=r);
            lable=all(sum(coverMat,2)>=1);
            if lable
                break
            else
                not_cover=find(~sum(coverMat,2)>=1);
                [~,bs_index]=sort(sum(neighborMat(not_cover,:),2),"descend");
                c_index=not_cover(bs_index(1));
                chargePos=[chargePos;nodes(c_index,:)];
            end
        end
        wrsn.spots=chargePos;
        wrsn.cover_mat=coverMat;
end
    % k-means++生成充电位置
    function wrsn=buildPositionByKmeans(wrsn)
        nodeNum=wrsn.node_num;
        r = max(wrsn.radius, wrsn.radius2);
        nodes=wrsn.points;
        chargePosNum=1;
        spots_index=randperm(nodeNum,chargePosNum);
        chargePos=nodes(spots_index,:);
        adjustTimes=1;
        loop=1;
        while true
            allPos=[nodes;chargePos];
            distanceMatrix=squareform(pdist(allPos));
            p2s=distanceMatrix(1:nodeNum,nodeNum+1:end);
            coverMatrix=(p2s<=r);
            if loop>nodeNum
                break;
            end
            fullCoverFlag=all(sum(coverMatrix,2)>=1);
            if fullCoverFlag
                if adjustTimes>20
                    break;
                end
                    for i=1:chargePosNum
                        spoints=nodes(coverMatrix(:,i),:);
                        if size(spoints,1)==1
                            continue ;
                        end
                        chargePos(i,:)=mean(spoints,1);
                    end
                    allPos=[nodes;chargePos];
                    distanceMatrix=squareform(pdist(allPos));
                    p2s=distanceMatrix(1:nodeNum,nodeNum+1:end);
                    coverMatrix=(p2s<=r);
                    b_cover=coverMatrix;
                    coverLable=true(1,chargePosNum);
                    for i=1:chargePosNum
                        i_cover=b_cover(:,i);
                        b_cover(:,i)=0;
                        c_cover=b_cover-i_cover;
                        if any(~any(c_cover<0,1))
                            coverLable(i)=0;
                        end
                    end
                    chargePos=chargePos(coverLable,:);
                    coverMatrix=coverMatrix(:,coverLable);
                    chargePosNum=size(chargePos,1);
                if ~all(sum(coverMatrix,2)>=1)
                    adjustTimes=adjustTimes+1;
                else
                    break;
                end
            else
                loop=loop+1;
                notCover=find(~sum(coverMatrix,2)>=1);
                n2s=min(p2s(notCover,:),[],2);
                probabilities=n2s./sum(n2s);
                n=length(probabilities);
                cumulativePro = 0;
                selected_object = 1;
                randomNumber = rand();
                for i = 1:n
                    cumulativePro = cumulativePro + probabilities(i);
                    if randomNumber <= cumulativePro
                        selected_object = i;
                        break;
                    end
                end
                nextSpot=notCover(selected_object);
                chargePos=[chargePos;nodes(nextSpot,:)];
                chargePosNum=size(chargePos,1);
            end
        end
        wrsn.spots=chargePos;
        wrsn.cover_mat=coverMatrix;
    end
    %网格化的方式生成充电位置
    function wrsn=buildPositionByGrid(wrsn)
        nodeNum=wrsn.node_num;
        r = max(wrsn.radius, wrsn.radius2);
        nodes=wrsn.points;
        xMax=wrsn.x_coords;
        yMax=wrsn.y_coords;
        zMax=wrsn.z_coords;
        d=5;
        X=d/2:d:xMax+d;
        Y=d/2:d:yMax+d;
        Z=d/2:d:zMax+d;
        [x,y,z]=meshgrid(X,Y,Z);
        b_spots=[reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
        allPos=[nodes;b_spots];
        distanceMatrix=squareform(pdist(allPos));
        p2s=distanceMatrix(1:nodeNum,nodeNum+1:end);
        coverMatrix=(p2s<=r);
        coverLable=sum(coverMatrix,1)>0;
        chargePos=b_spots(coverLable,:);
        coverMatrix=coverMatrix(:,coverLable);
        chargePosNum=size(chargePos,1);
        b_cover=coverMatrix;
        cover_lable=true(1,chargePosNum);
        for i=1:chargePosNum
            i_cover=b_cover(:,i);
            b_cover(:,i)=0;
            c_cover=b_cover-i_cover;
            if any(~any(c_cover<0,1))
                cover_lable(i)=0;
            end
        end
        coverMatrix=coverMatrix(:,cover_lable);
        chargePos=chargePos(cover_lable,:);
        wrsn.spots=chargePos;
        wrsn.cover_mat=coverMatrix;
    end
    function wrsn=buildPositionByNodes(wrsn)
        nodes=wrsn.points;
        r = max(wrsn.radius, wrsn.radius2);
        chargePos=nodes;
        nodeNum=size(nodes,1);
        allPos=[nodes;chargePos];
        distanceMatrix=squareform(pdist(allPos));
        p2s=distanceMatrix(1:nodeNum,nodeNum+1:end);
        coverMatrix=(p2s<=r);
        wrsn.spots=chargePos;
        wrsn.cover_mat=coverMatrix;
    end
end
end
function [spots,cover_mat,spots_n,points]=definePosition(wrsn,node_num,radius,g_num)
    d=sqrt(4*radius^2)/3;
    X=d/2:d:wrsn.x_coords+d;
    Y=d/2:d:wrsn.y_coords+d;
    Z=d/2:d:wrsn.z_coords+d;
    [x,y,z]=meshgrid(X,Y,Z);
    b_spots=[reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
    spots_n=node_num/g_num;
    list_spots=randperm(size(b_spots,1),spots_n);
    spots=b_spots(list_spots,:);
    points=[];
    for i=1:spots_n
        lable=true;
        spos=spots(i,:);
        spoints=[];
        while lable
            spoints=defineNodes(spos,g_num,radius);
            angle_mat=calcAngleMatrix(spoints-spos,g_num);
            angle_mat(angle_mat<=wrsn.angle)=1;
            angle_mat(angle_mat>wrsn.angle)=0;
            lable=~all(sum(angle_mat,2)>1);
        end
        points=[points;spoints];
    end
    all_pos=[points;spots];
    dist_mat=squareform(pdist(all_pos));
    p2s=dist_mat(1:node_num,node_num+1:end);
    r = max(wrsn.radius, wrsn.radius2);
    cover_mat=(p2s<=r);
end
%% 计算角度矩阵
function point_angle=calcAngleMatrix(pos_to_points,g_num)
    point_angle=zeros(g_num,g_num);
    for i=1:g_num
        for j=1:g_num
            if i==j
                point_angle(i,j)=0;
                continue;
            end
            dot_tmp=dot(pos_to_points(i,:),pos_to_points(j,:));
            magnitude_i = norm(pos_to_points(i,:));
            magnitude_j = norm(pos_to_points(j,:));
            angleRad = acos(dot_tmp / (magnitude_i * magnitude_j));
            point_angle(i,j) = rad2deg(angleRad);
        end
    end
end
function spoints=defineNodes(spot,num,radius)
    spoints=zeros(num,3);
    theta=rand(num,1)*pi;
    phi=rand(num,1)*2*pi;
    d=rand(num,1)*radius;
    spoints(:,1)=spot(1)+d.*sin(theta).*cos(phi);
    spoints(:,2)=spot(2)+d.*sin(theta).*sin(phi);
    spoints(:,3)=spot(3)+d.*cos(theta);
end
