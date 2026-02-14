% 蚁群解决TSP问题
function  [best_path,best_length]=acoTSP(cost_mat)
    ant_n=50;
    alpha=1;
    beta=5;
    rho=0.1;
    Q=1;
    iter_max=300;
    
    spot_n=size(cost_mat,1);
    
    if spot_n <= 1
        if spot_n == 1
            best_path = 1;
        else
            best_path = [];
        end
        best_length = 0;
        return;
    end

    ants=zeros(ant_n,spot_n);
    eta=10./cost_mat;
    eta(isinf(eta)) = 0;
    Tau=zeros(spot_n,spot_n);
    
    best_length=Inf;
    best_path=[];
    all_result=[];
    
    for iter=1:iter_max
        paths=zeros(1,ant_n);
        for ant=1:ant_n
            start=randi([2,spot_n]);
            visited=false(1,spot_n);
            current=start;
            visited(1)=true;
            visited(current)=true;
            ants(ant,1)=1;
            ants(ant,2)=start;
            row=3;
            
            while ~all(visited)
                unvisitedSpots = find(~visited);
                
                prob_vec = Tau(current,unvisitedSpots).^alpha .* eta(current,unvisitedSpots).^beta;
                
                if sum(prob_vec) == 0
                    prob_vec = eta(current,unvisitedSpots).^beta;
                end
                
                if sum(prob_vec) == 0
                    prob_vec = ones(1, length(unvisitedSpots));
                end
                
                probabilities = prob_vec / sum(prob_vec);

                r = rand();
                c_p = cumsum(probabilities);
                selected_object = find(r <= c_p, 1, 'first');
                
                if isempty(selected_object)
                    selected_object = 1;
                end

                nextSpot=unvisitedSpots(selected_object);
                current=nextSpot;
                ants(ant,row)=current;
                visited(current)=true;
                row=row+1;
            end
            
            current_path = ants(ant,:);
            len = 0;
            for k = 1:spot_n-1
                len = len + cost_mat(current_path(k), current_path(k+1));
            end
            len = len + cost_mat(current_path(end), current_path(1));
            paths(ant) = len;
        end
        
        [temp_length,index]=min(paths);
        
        if temp_length < best_length
            best_length=temp_length;
            best_path=ants(index,:);
        end
        
        delta_tau=zeros(spot_n,spot_n);
        for i=1:ant_n
            for j=1:spot_n-1
                a=ants(i,j);
                b=ants(i,j+1);
                delta_tau(a,b)=delta_tau(a,b)+Q/paths(i);
                delta_tau(b,a)=delta_tau(b,a)+Q/paths(i);
            end
            a=ants(i,end);
            b=ants(i,1);
            delta_tau(a,b)=delta_tau(a,b)+Q/paths(i);
            delta_tau(b,a)=delta_tau(b,a)+Q/paths(i);
        end
        Tau=(1-rho).*Tau+delta_tau;
        
        all_result(iter)=best_length;
    end
end
