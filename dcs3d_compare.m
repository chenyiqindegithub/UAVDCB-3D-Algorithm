function dcs3d_compare()
    angle=60;
    angle2=90;
    angle_between=60;
    nums=300:50:800;
    rate_list=0.1:0.1:1;
    iternum=4;
    baseDir="./data/sharetest/";
    rlist=6;
    radius2=6;
    nodeNum=400;
    % 显示开始时间
    startTime = datetime('now');
    disp(['程序开始时间: ', char(startTime)]);
%     saveExptsp(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum);
%     saveExpDir(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum);
%     saveExpPos(baseDir,angle,angle2,rlist,radius2,nums,iternum,angle_between);
%     saveExpHovFlyRate(baseDir,angle,angle2,rlist,radius2,angle_between,rate_list,nodeNum,iternum);
    saveExpEntire(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum);
    % 显示结束时间
    endTime = datetime('now');
    disp(['程序结束时间: ', char(endTime)]);
end
%对比不同的无人机参数
%对比不同的无人机参数
function saveExpHovFlyRate(baseDir,angle,angle2,rlist,radius2,angle_between,rate_list,n,iternum)
     rateDir=strcat(baseDir,"rate/");
     mkdir(rateDir);
     len=length(rate_list); 
     exp_list=[1,2,3,4];
     exp_num=length(exp_list);%实验组数
     for radius=rlist
        res=cell(exp_num,len,iternum);
        lable=1;
        for i=1:len
            rate=rate_list(i);
            for j=1:iternum
                wrsn=wrsn_config(n,radius,radius2,angle,angle2,angle_between);
                p=wrsn.move_cost_rate+wrsn.hover_cost_rate;
                switch lable
                     case 1   
                        wrsn.hover_cost_rate=p*(rate/(rate+1));
                        wrsn.move_cost_rate=p*(1/(rate+1));
                     case 2
                        wrsn.move_cost_rate=p*(rate/(rate+1));
                        wrsn.hover_cost_rate=p*(1/(rate+1));
                end
                for e=1:exp_num
                    wrsn.mode_p=exp_list(e);
                    res{e,i,j}=calcModelLoss(wrsn,j);
                end
            end 
            disp(['参数:',num2str(rate),'完成！！！']);
        end
         filename=strcat(rateDir,"node_r_",num2str(radius),".mat");
        save(filename,"res");
    end
end
%对比不同的方向选取算法
%对比不同的方向选取算法
function saveExpDir(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum)
    coneDir=strcat(baseDir,"cone/");
    mkdir(coneDir);
    len=length(nums); 
    exp_list=[1,2,3,4,5];
    exp_num=length(exp_list);%实验组数
    % 显示开始时间
    startTime = datetime('now');
    disp(['对比实验开始: ', char(startTime)]);
    for radius=rlist
        % 显示开始时间
        startTime = datetime('now');
        disp('参数');
        disp(radius);
        disp(['开始:', char(startTime)]);
        res=cell(exp_num,len,iternum);
        for i=1:len
            n=nums(i);
            % 显示开始时间
            startTime = datetime('now');
            disp('节点数 ');
            disp(n);
            disp(['开始迭代:', char(startTime)]);
            for j=1:iternum
                wrsn=wrsn_config(n,radius,radius2,angle,angle2,angle_between);
                %位置选取 k_means
                wrsn.mode_p=2;
                wrsn.isDense=1;
%                 wrsn.isDense=0;
                for e=1:exp_num
                   %方向选取
                   wrsn.mode_d=exp_list(e);
%                    [res{e,i,j},wrsn,~]=directionCompareMethod(wrsn,j);

                   [tmp_res, wrsn, exitflag] = directionCompareMethod(wrsn,j);
                   if exitflag == 1
                       res{e,i,j} = tmp_res;
                   else
                       res{e,i,j} = [];  % 防止未赋值
                       warning('LP failed at e=%d, i=%d, j=%d', e, i, j);
                   end

                end
            end
             % 显示结束时间
            endTime = datetime('now');
            disp('节点数 ');
            disp(nums(i));
            disp(['结束迭代:', char(endTime)]);
            disp(['参数:',num2str(n),'完成！！！']);
        end
        filename=strcat(coneDir,"node_r_",num2str(radius),".mat");
        save(filename,"res");
        % 显示结束时间
        endTime = datetime('now');
        disp('参数 ');
        disp(radius);
        disp(['结束:', char(endTime)]);
    end
    % 显示结束时间
    endTime = datetime('now');
    disp(['对比实验结束: ', char(endTime)]);
end
%对比不同TSP求解算法
function saveExptsp(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum)
     atspDir=strcat(baseDir,"tsp/");
     mkdir(atspDir);
     len=length(nums); 
     exp_list=[1,2,3];
     exp_num=length(exp_list);%实验组数
    for radius=rlist
        res=cell(exp_num,len,iternum);
        for i=1:len
            n=nums(i);
            for j=1:iternum
                wrsn=wrsn_config(n,radius,radius2,angle,angle2,angle_between);
                %位置选取 k_means
                wrsn.mode_p=2;
                wrsn.isDense=0;
                for e=1:exp_num
                   %TSP选取
                   wrsn.mode_tsp=exp_list(e);
                   [res{e,i,j},wrsn]=atspCompareMethod(wrsn,j);
                end
            end
            disp(['参数:',num2str(n),'完成！！！']);
        end
        filename=strcat(atspDir,"node_r_",num2str(radius),".mat");
        save(filename,"res");
    end
end
%对比不同的位置选取方式
function saveExpPos(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum)
     posDir=strcat(baseDir,"pos/");
     mkdir(posDir);
     len=length(nums); 
%      exp_list=[1,2,3,4,5];
     exp_list=[1,2,3,4];
     exp_num=length(exp_list);%实验组数
    for radius=rlist
        res=cell(exp_num,len,iternum);
        for i=1:len
            n=nums(i);
            for j=1:iternum
                wrsn=wrsn_config(n,radius,radius2,angle,angle2,angle_between);
                for e=1:exp_num
                    wrsn.mode_p=exp_list(e);
                    res{e,i,j}=calcModelLoss(wrsn,j);
                end
            end 
            disp(['参数:',num2str(n),'完成！！！']);
        end
     filename=strcat(posDir,"node_r_",num2str(radius),".mat");
     save(filename,"res");
    end
end
%对比整体算法

function saveExpEntire(baseDir,angle,angle2,rlist,radius2,nums,angle_between,iternum)
     rateDir=strcat(baseDir,"entire/");
     mkdir(rateDir);
     len=length(nums); 
     exp_num=3;%实验组数
    for radius=rlist
        res=cell(exp_num,len,iternum);
        for i=1:len
            n=nums(i);
            for j=1:iternum
                wrsn=wrsn_config(n,radius,radius2,angle,angle2,angle_between);
                wrsn.mode_tsp=3;
                wrsn.mode_p=1;
                wrsn.mode_d=5;%贪心（Group），固定方向 aco
                res{1,i,j}=calcModelLoss(wrsn,j);
%                 wrsn.mode_tsp=3;
%                 wrsn.mode_p=2;
%                 wrsn.mode_d=5;%k_means，固定方向 aco
%                 res{2,i,j}=calcModelLoss(wrsn,j);
%                 wrsn.mode_tsp=2;
%                 wrsn.mode_p=3;
%                 wrsn.mode_d=4;%网格，GCC a2sym
%                 res{3,i,j}=calcModelLoss(wrsn,j);
                wrsn.mode_tsp=2;
                wrsn.mode_p=3;
                wrsn.mode_d=3;%网格，ACC greedy
                res{2,i,j}=calcModelLoss(wrsn,j);
                wrsn.mode_tsp=1;
                wrsn.mode_p=4;
                wrsn.mode_d=1;%节点，cMFEDS lkh
                res{3,i,j}=calcModelLoss(wrsn,j);
            end
            disp(['参数:',num2str(n),'完成！！！']);
        end
    end
     filename=strcat(rateDir,"node_r_",num2str(radius),".mat");
     save(filename,"res");
end