% 建立导弹与目标攻击动态模型
% 毕设第五章 Chapter 5
% 随机模型 Stochastic Game Model

classdef DynamicMissileAndTarget
    % 导弹与目标模型 Model of Missiles and Target
    % 需要指定的属性 properties
    properties
        num_missiles; % 导弹数量 numbers of missiles
        num_targets; % 目标数量 numbers of targets
        target_require_num_list; % 目标所需导弹数量 required numbers of missiles to attack targets
        order_targets;
        
        % 智能体结构体，包含位置、速度、加速度、航向角
        Missiles;
        Targets;
        
    end
    
    % 常量
    properties
        dT = 0.02; % 采样时间0.02s
        tau_m = 0.04; % 导弹等效时间常数
        tau_t = 0.06; % 目标等效时间常数
            
        maxVOfMissile  = 1.5; % 导弹最大速度
        maxVOfTarget = 0.3; % 目标最大速度
        GOfMissile = 0.5; % 导弹最大加速度G值
        GOfTarget = 0.03; % 目标加速度G值
        
        max_inexcapable_radius = 15; % 最大不可逃逸区半径
        communication_radius = 30; % 通信半径
        radar_radius = 60; % 弹载雷达探测半径
    end
    
    % 计算得到的常量
    properties
        Tmax = 10000; % 参考最大攻击时间，这个时间由第一次攻击得到
        Time_to_go =[]; % 导弹预计剩余攻击时间
        extend_time_togo = []; % 扩展目标预计剩余攻击时间
        Energy_opt = []; % 最优控制律能量
        extend_energy = []; % 扩展目标预计能量
        extend_value = [];
        Adjacent = []; % 邻接矩阵
        distance_missiles_targets = []; % 导弹与目标距离矩阵
        oig_acc = []; % 最优控制制导律矩阵
        reassign_flag = 0; % 是否立即重分配
        maxtime = 0;
    end
    
    methods
        %% 构造函数
        function model = DynamicMissileAndTarget(num_missiles, num_targets)
            model.num_missiles = num_missiles;
            model.num_targets = num_targets;
        end
        
        %% 随机设置态势参数
        function model = RandomSetSituation(model)
            model.Missiles.p = 10 * rand(model.num_missiles, 2);
%             model.Missiles.p(model.num_missiles/2+1:model.num_missiles,2) = 40 + 10*rand(model.num_missiles/2, 1);
            %model.Missiles.p(:,2) = 20 * rand(model.num_missiles, 1);
            model.Missiles.v = 1 * ones(model.num_missiles, 1);
            model.Missiles.angle = 0.1 * pi * (rand(model.num_missiles, 1)-0.5);
            model.Missiles.alive = true(model.num_missiles,1);
            model.Missiles.alive_missiles = 1:model.num_missiles;
            model.Missiles.alive_num = model.num_missiles;
            model.Missiles.diff_R = -1*ones(model.num_missiles,1);
            
            
            model.Targets.p = 20*rand(model.num_targets, 2);
            model.Targets.p(:,1) = 30 + 5*rand(model.num_targets, 1);
            model.Targets.v = 0.6*ones(model.num_targets, 1);
            model.Targets.angle = 10/180*pi*rand(model.num_targets, 1);
            model.Targets.acc = zeros(model.num_targets, 1);
            model.Targets.alive = true(model.num_targets,1);
            model.Targets.num_attacked = zeros(model.num_targets,1);
            model.Targets.alive_targets = 1:model.num_targets;
            model.Targets.alive_num = model.num_targets;
            
            % 需要用两枚导弹攻击的目标序号
            index_important_targets = randperm(model.num_targets, model.num_missiles-model.num_targets);
            model.target_require_num_list = ones(model.num_targets,1);
            model.target_require_num_list(index_important_targets) = model.target_require_num_list(index_important_targets) +1;
            
            % 虚拟目标增广序列
            k=1;
            model.order_targets = zeros(model.num_targets,1);
            for i=1:model.num_targets
                for j=1:model.target_require_num_list(i)
                    model.order_targets(k) = i;
                    k = k+1;
                end
            end
            
            % 目标价值
            model.Targets.value = 200*ones(model.num_targets,1);
            model.Targets.value(index_important_targets) = 2*model.Targets.value(index_important_targets);
            model.extend_value = model.Targets.value(model.order_targets);
            
            % 计算矩阵邻接矩阵
            model.Missiles.distance_missiles = DistanceMatrixMissiles(model);
            model.Adjacent = GetAdjacentMatrix(model);
            
            % 计算导弹与目标距离矩阵
            model.distance_missiles_targets = DistanceMatrixMissileAndTarget(model);
            model.Missiles.target_set = DetectTargetSet(model);
            
            % 共享目标信息
            model = TargetShared(model);
            
            % 计算剩余攻击时间矩阵
            model.Time_to_go = AttackTimeToGo(model);
            model.extend_time_togo = model.Time_to_go(:,model.order_targets);
            model.maxtime = max(max(model.Time_to_go));
            
            % 计算最优制导律及预计消耗能量
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model,0);
            model.extend_energy = model.Energy_opt(:,model.order_targets);
            
            
        end
        
        %% 用户指定态势参数
        function model = SetSituation(model, Missiles, Targets, index_important_targets, target_value)
            model.Missiles = Missiles;
            model.Targets = Targets;
            
            model.target_require_num_list = ones(model.num_targets,1);
            model.target_require_num_list(index_important_targets) = model.target_require_num_list(index_important_targets) +1;
            
            % 目标价值
            model.Targets.value = target_value;
            
            % 计算矩阵邻接矩阵
            model.Missiles.distance_missiles = DistanceMatrixMissiles(model);
            model.Adjacent = GetAdjacentMatrix(model);
            
            % 计算导弹与目标距离矩阵
            model.distance_missiles_targets = DistanceMatrixMissileAndTarget(model);
            model.Missiles.target_set = DetectTargetSet(model);
            
            % 共享目标信息
            model = TargetShared(model);
            
            % 计算剩余攻击时间矩阵
            model.Time_to_go = AttackTimeToGo(model);
            
            % 计算最优制导律及预计消耗能量
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model);
            model.maxtime = max(max(model.Time_to_go));
        end
        
        %% 根据可选目标集获得初始解
        function initial_plan = GetInitialPlan(model)
            initial_plan = zeros(model.Missiles.alive_num,1);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                select_targets = find(model.Missiles.target_set(m,:)==1);
                select_index = randi(length(select_targets));
                initial_plan(i) = select_targets(select_index);
            end
        end
        
        %% 获得尚未完成攻击的导弹集合
        function model = GetAliveMissiles(model)
            index = 1:model.num_missiles;
            model.Missiles.alive_missiles = index(model.Missiles.alive==1);
            model.Missiles.alive_num = sum(model.Missiles.alive);
        end
        
        %% 获得尚未被击落的目标集合
        function model = GetAliveTarget(model)
            index = 1:model.num_targets;
            model.Targets.alive_targets = index(model.Targets.alive==1);
            model.Targets.alive_num = sum(model.Targets.alive);
        end
        
        %% 取出生存导弹和目标的分配子集
        function sub_plan = GetAliveAssignPlan(model,assign_plan)
            sub_plan = zeros(model.Missiles.alive_num,1);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                if t == 0
                    sub_plan(i) = 0;
                else
                    sub_plan(i) = find(model.Targets.alive_targets==t);
                end
            end
        end
        
        %% 由分配子集反推得到全体分配
        function new_assign_plan = GetAllAssignPlan(model,assign_plan,sub_plan)
            new_assign_plan = assign_plan;
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                new_assign_plan(m) = model.Targets.alive_targets(sub_plan(i));
            end
        end
        
        %% 计算导弹邻接矩阵
        function adjacent = GetAdjacentMatrix(model)
            n = model.Missiles.alive_num;
            adjacent = false(n);
            for i=1:n
                for j=1:n
                    if model.Missiles.distance_missiles(i,j) < model.communication_radius
                        adjacent(i,j) = true;
                    end
                end
            end
        end
        
        %% 计算导弹探测到的目标集合
        function target_set = DetectTargetSet(model)
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            target_set = false(Nm, Nt);
            for i=1:Nm
                for j=1:Nt
                    if model.distance_missiles_targets(i,j) < model.radar_radius
                        target_set(i,j) = true;
                    end
                end
            end
        end
        
        %% 导弹共享目标集合 Missiles in same group share targets
        function model = TargetShared(model)
            %target_list = false(model.num_missiles,model.num_targets);
            for i=1:model.Missiles.alive_num
                for j=1:model.Missiles.alive_num
                    if model.Adjacent(i,j) == true
                        model.Missiles.target_set(j,:) = model.Missiles.target_set(j,:) | model.Missiles.target_set(i,:);
                    end
                end
                %target_list(model.Adjacent(i,:)==true,:) = ones(num_neighbors,1) * target;
            end
        end
        
        %% 计算剩余攻击时间矩阵
        function TimeToGo = AttackTimeToGo(model)
            N = 3; % 比例系数
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            TimeToGo = zeros(Nm, Nt);
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                for j=1:Nt
                    t = model.Targets.alive_targets(j);
                    d = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(t))]; % 导弹方向向量
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
%                     relative_theta = (model.Missiles.v(m)*sin(alpha) - model.Targets.v(t)*sin(beta))/d;
                    relative_v = model.Missiles.v(m) - model.Targets.v(t)*cos(beta);
                    TimeToGo(i,j) = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                end
            end
        end
        
        %% 计算最优制导律及其预计消耗能量
        function [acc, energy] = OptimalInterceptGuidance(model,step)
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            acc = zeros(Nm, Nt);
            energy = zeros(Nm, Nt);
            t_go = model.maxtime - step*model.dT;
            T_m = t_go/model.tau_m;
            T_t = t_go/model.tau_t;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                for j=1:Nt
                    t = model.Targets.alive_targets(j);
                    
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
                    R = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(t))]; % 导弹方向向量
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
%                     dotR = model.Targets.v(t)*cos(beta) - model.Missiles.v(m)*cos(alpha); % 弹目距离变化率
                    
                    
                    % 导弹前一时刻位置
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % 目标前一时刻位置
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % 前一时刻弹目距离
                    dotR = (R-pre_R)/model.dT; % 弹目距离变化率
                    
                    acc(i,j) = Noig/t_go^2 * (R + dotR*t_go + model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i) - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(j));
%                     acc(i,j) = sign(acc(i,j)) * min(abs(acc(i,j)),model.GOfMissile);
                    
                    P = [1;t_go;model.tau_m^2*(exp(-T_m)-1+T_m)];
                    X_state = [R;dotR;model.Missiles.acc(i)];
                    energy(i,j) = 0.5*X_state'*(P*P')*X_state;
                end
            end
        end
        
        
        
        %% 导弹根据最优制导律OIG移动
        function [model,assign_plan] = MissilesMoveByOIG(model,assign_plan,step)
            N=3;
            Nm = model.Missiles.alive_num;
            time_to_go = model.maxtime - step*model.dT;
            T_m = time_to_go/model.tau_m;
            T_t = time_to_go/model.tau_t;
%             Nt = model.Targets.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                %target = model.Targets.alive_targets(t);
                if t ~= 0 % 导弹被分配到目标
%                     d = model.distance_missiles_targets(i,model.Targets.alive_targets(t));
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))]; % 导弹方向向量
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    R = norm(b);
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
%                     aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % 目标方向向量
%                     beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
%                     relative_v = model.Missiles.v(m) - model.Targets.v(t)*cos(beta);
                    

                    % 导弹前一时刻位置
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % 目标前一时刻位置
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % 前一时刻弹目距离
                    diff_R = (R-pre_R)/model.dT; % 弹目距离变化率
                    
%                     time_to_go = R/(-diff_R)*(1+0.5*alpha^2/(2*N-1));
                    
                    
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
%                     TimeToGo = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                    
                    %                     dotR = model.Targets.v(target)*cos(beta) - model.Missiles.v(i)*cos(alpha); % 弹目距离变化率
                    acc = Noig/time_to_go^2 * (R + diff_R*time_to_go - model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i)...
                        - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(t));
                    model.Missiles.acc(i) = sign(acc)*min(model.GOfMissile,abs(acc));
                    if R<0.05 || diff_R>0% 击落距离30m
                        model.Targets.num_attacked(t) = model.Targets.num_attacked(t) + 1;
                        model.Missiles.alive(m) = 0;
                        model.Missiles.v(m) = 0; % 击中目标的导弹停止移动
                        model.Missiles.target_set(i,:) = false;
                        
                        
                        if model.Targets.num_attacked(t) == model.target_require_num_list(t)
                            
                            model.Targets.v(t) = 0; % 目标被击中，停止移动
                            model.Targets.alive(t) = 0;
                            model.Missiles.target_set(:,t) = false; % 被击中的目标不再成为可选目标
                            if sum(assign_plan==t)>1
                                co_missiles = assign_plan == t;
                                co_m = co_missiles & model.Missiles.alive;
                                assign_plan(co_m==1) = 0;
                                model.reassign_flag = 1;
                            end
                            
                                
                        end
                    else
                        
                        accMX = model.Missiles.acc(m)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(m)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                        model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                        model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                        model.Missiles.p(m,:) = model.Missiles.p(m,:) + model.Missiles.v(m)*model.dT*...
                            [cos(model.Missiles.angle(m)), sin(model.Missiles.angle(m))];
%                         model.Missiles.angle(m) = r;
                    end
                end
            end
        end
        
        %% 导弹以PNG方式移动
        function [model,assign_plan] = MissilesMoveByPNG(model,assign_plan)
            N=4;
            Nm = model.Missiles.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                %t_i = find(model.Targets.alive_targets==t);
                if t ~= 0 % 导弹被分配到目标
                    %R = model.distance_missiles_targets(i,t_i);
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    R = norm(b);
                    v = [-model.Missiles.v(m)*cos(model.Missiles.angle(m)) + model.Targets.v(t)*cos(model.Targets.angle(t)),...
                        -model.Missiles.v(m)*sin(model.Missiles.angle(m)) + model.Targets.v(t)*sin(model.Targets.angle(t))];
                    
                    % 导弹前一时刻位置
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % 目标前一时刻位置
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % 前一时刻弹目距离
                    diff_R = (R-pre_R)/model.dT; % 弹目距离变化率
                    pre_diff_R = model.Missiles.diff_R(m);
                    model.Missiles.diff_R(m) = diff_R;
%                   
%                     lambda = atan2(b(2),b(1));
%                     pre_lambda = atan2(pre_target_p(2)-pre_missiles_p(2),pre_target_p(1)-pre_missiles_p(1));
%                     Qtm = (lambda-pre_lambda)/model.dT;
                    
                    
                    Vtm = -(b(1)*v(1)+b(2)*v(2))/R;
                    Qtm = (b(1)*v(2)-b(2)*v(1))/(R^2);
                    if Qtm<-pi
                        Qtm = Qtm + 2*pi;
                    elseif Qtm>pi
                        Qtm = Qtm - 2*pi;
                    end
                    acc = N*Vtm*Qtm;
                    model.Missiles.acc(m) = sign(acc)*min(abs(acc),model.GOfMissile);
%                     r = model.Missiles.angle(m) + model.Missiles.acc(m)*model.dT;
                    if R<0.05 || diff_R>0% 击落距离30m
                        model.Targets.num_attacked(t) = model.Targets.num_attacked(t) + 1;
                        model.Missiles.alive(m) = 0;
                        model.Missiles.v(m) = 0; % 击中目标的导弹停止移动
                        model.Missiles.target_set(i,:) = false;
                        
                        
                        if model.Targets.num_attacked(t) == model.target_require_num_list(t)
                            
                            model.Targets.v(t) = 0; % 目标被击中，停止移动
                            model.Targets.alive(t) = 0;
                            model.Missiles.target_set(:,t) = false; % 被击中的目标不再成为可选目标
                            if sum(assign_plan==t)>1
                                co_missiles = assign_plan == t;
                                co_m = co_missiles & model.Missiles.alive;
                                assign_plan(co_m==1) = 0;
                                model.reassign_flag = 1;
                            end
                            
                                
                        end
                    else
                        
                        accMX = model.Missiles.acc(m)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(m)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                        model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                        model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                        model.Missiles.p(m,:) = model.Missiles.p(m,:) + model.Missiles.v(m)*model.dT*...
                            [cos(model.Missiles.angle(m)), sin(model.Missiles.angle(m))];
%                         model.Missiles.angle(m) = r;
                    end
                end
               

            end
        end

        
        %% 判断导弹是否可以改变目标
        function flag = MissileCanChangeTarget(model)
            flag = false(model.Missiles.alive_num,model.Targets.alive_num);
            N=4;
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            for i=1:Nm
                for j=1:Nt
                    m = model.Missiles.alive_missiles(i);
                    t = model.Targets.alive_targets(j);
                    %t_i = find(model.Targets.alive_targets==t);
                    if t ~= 0 % 导弹被分配到目标
                        %R = model.distance_missiles_targets(i,t_i);
                        b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                        R = norm(b);
                        v = [-model.Missiles.v(m)*cos(model.Missiles.angle(m)) + model.Targets.v(t)*cos(model.Targets.angle(t)),...
                            -model.Missiles.v(m)*sin(model.Missiles.angle(m)) + model.Targets.v(t)*sin(model.Targets.angle(t))];
                        
                        % 导弹前一时刻位置
                        pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                        % 目标前一时刻位置
                        pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
                        %
                        pre_R = norm(pre_missiles_p - pre_target_p); % 前一时刻弹目距离
                        Vtm = (R-pre_R)/model.dT; % 弹目距离变化率
                        %pre_diff_R = model.Missiles.diff_R(m);
                        
                        lambda = atan2(b(2),b(1));
                        pre_lambda = atan2(pre_target_p(2)-pre_missiles_p(2),pre_target_p(1)-pre_missiles_p(1));
                        Qtm = (lambda-pre_lambda)/model.dT;
                        
                        acc = N*Vtm*Qtm;
                        
                        if abs(acc) <= model.GOfMissile
                            flag(i,j) = true & model.Missiles.target_set(i,j);
                        end
                    end
                    
                    %                     Vtm = -(b(1)*v(1)+b(2)*v(2))/R;
                    %                     Qtm = (b(1)*v(2)-b(2)*v(1))/(R^2);
                    
                end
            end
        end
        
        %% 目标移动
        function model = TargetMove(model)
            model.Targets.p = model.Targets.p + ...
                model.dT*[model.Targets.v .* cos(model.Targets.angle), model.Targets.v .* sin(model.Targets.angle)];
        end
        
%         %% 生成每次迭代的分配初始解
%         function initial_plan = GetInitialAssignment(model,pre_plan)
%             % pre_plan为前一时刻的分配解，第一次为随机生成
%             for i=1:model.num_missiles
%                 
%         end
        
        %% 计算导弹目标距离矩阵 Get the distance matrix
        function dismatrix = DistanceMatrixMissileAndTarget(model)
            dismatrix = zeros(model.Missiles.alive_num, model.Targets.alive_num);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                for j=1:model.Targets.alive_num
                    t = model.Targets.alive_targets(j);
                    dismatrix(i,j) = norm(model.Missiles.p(m,:)-model.Targets.p(t,:));
                end
            end
            
        end
        
         %% 计算导弹目标距离矩阵 Get the distance matrix
        function dismatrix = DistanceMatrixMissiles(model)
            dismatrix = zeros(model.Missiles.alive_num, model.Targets.alive_num);
            for i=1:model.Missiles.alive_num
                m1 = model.Missiles.alive_missiles(i);
                for j=1:model.Missiles.alive_num
                    m2 = model.Missiles.alive_missiles(j);
                    if i==j
                        dismatrix(i,i) = 0;
                    elseif i<j
                        dismatrix(i,j) = norm(model.Missiles.p(m1,:)-model.Missiles.p(m2,:));
                    else
                        dismatrix(i,j) = dismatrix(j,i);
                    end
                end
            end
            
        end
        
        %% 解码
        function plan = DecodeFromExtendPlan(model,extend_plan)
            plan = zeros(model.num_missiles,1);
            for i=1:model.num_missiles
                plan(i) = model.order_targets(extend_plan(i));
            end
        end
        
    end
end



