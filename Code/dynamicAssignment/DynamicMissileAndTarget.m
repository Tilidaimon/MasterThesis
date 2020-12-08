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
        
        % 智能体结构体，包含位置、速度、加速度、航向角
        Missiles;
        Targets;
        
    end
    
    % 常量
    properties
        dT = 0.02; % 采样时间0.02s
        tau_m = 0.4; % 导弹等效时间常数
        tau_t = 0.6; % 目标等效时间常数
            
        maxVOfMissile  = 1.5; % 导弹最大速度
        maxVOfTarget = 0.3; % 目标最大速度
        GOfMissile = 0.06 * 9.8; % 导弹最大加速度G值
        GOfTarget = 0.03 * 9.8; % 目标加速度G值
        
        max_inexcapable_radius = 25; % 最大不可逃逸区半径
        communication_radius = 30; % 通信半径
        radar_radius = 50; % 弹载雷达探测半径
    end
    
    % 计算得到的常量
    properties
        Tmax = 10000; % 参考最大攻击时间，这个时间由第一次攻击得到
        Time_to_go =[]; % 导弹预计剩余攻击时间
        Energy_opt = []; % 最优控制律能量
        Adjacent = []; % 邻接矩阵
        distance_missiles_targets = []; % 导弹与目标距离矩阵
        oig_acc = []; % 最优控制制导律矩阵
    end
    
    methods
        %% 构造函数
        function model = DynamicMissileAndTarget(num_missiles, num_targets)
            model.num_missiles = num_missiles;
            model.num_targets = num_targets;
        end
        
        %% 随机设置态势参数
        function model = RandomSetSituation(model)
            model.Missiles.p = 20 * rand(model.num_missiles, 2);
            model.Missiles.p(:,2) = 30 * rand(model.num_missiles, 1);
            model.Missiles.v = model.maxVOfMissile * ones(model.num_missiles, 1);
            model.Missiles.angle = 0.1 * pi * (rand(model.num_missiles, 1)-0.5);
            model.Missiles.alive = true(model.num_missiles,1);
            model.Missiles.alive_missiles = 1:model.num_missiles;
            model.Missiles.alive_num = model.num_missiles;
            
            
            model.Targets.p = 20*rand(model.num_targets, 2);
            model.Targets.p(:,1) = 60 + 5*rand(model.num_targets, 1);
            model.Targets.v = 0.3*ones(model.num_targets, 1);
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
            
            % 目标价值
            model.Targets.value = 200*ones(model.num_targets,1);
            model.Targets.value(index_important_targets) = 2*model.Targets.value(index_important_targets);
            
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
                sub_plan(i) = find(model.Targets.alive_targets==t);
            end
        end
        
        %% 由分配子集反推得到全体分配
        function new_assign_plan = GetAllAssignPlan(model,assign_plan,sub_plan)
            new_assign_plan = assign_plan;
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                t_i = sub_plan(i);
                assign_plan(m) = model.Targets.alive_targets(t_i);
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
                    relative_v = model.Missiles.v(m) - model.Targets.v(t)*cos(beta);
                    TimeToGo(i,j) = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                end
            end
        end
        
        %% 计算最优制导律及其预计消耗能量
        function [acc, energy] = OptimalInterceptGuidance(model)
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            acc = zeros(Nm, Nt);
            energy = zeros(Nm, Nt);
            
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                for j=1:Nt
                    t = model.Targets.alive_targets(j);
                    t_go = model.Time_to_go(i,j);
                    T_m = t_go/model.tau_m;
                    T_t = t_go/model.tau_t;
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
                    R = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(t))]; % 导弹方向向量
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
                    dotR = model.Targets.v(t)*cos(beta) - model.Missiles.v(m)*cos(alpha); % 弹目距离变化率
                    acc(i,j) = Noig/t_go^2 * (R + dotR*t_go + model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i)) - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(j);
                    acc(i,j) = sign(acc(i,j)) * min(abs(acc(i,j)),model.GOfMissile);
                    
                    P = [1;t_go;model.tau_m^2*(exp(-T_m)-1+T_m)];
                    X_state = [R;dotR;model.Missiles.acc(i)];
                    energy(i,j) = 0.5*X_state'*(P*P')*X_state;
                end
            end
        end
        
        
        
        %% 导弹根据最优制导律OIG移动
        function model = MissilesMoveByOIG(model,assign_plan)
            N=3;
            Nm = model.Missiles.alive_num;
%             Nt = model.Targets.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                target = model.Targets.alive_targets(t);
                if target ~= 0 % 导弹被分配到目标
                    d = model.distance_missiles_targets(i,t);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))]; % 导弹方向向量
                    b = [model.Targets.p(target,1)-model.Missiles.p(m,1), model.Targets.p(target,2)-model.Missiles.p(m,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(target)),sin(model.Targets.angle(target))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
                    relative_v = model.Missiles.v(m) - model.Targets.v(target)*cos(beta);
                    TimeToGo = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                    
                    R=d;
                    T_m = TimeToGo/model.tau_m;
                    T_t = TimeToGo/model.tau_t;
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
                    dotR = -relative_v;
                    %                     dotR = model.Targets.v(target)*cos(beta) - model.Missiles.v(i)*cos(alpha); % 弹目距离变化率
                    model.Missiles.acc(i) = Noig/TimeToGo^2 * (R + dotR*TimeToGo - model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i))...
                        - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(target);
                    model.Missiles.acc(i) = sign(model.Missiles.acc(i))*min(model.GOfMissile,abs(model.Missiles.acc(i)));
                    if model.distance_missiles_targets(i,t)<0.03 % 击落距离30m
                        model.Targets.num_attacked(target) = model.Targets.num_attacked(target) + 1;
                        model.Missiles.alive(i) = 0;
                        model.Missiles.v(m) = 0; % 击中目标的导弹停止移动
                        model.Missiles.target_set(i,:) = 0;
                        
                        if model.Targets.num_attacked(target) == model.target_require_num_list(target)
                            model.Targets.v(target) = 0; % 目标被击中，停止移动
                            model.Targets.alive(target) = false;
                            model.Missiles.target_set(:,target) = false; % 被击中的目标不再成为可选目标
                            
                        end
                    else
                        accMX = model.Missiles.acc(i)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(i)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                    end
                end
                model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                model.Missiles.p(i,:) = model.Missiles.p(i,:) + model.Missiles.v(i)*model.dT*...
                    [cos(model.Missiles.angle(i)), sin(model.Missiles.angle(i))];
            end
        end
        
        %% 导弹以PNG方式移动
        function model = MissilesMoveByPNG(model,assign_plan)
            N=3;
            Nm = model.Missiles.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                t_i = find(model.Targets.alive_targets==t);
                if t ~= 0 % 导弹被分配到目标
                    R = model.distance_missiles_targets(i,t_i);
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % 视线角向量
                    v = [-model.Missiles.v(m)*cos(model.Missiles.angle(m)) + model.Targets.v(t)*cos(model.Targets.angle(t)),...
                        -model.Missiles.v(m)*sin(model.Missiles.angle(m)) + model.Targets.v(t)*sin(model.Targets.angle(t))];
                    Vtm = -(b(1)*v(1)+b(2)*v(2))/R;
                    Rtm = (b(1)*v(2)-b(2)*v(1))/(R^2);
                    model.Missiles.acc(m) = N*Vtm*Rtm;
                    if model.distance_missiles_targets(i,t_i)<0.03 % 击落距离30m
                        model.Targets.num_attacked(t) = model.Targets.num_attacked(t) + 1;
                        model.Missiles.alive(m) = 0;
                        model.Missiles.v(m) = 0; % 击中目标的导弹停止移动
                        model.Missiles.target_set(i,:) = false;
                        
                        
                        if model.Targets.num_attacked(t) == model.target_require_num_list(t)
                            model.Targets.v(t) = 0; % 目标被击中，停止移动
                            model.Targets.alive(t) = 0;
                            model.Missiles.target_set(:,t) = false; % 被击中的目标不再成为可选目标
                        end
                    else
                        accMX = model.Missiles.acc(m)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(m)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                    end
                    model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                    model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                    model.Missiles.p(m,:) = model.Missiles.p(m,:) + model.Missiles.v(m)*model.dT*...
                        [cos(model.Missiles.angle(m)), sin(model.Missiles.angle(m))];
                end
            end
        end

        
        %% 判断导弹是否可以改变目标
        function flag = MissileCanChangeTarget(model)
            acc_flag = ~(abs(model.oig_acc) > model.GOfMissile);
            target_flag = acc_flag & model.Missiles.target_set;
            flag = logical(sum(target_flag,2));
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
        
    end
end



