% 建立导弹与目标攻击动态模型
% 毕设第五章 Chapter 5
% 随机模型 Stochastic Game Model

classdef MissileAndTarget
    % 导弹与目标模型 Model of Missiles and Target
    % 需要指定的属性 properties
    properties
        num_missiles; % 导弹数量 numbers of missiles
        num_targets; % 目标数量 numbers of targets
        target_require_num_list; % 目标所需导弹数量 required numbers of missiles to attack targets
        assign_plan; % 分配结果
        
        % 智能体结构体，包含位置、速度、加速度、航向角
        Missiles;
        Targets;
    end
    
    % 常量
    properties
        dT = 0.02; % 采样时间0.02s
        maxVOfMissile  = 1.5; % 导弹最大速度
        maxVOfTarget = 0.3; % 目标最大速度
        GOfMissile = 0.06 * 9.8; % 导弹的加速度G值
        GOfTarget = 0.03 * 9.8; % 导弹加速度G值
        
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
        function model = MissileAndTarget(num_missiles, num_targets)
            model.num_missiles = num_missiles;
            model.num_targets = num_targets;
        end
        
        %% 随机设置态势参数
        function model = RandomSet(model)
            model.Missiles.p = 100 * rand(model.num_missiles, 2);
            model.Missiles.v = 0.4 * rand(model.num_missiles, 1);
            model.Missiles.angle = 0.1 * pi * rand(model.num_missiles, 1);
            
            
            model.Targets.p = 20*rand(model.num_targets, 2);
            model.Targets.p(:,1) = 60 + 5*rand(model.num_targets, 1);
            model.Targets.v = 0.3*ones(model.num_targets, 1);
            model.Targets.angle = 10/180*pi*rand(model.num_targets, 1);
            
            index_important_targets = randperm(model.)
            model.Targets.value = 
            
            % 计算矩阵邻接矩阵
            model.Missiles.distance_missiles = DistanceMatrix(model.Missiles.p, model.Missiles.p);
            model.Adjacent = false(model.num_missiles);
            for i=1:model.num_missiles
                for j=1:model.num_missiles
                    if model.Missiles.distance_missiles(i,j) < model.communication_radius
                        model.Adjacent(i,j) = true;
                    end
                end
            end
            
            % 计算导弹与目标距离矩阵
            model.distance_missiles_targets = DistanceMatrix(model.Missiles.p, model.Targets.p);
            model.Missiles.target_set = false(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    if model.distance_missiles_targets(i,j) < model.radar_radius
                        model.Missiles.target_set(i,j) = true;
                    end
                end
            end
            
            % 共享目标信息
            model.Missiles.target_set = TargetShared(model.Missiles.target_set, model.Adjacent);
            
            % 计算剩余攻击时间矩阵
            model.Time_to_go = AttackTimeToGo(model);
            
            % 计算最优制导律及预计消耗能量
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model);
            
        end
        
        %% 导弹共享目标集合 Missiles in same group share targets
        function target_list = TargetShared(model)
            target_list = false(model.num_missiles);
            for i=1:numMissiles
                target = model.Missiles.target_set(i,:);
                num_neighbors = 0;
                for j=1:numMissiles
                    if model.Adjacent(i,j) == true
                        target = target | model.Missiles.target_set(j,:);
                        num_neighbors = num_neighbors+1;
                    end
                end
                target_list(model.Adjacent(i,:)==true,:) = ones(num_neighbors,1) * target;
            end
        end
        
        %% 计算剩余攻击时间矩阵
        function TimeToGo = AttackTimeToGo(model)
            N = 4; % 比例系数
            TimeToGo = zeros(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    d = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(i)),sin(model.Missiles.angle(i))]; % 导弹方向向量
                    b = [model.Targets.p(j,1)-model.Missiles.p(i,1), model.Targets.p(j,2)-model.Missiles.p(i,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(j)),sin(model.Targets.angle(j))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
                    relative_v = model.Missiles.v(i) - model.Targets.v(j)*cos(beta);
                    model.Time_to_go(i,j) = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                end
            end
        end
        
        %% 计算最优制导律及其预计消耗能量
        function [acc, energy] = OptimalInterceptGuidance(model)
            tau = 0.4; % 导弹等效时间常数
            acc = zeros(model.num_missiles, model.num_targets);
            energy = zeros(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    t_go = model.Time_to_go(i,j);
                    T = t_go/tau;
                    Noig = T^2*(exp(-T)-1+T)/(-0.5*exp(-2*T)-2*T*exp(-T)+T^3/3-T^2+T+0.5);
                    R = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(i)),sin(model.Missiles.angle(i))]; % 导弹方向向量
                    b = [model.Targets.p(j,1)-model.Missiles.p(i,1), model.Targets.p(j,2)-model.Missiles.p(i,2)]; % 视线角向量
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
                    aT = [cos(model.Targets.angle(j)),sin(model.Targets.angle(j))]; % 目标方向向量
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
                    dotR = model.Targets.v(j)*cos(beta) - model.Missiles.v(i)*cos(alpha); % 弹目距离变化率
                    acc(i,j) = Noig/t_go^2 * (R+dotR*t_go - tau^2*(exp(-T)-1+T)*model.Missiles.acc(i)-t_go^2*Target.acc(j)/2);
                    
                    P = [1;t_go;tau^2*(exp(-T)-1+T)];
                    X_state = [R;dotR;A_c(i,j)];
                    energy(i,j) = 0.5*X_state'*(P*P')*X_state;
                end
            end
        end
        
        
        
    end
end

%% 计算距离矩阵 Get the distance matrix
function [dismatrix] = DistanceMatrix(pos1, pos2)

[numMissiles,~] = size(pos1);
[numTargets,~] = size(pos2);

dismatrix = zeros(numMissiles, numTargets);
for i=1:pos1
    for j=1:pos2
        dismatrix(i,j) = norm(pos1(i,:)-pos2(j,:));
    end
end

end

