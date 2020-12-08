clear
clc
file = fopen('log.txt','w');

% 初始化模型
num_missiles = 12;
num_targets = 10;
model = DynamicMissileAndTarget(num_missiles,num_targets);
model = RandomSetSituation(model);

stepMissiles = 6000; % 最大迭代步数
stepTime = 1500; % 最大分配算法迭代步数

% 数据存储
missiles_save.p(1,:,:) = model.Missiles.p;
targets_save.p(1,:,:) = model.Targets.p;

missiles_save.angle(1,:) = model.Missiles.angle;
targets_save.angle(1,:) = model.Targets.angle;
initial_plan = randi(model.num_targets,model.num_missiles,1);
assign_plan = SAPLagrangeFunction(model,stepTime,initial_plan);
fprintf(file,"missileplan0=");
fprintf(file,num2str(assign_plan'));
fprintf(file,"\n");

%% 动态过程
stopSteps = stepMissiles;
model.Targets.acc =0*(rand(model.num_targets,1)-0.5);  % 目标随机加速度引起的角度变化
for i=2:stepMissiles
    pre_assign_plan = assign_plan; % 前一时刻分配
    pre_missiles_alive = model.Missiles.alive;
    pre_targets_attack = model.Targets.num_attacked;
    if rem(i,50)==0 && model.Missiles.alive_num>1
        % 计算剩余攻击时间矩阵
        model.Time_to_go = AttackTimeToGo(model);
        
        % 计算最优制导律及预计消耗能量
        [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model);
        
        % 进行分配
        sub_plan = GetAliveAssignPlan(model,assign_plan);
        new_sub_plan = SAPLagrangeFunction(model,stepTime,sub_plan);
        new_assign_plan = GetAllAssignPlan(model,assign_plan,new_sub_plan);
        sub_plan = new_sub_plan;
        if sum(new_assign_plan ~= pre_assign_plan) > 0
            flag = MissileCanChangeTarget(model);
            for k=1:model.num_missiles
                if flag(k)==0
                    assign_plan(k) = pre_assign_plan(k);
                else
                    assign_plan(k) = new_assign_plan(k);
                end
                
            end
            if sum(flag)>0
                fprintf(file,"Change");fprintf(file,"\n");
                fprintf(file,"currentPlan=");
                fprintf(file,num2str(assign_plan'));fprintf(file,"\n");
            end
        end
    end
    
    % 统计还存活的导弹和目标
    model = GetAliveMissiles(model);
    model = GetAliveTarget(model);
    
    % 目标移动
    model.Targets.angle = model.Targets.angle + model.Targets.acc*model.dT;
    
    % 计算矩阵邻接矩阵
    model.Missiles.distance_missiles = DistanceMatrixMissiles(model);
    model.Adjacent = GetAdjacentMatrix(model);
    
    % 计算导弹与目标距离矩阵
    model.distance_missiles_targets = DistanceMatrixMissileAndTarget(model);
    model.Missiles.target_set = DetectTargetSet(model);
    
    % 共享目标信息
    model = TargetShared(model);
    
    model = MissilesMoveByPNG(model,assign_plan);
%     model.Missiles.acc
    model = TargetMove(model);
    
    
    
    if sum(pre_targets_attack ~= model.Targets.num_attacked)>0
        for k=1:model.num_targets
            if pre_targets_attack(k) ~= model.Targets.num_attacked(k)
                fprintf(file,"Target ");
                fprintf(file,num2str(k));
                fprintf(file," has been hit ");
                fprintf(file,num2str(model.Targets.num_attacked(k)));
                fprintf(file," of ");
                fprintf(file,num2str(model.target_require_num_list(k)));
                fprintf(file," times by Missile");
                missile_indexs = find(assign_plan==k);
                if length(missile_indexs)>1
                    for p=1:length(missile_indexs)
                        m = missile_indexs(p);
                        if model.Missiles.alive(m) == 0 && pre_missiles_alive(m) == 1
                            attack_index = p;
                        end
                    end
                else
                    attack_index = missile_indexs;
                end
                fprintf(file,num2str(attack_index));
                fprintf(file," \n");
            end
        end
    end
    
    missiles_save.p(i,:,:) = model.Missiles.p;
    targets_save.p(i,:,:) = model.Targets.p;
    
    missiles_save.angle(i,:) = model.Missiles.angle;
    targets_save.angle(i,:) = model.Targets.angle;
    
    % 检测是否停止攻击
    stopFlag = logical(sum(model.Targets.alive));
    if ~stopFlag
        stopsteps = i;
        break;
    end
end


%% 绘图 Plot
figure(1);


% 绘制导弹位置
theta = 0:0.1:2*pi-0.1;
r1 = model.communication_radius;
r2 = model.radar_radius;
x = missiles_save.p(1,1,1);
y = missiles_save.p(1,1,2);
scatter(x,y,'b'); hold on;
quiver(x,y, cos(missiles_save.angle(1,1))',sin(missiles_save.angle(1,1))',5,'color','b');
text(x+1,y+1,num2str(1));
% plot(x+r1*cos(theta),y+r1*sin(theta),'g--');
% plot(x+r2*cos(theta),y+r2*sin(theta),'b--');

for i=1:model.num_missiles
    x = missiles_save.p(1,i,1);
    y = missiles_save.p(1,i,2);
    scatter(x,y,'b');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(missiles_save.angle(1,i))',sin(missiles_save.angle(1,i))',5,'color','b');
%     plot(x+r1*cos(theta),y+r1*sin(theta),'g--');
%     plot(x+r2*cos(theta),y+r2*sin(theta),'b--');
end

% for i=1:model.num_missiles-1
%     for j=(i+1):num_missiles
%         if model.Adjacent(i,j) == true
%             plot([missiles_save.p(1,i,1),missiles_save.p(1,j,1)],[missiles_save.p(1,i,2),missiles_save.p(1,j,2)],'g');
%         end
%     end
% end

% 绘制目标位置
for i=1:model.num_targets
    x = targets_save.p(1,i,1);
    y = targets_save.p(1,i,2);
    scatter(x,y,'r');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(targets_save.angle(1,i))',sin(targets_save.angle(1,i))',5,'color','r');
end

% for i=1:model.num_missiles
%     plot([missiles_save.p(1,i,1),targets_save.p(1,assign_plan(i),1)],...
%         [missiles_save.p(1,i,2),targets_save.p(1,assign_plan(i),2)],'m');
% end

for i=1:model.num_missiles
    plot(missiles_save.p(:,i,1),missiles_save.p(:,i,2),'b');
    plot(targets_save.p(:,assign_plan(i),1),targets_save.p(:,assign_plan(i),2),'r');
end