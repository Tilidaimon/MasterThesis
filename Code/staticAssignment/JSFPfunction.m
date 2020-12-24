function [AssignTarget, Ug, time] = JSFPfunction(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;

AssignTarget = randi(TargetNum,AgentNum,1);%分配任务表
%% 迭代更新
T = Epis;%总迭代步数
% Ug = zeros(1,T);
MeanUtility = zeros(AgentNum,TargetNum);%后悔矩阵
Ug = zeros(T,1);
Ug(1) = GlobalUtility(model,AssignTarget);
% 在第k步迭代时，计算每个智能体的后悔值，即前k-1步没有选择任务l的后悔程度，后悔程度以任务l带来的效用提高衡量

satisfied = zeros(AgentNum,1);
lambda = max(model.Targets.value);
tic
for k=1:T
    newAssignTarget = AssignTarget;
    for i = 1:AgentNum
        % 假设第i个智能体不选任务
        Assign2 = AssignTarget;
        Assign2(i) = 0;
        
        % 计算k+1时刻不同任务的效用值和后悔值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        %Regret = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的后悔值
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,model,Assign1,lambda)...
                - TargetUtility(j,model,Assign2,lambda);% 计算第i个智能体的边际收益
            %计算k+1时刻的后悔值
            MeanUtility(i,j) = ((k-1)*MeanUtility(i,j) + 1*AgentUtility(j))/(k);
        end
        % 计算智能体i在k+1时刻执行不同任务的概率分布
        %RM = Regret/sum(Regret);% 后悔匹配函数
        
        
        % 任务分配
        select_target = find(MeanUtility(i,:)==max(MeanUtility(i,:)));
        prob = rand;
        acp = 0.5;
        if length(select_target)>1
            select_index = randi(length(select_target));
            newassign = select_target(select_index);
        else
            newassign = select_target;
        end
        if prob>acp
            newAssignTarget(i) = newassign;
        end
        if newAssignTarget(i)==AssignTarget(i)
            satisfied(i) = satisfied(i) + 1;
        else
            satisfied(i) = 0;
        end
    end
    AssignTarget = newAssignTarget;
    %计算总效用
    Ug(k) = GlobalUtility(model,AssignTarget);
    if sum(satisfied>20)==AgentNum
        time = toc
        Ug(k+1:T) = Ug(k)*ones(T-k,1);
        break
    end
end
time = toc

end

%% 计算效用函数
function [U] = TargetUtility(target,model,AssignTarget,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = model.target_require_num_list(target);
if num_missiles == 0
%     con = lambda*(num_require-num_missiles);
%     task_cost = 0;
    U = 0;
else
    Time_to_go_max = sum(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = model.Targets.value(target)*num_missiles - task_cost + con;
end


end

%% 全局效用
function [Ug] = GlobalUtility(model,plan)
Nt = model.num_targets;
Ut = zeros(Nt,1);
flag = 1;
for j=1:Nt
    part_missiles = find(plan == j);
    num_missiles = length(part_missiles);
    if num_missiles > model.target_require_num_list(j)
        flag = 0;
        break
    end
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = sum(model.Time_to_go(part_missiles,j));
        J_sum = sum(model.Energy_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,model.Targets.value(j)*num_missiles/model.target_require_num_list(j) - cost);
    end
    
end
if flag == 0
    Ug = 0;
else
    Ug = sum(Ut);
end
end