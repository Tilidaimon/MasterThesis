function [AssignTarget, MaxGlobalUtility, GlobalUtility] = SAP_lagrange_function(TargetValue,RequireNum,Time_to_go,J_opt,Epis)

[AgentNum,TargetNum] = size(Time_to_go);


ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = zeros(1,AgentNum);%分配任务表
% 初始任务分配
% AssignTarget = randperm(TargetNum);

for i = 1:AgentNum
    ProbDist(i,:) = ProbDist(i,:)/sum(ProbDist(i,:));
    e = rand;
    cp = 0;
    for j = 1:TargetNum
        cp = cp + ProbDist(i,j);
        if e <= cp
            AssignTarget(i) = j;
            break
        end
    end
end
AssignMatrix = zeros(AgentNum,TargetNum);
for i = 1:AgentNum
    AssignMatrix(i,AssignTarget(i)) = 1;
end

%% 迭代更新
T = Epis;%总迭代步数
GlobalUtility = zeros(1,T);
lambda = 10;


for j = 1:TargetNum
        GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue(j), AssignTarget, RequireNum, Time_to_go, J_opt, lambda);
end

for k = 2:T
    %for i = 1:AgentNum
        i = randi(AgentNum);
        % 假设第i个智能体不选任务
        Assign2 = AssignTarget;
        Assign2(i) = 0;
        
        tau = 10/k^2;
        % 计算k+1时刻不同任务的效用值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,TargetValue(j),Assign1,RequireNum,Time_to_go,J_opt,lambda)...
                - TargetUtility(j,TargetValue(j),Assign2,RequireNum,Time_to_go,J_opt,lambda);% 计算第i个智能体的边际收益
        end
        %计算概率分布，此处是以熵最大化为原则选择概率分布
        AU = AgentUtility/tau;
        newProb = exp(AU)/sum(exp(AU));
        flag = 0;
        for j = 1:TargetNum
            if isnan(newProb(j))
                flag = 1;
                break;
            end
        end
        if flag == 0
            ProbDist(i,:) = newProb;
        else
            AU = AU/max(abs(AU))*k/2;
            ProbDist(i,:) = exp(AU)/sum(exp(AU));
        end
        
        % 任务分配
        e = rand;
        cp = 0;
        for j = 1:TargetNum
            cp = cp + ProbDist(i,j);
            if e <= cp
                AssignTarget(i) = j;
                AssignMatrix(i,j) = 1;
                break
            end
        end
    %end
    lagrange_multipler = 0;
    for p=1:TargetNum
        part_missiles = find(AssignTarget==p);
        num_missiles = length(part_missiles);
        lagrange_multipler = lagrange_multipler + min(RequireNum(p)-num_missiles,0);
    end
    %计算总效用
    for j = 1:AgentNum
        s = AssignTarget(j);
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(s, TargetValue(s), AssignTarget, RequireNum, Time_to_go, J_opt, lambda);
        GlobalUtility(k) = GlobalUtility(k) - lagrange_multipler;
    end
end

MaxGlobalUtility = GlobalUtility(Epis);%-lambda*(sum(RequireNum)-AgentNum);

end

%计算效用函数
function [U] = TargetUtility(target,target_value,AssignTarget,require_num,Time_to_go,J_opt,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = require_num(target);
if num_missiles == 0
    task_cost = 0;
    U = 0;
else
    Time_to_go_max = max(Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(J_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = target_value - task_cost + con;
end



end