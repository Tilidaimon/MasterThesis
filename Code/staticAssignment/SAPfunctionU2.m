function [AssignTarget, MaxGlobalUtility, GlobalUtility] = SAPfunctionU2(AgentNum, TargetNum, TargetValue, Distance, CompleteProb)



ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randperm(TargetNum);%分配任务表

%% 迭代更新
T = 500;%总迭代步数
GlobalUtility = zeros(1,T);


for j = 1:TargetNum
    GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, Distance, AssignTarget, CompleteProb);
end

for k = 1:T
    i = randi(AgentNum);
    % 假设第i个智能体不选任务
    Assign2 = AssignTarget;
        % 假设第i个智能体不选任务
        tau = 10/k^2;
        % 计算k+1时刻不同任务的效用值
        AgentUdist = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        for j = 1:TargetNum
            Assign1 = Assign2;
            exindex = find(Assign2==j);
            extarget = Assign1(i);
            Assign1(i) = j;
            Assign1(exindex) = extarget;
            AgentU1 = AgentUtility(i,TargetValue,Distance,Assign1,CompleteProb)-AgentUtility(i,TargetValue,Distance,Assign2,CompleteProb);
            AgentU2 = AgentUtility(exindex,TargetValue,Distance,Assign1,CompleteProb)-AgentUtility(exindex,TargetValue,Distance,Assign2,CompleteProb);
            AgentUdist(j) = AgentU1+AgentU2;% 计算第i个智能体的边际收益
        end
        %计算概率分布，此处是以熵最大化为原则选择概率分布
        AU = AgentUdist/tau;
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
            extarget = AssignTarget(i);
            AssignTarget(AssignTarget==j) = extarget;
            AssignTarget(i) = j;
            break
        end
    end
    
    %计算总效用
    for s = 1:TargetNum
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(s, TargetValue, Distance, AssignTarget, CompleteProb);
    end
end

MaxGlobalUtility = GlobalUtility(T);

end

%计算效用函数
function [U] = TargetUtility(target,targetValue,distance,AssignTarget,completeProb)

Agents = find(AssignTarget==target);
nAgent = length(Agents);
totalProb = 1;
cost = 0;
if nAgent>0
    for i = 1:nAgent
        totalProb = totalProb * (1-completeProb(Agents(i),target));
        cost = cost + distance(i,target);
    end
else
    totalProb = 0;
end
% lambda = 0.1;
U = (targetValue(target)-cost)*(1-totalProb);%+lambda*(1-nAgent);

end

function [U] = AgentUtility(agent,targetValue,distance,AssignTarget,completeProb)

Assign = AssignTarget;
Assign(agent) = 0;
target = AssignTarget(agent);
U = TargetUtility(target,targetValue,distance,AssignTarget,completeProb)...
    -TargetUtility(target,targetValue,distance,Assign,completeProb);
end