function [AssignTarget, MaxGlobalUtility, GlobalUtility] = SAPfuncSimple(AgentNum,TargetNum,TargetValue,CompleteProb)


ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randperm(TargetNum);%分配任务表
AssignTarget = AssignTarget(1:AgentNum);


%% 迭代更新
T = 500;%总迭代步数
GlobalUtility = zeros(1,T);


for j = 1:TargetNum
    GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, AssignTarget, CompleteProb);
end

for k = 2:T
%     selectAgentNum = 100;
%     reorderAgent = randperm(AgentNum);
%     selectAgent = reorderAgent(1:selectAgentNum);
%     for p = 1:selectAgentNum
%     for i = 1:AgentNum
        i = randi(AgentNum);
%         i = selectAgent(p);
       if i==1
            Assign2 = [0,AssignTarget(2:AgentNum)];
        elseif i == AgentNum
            Assign2 = [AssignTarget(1:AgentNum-1),0];
        else
            Assign2 = [AssignTarget(1:i-1),0,AssignTarget(i+1:AgentNum)];
        end
        tau = 10/k^2;
        % 计算k+1时刻不同任务的效用值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,TargetValue,Assign1,CompleteProb)...
                - TargetUtility(j,TargetValue,Assign2,CompleteProb);% 计算第i个智能体的边际收益
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
%     end
    
    %计算总效用
    for j = 1:TargetNum
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(j, TargetValue, AssignTarget, CompleteProb);
    end
end

MaxGlobalUtility = GlobalUtility(T);

end

function [U] = TargetUtility(target,targetValue,AssignTarget,completeProb)

nAgent = length(find(AssignTarget==target));
U = targetValue(target)*(1-(1-completeProb(target))^nAgent);

end