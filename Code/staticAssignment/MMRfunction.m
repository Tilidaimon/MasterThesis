function [AssignRRM,UtilityRRM] = RRMfunction(AgentNum,TargetNum,TargetValue,CompleteProb)

ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randperm(TargetNum);%分配任务表
AssignTarget = AssignTarget(1:AgentNum);

%% RRM对比
AssignRRM = zeros(1,AgentNum);
AssignedAgent = zeros(1,TargetNum);
SurviveProb = zeros(1,TargetNum);

for j = 1:TargetNum
    SurviveProb(j) = TargetValue(j)*(1-CompleteProb(j))^AssignedAgent(j);
end
for i = 1:AgentNum
    [~,k] = max(SurviveProb.*CompleteProb);
    AssignedAgent(k) = AssignedAgent(k) + 1;
    AssignRRM(i) = k;
    SurviveProb(k) = SurviveProb(k)*(1-CompleteProb(k));
end
 UtilityRRM = 0;
for j = 1:TargetNum
    UtilityRRM = UtilityRRM + TargetUtility(j, TargetValue(j), AssignRRM, CompleteProb(j));
end

end

function [U] = TargetUtility(target,targetValue,AssignTarget,completeProb)

nAgent = length(find(AssignTarget==target));
U = targetValue*(1-(1-completeProb)^nAgent);

end