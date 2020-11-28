function [AssignRRM,UtilityRRM] = RRMfunction(AgentNum,TargetNum,TargetValue,CompleteProb)

ProbDist = rand(AgentNum, TargetNum);%������ѡ��������ʷֲ���洢
AssignTarget = randperm(TargetNum);%���������
AssignTarget = AssignTarget(1:AgentNum);

%% RRM�Ա�
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