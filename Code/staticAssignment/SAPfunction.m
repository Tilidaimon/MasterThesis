function [AssignTarget, MaxGlobalUtility, GlobalUtility] = SAPfunction(TargetValue,RequireNum,Time_to_go,J_opt,Epis)

[AgentNum,TargetNum] = size(Time_to_go);

ProbDist = rand(AgentNum, TargetNum);%������ѡ��������ʷֲ���洢
assigntarget = [];
for j=1:TargetNum
    assigntarget = [assigntarget, j*ones(1,RequireNum(j))];
end

AssignTarget = assigntarget;

%% ��������
T = Epis;%�ܵ�������
GlobalUtility = zeros(1,T);


for j = 1:TargetNum
    GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, AssignTarget, Time_to_go, J_opt);
end

for k = 2:T
%     selectAgentNum = 100;
%     reorderAgent = randperm(AgentNum);
%     selectAgent = reorderAgent(1:selectAgentNum);
%     for p = 1:selectAgentNum
%     for i = 1:AgentNum
        i = randi(AgentNum);
%         i = selectAgent(p);
        Assign2 = AssignTarget;
        % �����i�������岻ѡ����
        tau = 10/k^2;
        % ����k+1ʱ�̲�ͬ�����Ч��ֵ
        AgentUdist = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ�����Ч��ֵ
        for j = 1:TargetNum
            Assign1 = Assign2;
            exindex = find(Assign2==j);
            if length(exindex)>1
                exindex = exindex(randi(length(exindex)));
            end
            extarget = Assign1(i);
            Assign1(i) = j;
            Assign1(exindex) = extarget;
            AgentU1 = AgentUtility(i,TargetValue,Assign1,Time_to_go,J_opt)-AgentUtility(i,TargetValue,Assign2,Time_to_go,J_opt);
            AgentU2 = AgentUtility(exindex,TargetValue,Assign1,Time_to_go,J_opt)-AgentUtility(exindex,TargetValue,Assign2,Time_to_go,J_opt);
            AgentUdist(j) = AgentU1+AgentU2;% �����i��������ı߼�����
        end
        %������ʷֲ����˴����������Ϊԭ��ѡ����ʷֲ�
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
        
        % �������
        e = rand;
        cp = 0;
        for j = 1:TargetNum
            cp = cp + ProbDist(i,j);
            if e <= cp
                extarget = AssignTarget(i);
                exindexs = find(AssignTarget==j);
                if length(exindex)>1
                    exindex = exindexs(randi(length(enindexs)));
                else
                    exindex = exindexs;
                end
                AssignTarget(exindex) = extarget;
                AssignTarget(i) = j;
                break
            end
        end
%     end
    
    %������Ч��
    for j = 1:TargetNum
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(j, TargetValue, AssignTarget, Time_to_go, J_opt);
    end
end

MaxGlobalUtility = GlobalUtility(T);

end

%����Ч�ú���
function [U] = TargetUtility(target,target_value,AssignTarget,Time_to_go,J_opt)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
if num_missiles == 0
    task_cost = 0;
else
    Time_to_go_max = max(Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(J_opt(part_missiles,target));
end

U = target_value(target) - task_cost;

end

function [U] = AgentUtility(agent,target_value,AssignTarget,Time_to_go,J_opt)

Assign = AssignTarget;
Assign(agent) = 0;
target = AssignTarget(agent);
U = TargetUtility(target,target_value,AssignTarget,Time_to_go,J_opt)...
    -TargetUtility(target,target_value,Assign,Time_to_go,J_opt);
end