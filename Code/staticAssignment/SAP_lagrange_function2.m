function [AssignTarget, MaxGlobalUtility, GlobalUtility] = SAP_lagrange_function2(AgentNum, TargetNum, Adjacent, TargetValue,RequireNum,Time_to_go,J_opt,Epis)

ProbDist = ones(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randi(TargetNum,AgentNum,1);%分配任务表
% % 初始任务分配
% % AssignTarget = randperm(TargetNum);
% 
% for i = 1:AgentNum
%     ProbDist(i,:) = ProbDist(i,:)/sum(ProbDist(i,:));
%     e = rand;
%     cp = 0;
%     for j = 1:TargetNum
%         cp = cp + ProbDist(i,j);
%         if e <= cp
%             AssignTarget(i) = j;
%             break
%         end
%     end
% end
% assigntarget = [];
% for j=1:TargetNum
%     assigntarget = [assigntarget;j*ones(RequireNum(j),1)];
% end
% 
% AssignTarget = assigntarget;
AssignMatrix = zeros(AgentNum,TargetNum);
for i = 1:AgentNum
    AssignMatrix(i,AssignTarget(i)) = 1;
end

%% 迭代更新
T = Epis;%总迭代步数
GlobalUtility = zeros(1,T);
lambda = 300;

sap_satisfied = zeros(AgentNum,1);
ex_satisfied = zeros(AgentNum,1);

for j = 1:TargetNum
        GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, AssignTarget, RequireNum, Time_to_go, J_opt, lambda);
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
            AgentUtility(j) = TargetUtility(j,TargetValue,Assign1,RequireNum,Time_to_go,J_opt,lambda)...
                - TargetUtility(j,TargetValue,Assign2,RequireNum,Time_to_go,J_opt,lambda);% 计算第i个智能体的边际收益
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
                if j==AssignTarget(i)
                    sap_satisfied(i) = 1;
                else
                    sap_satisfied(i) = 0;
                end
                AssignTarget(i) = j;
                AssignMatrix(i,j) = 1;
                break
            end
        end
    %end
   neighbors = find(Adjacent(i,:)==1);
   num_neighbors = length(neighbors);
   if sum(sap_satisfied(neighbors))==num_neighbors
       for p=1:num_neighbors
           if p==i
               ex_satisfied(p) = 1;
               continue
           end
           ex_satisfied(neighbors) = 1;
           AssignEx1 = AssignTarget;
           AssignEx2 = AssignTarget;
           task1 = AssignEx2(i);
           AssignEx2(i) = AssignEx2(p);
           AssignEx2(p) = task1;
           AgentU1 = Agent_Utility(i,TargetValue,AssignEx2,RequireNum,Time_to_go,J_opt,lambda)-Agent_Utility(i,TargetValue,AssignEx1,RequireNum,Time_to_go,J_opt,lambda);
           AgentU2 = Agent_Utility(p,TargetValue,AssignEx2,RequireNum,Time_to_go,J_opt,lambda)-Agent_Utility(p,TargetValue,AssignEx1,RequireNum,Time_to_go,J_opt,lambda);
           Net_utility = AgentU1+AgentU2;
           if Net_utility>0
               extask = AssignTarget(i);
               AssignTarget(i) = AssignTarget(p);
               AssignTarget(p) = extask;
               sap_satisfied(i) = 0;
               sap_satisfied(p) = 0;
               ex_satisfied(i) = 0;
               ex_satisfied(p) = 0;
           end
           
       end
   end
    
    
    lagrange_multipler = 0;
    for p=1:TargetNum
        part_missiles = find(AssignTarget==p);
        num_missiles = length(part_missiles);
        lagrange_multipler = lagrange_multipler + min(RequireNum(p)-num_missiles,0);
    end
    %计算总效用
    for j = 1:AgentNum
        s = AssignTarget(j);
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(s, TargetValue, AssignTarget, RequireNum, Time_to_go, J_opt, lambda);
        GlobalUtility(k) = GlobalUtility(k) - lagrange_multipler;
    end
    
    if (sum(sap_satisfied) == AgentNum)&&(sum(ex_satisfied)==AgentNum)
        MaxGlobalUtility = GlobalUtility(k);
        GlobalUtility(k+1:T) = MaxGlobalUtility*ones(1,T-k);
        break
    end
    if k==T
        MaxGlobalUtility = GlobalUtility(T);
    end
end

%-lambda*(sum(RequireNum)-AgentNum);

end

%% 计算效用函数
function [U] = TargetUtility(target,target_value,AssignTarget,require_num,Time_to_go,J_opt,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = require_num(target);
if num_missiles == 0
    con = lambda*min(num_require-num_missiles,0);
    U = max(0,con);
else
    Time_to_go_max = max(Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(J_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = max(target_value(target) - task_cost + con,0);
end

end

function [U] = Agent_Utility(agent,target_value,AssignTarget,RequireNum,Time_to_go,J_opt,lambda)

Assign = AssignTarget;
Assign(agent) = 0;
target = AssignTarget(agent);
U = TargetUtility(target,target_value,AssignTarget,RequireNum,Time_to_go,J_opt,lambda)...
    -TargetUtility(target,target_value,Assign,RequireNum,Time_to_go,J_opt,lambda);
end