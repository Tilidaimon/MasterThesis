function [AssignTarget,Ug,time] = SAP_lagrange_function2(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;
Adjacent = model.Adjacent;

ProbDist = ones(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randi(TargetNum,AgentNum,1);%分配任务表


%% 迭代更新
T = Epis;%总迭代步数
%Ug = zeros(1,T);
lambda = max(model.Targets.value);

sap_satisfied = zeros(AgentNum,1);
ex_satisfied = zeros(AgentNum,1);
Ug = zeros(T,1);
Ug(1) = GlobalUtility(model,AssignTarget);
tic
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
            AgentUtility(j) = TargetUtility(j,model,Assign1,lambda)- TargetUtility(j,model,Assign2,lambda);% 计算第i个智能体的边际收益
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
                    sap_satisfied(i) = sap_satisfied(i) + 1;
                else
                    sap_satisfied(i) = 0;
                end
                AssignTarget(i) = j;
                break
            end
        end
    %end
   neighbors = find(Adjacent(i,:)==1);
   num_neighbors = length(neighbors);
   if sum(sap_satisfied(neighbors)~=0)==num_neighbors
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
           AgentU1 = Agent_Utility(i,model,AssignEx2,lambda)-Agent_Utility(i,model,AssignEx1,lambda);
           AgentU2 = Agent_Utility(p,model,AssignEx2,lambda)-Agent_Utility(p,model,AssignEx1,lambda);
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
   Ug(k) = GlobalUtility(model,AssignTarget);
   if(sum(sap_satisfied>20)==model.num_missiles && sum(ex_satisfied)==model.num_missiles)
       time = toc
       Ug(k+1:T) = Ug(k)*ones(T-k,1);
       break
   end
end
time = toc
%-lambda*(sum(RequireNum)-AgentNum);

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

function [U] = Agent_Utility(agent,model,AssignTarget,lambda)

Assign = AssignTarget;
Assign(agent) = 0;
target = AssignTarget(agent);
U = TargetUtility(target,model,AssignTarget,lambda)-TargetUtility(target,model,Assign,lambda);
end

%% 分配结果效用
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