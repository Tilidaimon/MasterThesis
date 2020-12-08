% 带有Lagrange乘子法的SAP算法
function assign_plan = SAPLagrangeFunction(model,Epis,initial_plan)
% 输入： model:模型
%        Epis:迭代次数
% 输出： assign_plan：分配方案

AgentNum = model.Missiles.alive_num;
TargetNum = model.Targets.alive_num;
ProbDist = ones(AgentNum, TargetNum);%智能体选择任务概率分布表存储
assign_plan = initial_plan;%randi(model.num_targets,model.num_missiles,1);%分配任务表

%% 迭代更新
T = Epis;%总迭代步数
GlobalUtility = zeros(1,T);
lambda = 300;

sap_satisfied = zeros(AgentNum,1);
ex_satisfied = zeros(AgentNum,1);

for j = 1:TargetNum
    GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, model, assign_plan, lambda);
end

for k = 2:T
    %for i = 1:AgentNum
        
        i = randi(model.Missiles.alive_num);
        % 假设第i个智能体不选任务
        Assign2 = assign_plan;
        Assign2(i) = 0;
        
        tau = 10/k^2;
        % 计算k+1时刻不同任务的效用值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        for j = 1:TargetNum
            if model.Missiles.target_set(i,j) == true
                Assign1 = Assign2;
                Assign1(i) = j;
                % 计算第i个智能体的边际收益
                AgentUtility(j) = TargetUtility(j,model,Assign1,lambda) - TargetUtility(j,model,Assign2,lambda);
            else
                AgentUtility(j) = 0;
            end
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
                if j==assign_plan(i)
                    sap_satisfied(i) = 1;
                else
                    sap_satisfied(i) = 0;
                end
                assign_plan(i) = j;
                break
            end
        end
    %end
   neighbors = find(model.Adjacent(i,:)==true);
   num_neighbors = length(neighbors);
   if sum(sap_satisfied(neighbors))==num_neighbors
       for p=1:num_neighbors
           if p==i
               ex_satisfied(p) = 1;
               continue
           end
           ex_satisfied(neighbors) = 1;
           AssignEx1 = assign_plan;
           AssignEx2 = assign_plan;
           task1 = AssignEx2(i);
           AssignEx2(i) = AssignEx2(p);
           AssignEx2(p) = task1;
           AgentU1 = Agent_Utility(i,model,AssignEx2,lambda)-Agent_Utility(i,model,AssignEx1,lambda);
           AgentU2 = Agent_Utility(p,model,AssignEx2,lambda)-Agent_Utility(p,model,AssignEx1,lambda);
           Net_utility = AgentU1+AgentU2;
           if Net_utility>0
               extask = assign_plan(i);
               assign_plan(i) = assign_plan(p);
               assign_plan(p) = extask;
               sap_satisfied(i) = 0;
               sap_satisfied(p) = 0;
               ex_satisfied(i) = 0;
               ex_satisfied(p) = 0;
           end
           
       end
   end
    
    
    lagrange_multipler = 0;
    for p=1:TargetNum
        part_missiles = find(assign_plan==p);
        num_missiles = length(part_missiles);
        lagrange_multipler = lagrange_multipler + min(model.target_require_num_list(p)-num_missiles,0);
    end
    %计算总效用
    for j = 1:AgentNum
        s = assign_plan(j);
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(s, model, assign_plan, lambda);
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
function [U] = TargetUtility(target, model, assign_plan, lambda)

part_missiles = find(assign_plan == target);
num_missiles = length(part_missiles);
num_require = max(0,model.target_require_num_list(target) - model.Targets.num_attacked(target));
if num_missiles == 0
    con = lambda*min(num_require-num_missiles,0);
    U = con;
else
    Time_to_go_max = max(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = model.Targets.value(target) - task_cost + con;
end
end

function [U] = Agent_Utility(agent,model,assign_plan,lambda)

Assign = assign_plan;
Assign(agent) = 0;
target = assign_plan(agent);
U = TargetUtility(target,model,assign_plan,lambda)...
    -TargetUtility(target,model,Assign,lambda);
end