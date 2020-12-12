function [AssignTarget, Ug, time] = GRMfunction(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;


ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randi(TargetNum,AgentNum,1);%分配任务表
%% 迭代更新
T = Epis;%总迭代步数
% Ug = zeros(1,T);
RegretMatrix = zeros(AgentNum,TargetNum);%后悔矩阵
Ug(1) = GlobalUtility(model,AssignTarget);
% 在第k步迭代时，计算每个智能体的后悔值，即前k-1步没有选择任务l的后悔程度，后悔程度以任务l带来的效用提高衡量

rho = 0.1;% 遗忘因子
alpha = 0.5;% 惰性因子
lambda = 30;
satisfied = false(AgentNum,1);
tic
for k = 2:T
    for i = 1:AgentNum
        % 假设第i个智能体不选任务
        Assign2 = AssignTarget;
        Assign2(i) = 0;
        
        % 计算k+1时刻不同任务的效用值和后悔值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        Regret = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的后悔值
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,model,Assign1,lambda)...
                - TargetUtility(j,model,Assign2,lambda);% 计算第i个智能体的边际收益
            %计算k+1时刻的后悔值
            Regret(j) = (1-rho)*RegretMatrix(i,j) + rho*(AgentUtility(j)-AgentUtility(AssignTarget(i)));
        end
        RegretMatrix(i,:) = Regret;
        % 计算智能体i在k+1时刻执行不同任务的概率分布
        Regret(Regret<0) = 0;% 后悔值小于0的值取0
        %RM = Regret/sum(Regret);% 后悔匹配函数
        RM = zeros(1,TargetNum);
        tau = 10/k^2;
        deno = 0;
        flag = 0;
        for s = 1:TargetNum
            if Regret(s) ~= 0
                deno = deno + exp(Regret(s)/tau);
                
            end
        end
        if isnan(deno)
            flag = 1;
        end
%         if flag==1
%             
%             Rm = Regrets/tau;
%             maxR = max(abs(Rm));
%             modiRM = Rm/maxR*k/2;
%             for s = 1:TargetNum
%                 if Regret(s) ~= 0
%                     deno = deno + exp(modiRM);
%                 end
%             end
%             for s = 1:TargetNum
%                 if Regret(s)>0
%                     RM(s) = exp(modiRM)/deno;
%                 else
%                     RM(s) = 0;
%                 end
%             end
%         else
            for s = 1:TargetNum
                if Regret(s)>0
                    RM(s) = exp(Regret(s)/tau)/deno;
                else
                    RM(s) = 0;
                end
            end
%         end
        vertex = zeros(1,TargetNum);
        vertex(AssignTarget(i)) = 1;
        ProbDist(i,:) = alpha*RM + (1-alpha)*vertex;% 归一化值作为k+1时刻选择任务概率分布
        
        % 任务分配
        e = rand;
        cp = 0;
        for j = 1:TargetNum
            cp = cp + ProbDist(i,j);
            if e <= cp
                if j==AssignTarget(i)
                    satisfied(i) = 1;
                else
                    satisfied(i) = 0;
                end
                AssignTarget(i) = j;
                
                break
            end
        end
    end
    
    %计算总效用
    Ug(k) = GlobalUtility(model,AssignTarget);
%     if sum(satisfied)==AgentNum
%         break
%     end
end
time = toc;

end

%计算效用函数
function [U] = TargetUtility(target,model,AssignTarget,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = model.target_require_num_list(target);
if num_missiles == 0
    con = lambda*min(num_require-num_missiles,0);
    U = max(0,con);
else
    Time_to_go_max = max(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = max(model.Targets.value(target) - task_cost + con,0);
end

end

%% 全局效用
function [Ug] = GlobalUtility(model,plan)
Nt = model.num_targets;
Ut = zeros(Nt,1);

for j=1:Nt
    part_missiles = find(plan == j);
    num_missiles = length(part_missiles);
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = max(model.Time_to_go(part_missiles,j));
        J_sum = sum(model.Energy_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,model.Targets.value(j) - cost);
    end
    
end
Ug = sum(Ut);
end