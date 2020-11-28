function [AssignTarget, MaxGlobalUtility, GlobalUtility] = GRMfuncSimple(AgentNum,TargetNum,TargetValue,CompleteProb)

ProbDist = rand(AgentNum, TargetNum);%智能体选择任务概率分布表存储
AssignTarget = randperm(TargetNum);%分配任务表
AssignTarget = AssignTarget(1:AgentNum);

%% 迭代更新
T = 500;%总迭代步数
GlobalUtility = zeros(1,T);
RegretMatrix = zeros(AgentNum,TargetNum);%后悔矩阵
for j = 1:TargetNum
        GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, AssignTarget, CompleteProb);
end
% 在第k步迭代时，计算每个智能体的后悔值，即前k-1步没有选择任务l的后悔程度，后悔程度以任务l带来的效用提高衡量

rho = 0.1;% 遗忘因子
alpha = 0.5;% 惰性因子

for k = 2:T
    for i = 1:AgentNum
        % 假设第i个智能体不选任务
        if i==1
            Assign2 = [0,AssignTarget(2:AgentNum)];
        elseif i == AgentNum
            Assign2 = [AssignTarget(1:AgentNum-1),0];
        else
            Assign2 = [AssignTarget(1:i-1),0,AssignTarget(i+1:AgentNum)];
        end
        
        % 计算k+1时刻不同任务的效用值和后悔值
        AgentUtility = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的效用值
        Regret = zeros(1,TargetNum);%智能体i在k+1时刻执行不同任务的后悔值
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,TargetValue,Assign1,CompleteProb)...
                - TargetUtility(j,TargetValue,Assign2,CompleteProb);% 计算第i个智能体的边际收益
            %计算k+1时刻的后悔值
            Regret(j) = (1-rho)*RegretMatrix(i,j) + rho*(AgentUtility(j)-AgentUtility(AssignTarget(i)));
        end
        RegretMatrix(i,:) = Regret;
        % 计算智能体i在k+1时刻执行不同任务的概率分布
        Regret(Regret<0) = 0;% 后悔值小于0的值取0
        %RM = Regret/sum(Regret);% 后悔匹配函数
        RM = zeros(1,TargetNum);
%         RM = Regret/sum(Regret);
        tau = 10/k^2;
        deno = 0;
        for s = 1:TargetNum
            if Regret(s) ~= 0
                deno = deno + exp(Regret(s)/tau);
            end
        end
        for s = 1:TargetNum
            if Regret(s)>0
                if isinf(deno)
                    deno = 0;
                    for p = 1:TargetNum
                        if Regret(p) ~= 0
                            deno = deno + exp(Regret(p)/tau-k);
                        end
                    end
                    
                end
                RM(s) = exp(Regret(s)/tau-k)/deno;
            else
                RM(s) = 0;
            end
        end
        vertex = zeros(1,TargetNum);
        vertex(AssignTarget(i)) = 1;
        ProbDist(i,:) = alpha*RM + (1-alpha)*vertex;% 归一化值作为k+1时刻选择任务概率分布
        
        % 任务分配
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
