function [FinalAssignList, globalutility, time] = PureHCGFunction(epsilon, model, Epis)

%% Parameter Initialisation
AgentNum = model.num_missiles;
TargetNum = model.num_targets;
Adjacent = model.Adjacent;

AssignList = randi(TargetNum,AgentNum,1); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(TargetNum,1); % Initial Partition
for i = 1:TargetNum
    coAgentNums(i) = sum(AssignList==i);
end
% coAgentNums(:,1) = AgentNum*ones(AgentNum,1); % Initial numbers distribution of cooperated agents


%% Decision-making process begins
T = Epis; 
hcg_satisfied = false(AgentNum,1);
ex_satisfied = false(AgentNum,1);

globalutility = zeros(1,T);
tic
for t = 1:T
    for i=1:AgentNum
    best_task = AssignList(i); % current task of agent i
    curUtility = AuxiliaryUtility(i, best_task, coAgentNums(best_task), model, epsilon);
    % update task of agent i
    % choose best task for agent i
    tmpUti = zeros(TargetNum,1);
    tau = 10/t^2;
    for j = 1:TargetNum
        if AssignList(i) == j
            tmpUti(j) = curUtility;
        elseif Adjacent(i,j) == 1
            %                     tmpUti = AgentTaskUtility(j, TaskCat(j), coAgentNums(i,j)+1, Distance(i,j), Rmax(j), DesiredNum(j), Rmin(j), epi);
            tmpUti(j) = AuxiliaryUtility(i, j, coAgentNums(j)+1, model, epsilon);
        end
    end
    
    new_task = find(tmpUti==max(tmpUti));
    if length(new_task)>1
        best_task = new_task(randi(length(new_task)));
    else
        best_task = new_task;
    end
    if best_task == AssignList(i)
        hcg_satisfied(i) = hcg_satisfied(i) +1;
    else
        hcg_satisfied(i) = 0;
    end
    coAgentNums(AssignList(i)) = coAgentNums(AssignList(i)) - 1;
    AssignList(i) = best_task;
    coAgentNums(best_task) = coAgentNums(best_task) + 1;
    end
    globalutility(t) = GlobalUtility(model,AssignList);
    if(sum(hcg_satisfied==20)==model.num_missiles)% && sum(ex_satisfied)==model.num_missiles)
       globalutility(k+1:T) = globalutility(k)*ones(T-k,1);
       time=toc
       FinalAssignList = AssignList;
       break
   end
end
time = toc
FinalAssignList = AssignList;
end


%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(agent, target, coAgentNum, model, e)

rmin = model.Targets.value(target);
if coAgentNum > model.target_require_num_list(target)
    AU = 0;
else
    %R = (task_value/desired_num) * exp(-coAgentNum/desired_num+1);
    R = rmin*log(coAgentNum+e-1)/coAgentNum;
    cost = model.Time_to_go(agent,target) + model.Energy_opt(agent,target);
    AU = max(0,R - cost);
end
end

%% Distributed Mutual Exclusion Subroutine Function
function [evolved, timestamp, coAgentNums, AssignList, satisfied] ...
    = DMutex(satisfied, adjacent, evolved, timestamp, coAgentNums, AssignList)

[agentnum, tasknum] = size(coAgentNums);
newevol = evolved;
newts = timestamp;
newlist = AssignList;

for i = 1:agentnum
    satisfied(i) = 1;
    for j = 1:agentnum
        if adjacent(i,j) == 1
            if (evolved(j) > evolved(i))||((evolved(j)==evolved(i))&&(timestamp(j)>timestamp(i)))
                newevol(i) = evolved(j);
                newts(i) = timestamp(j);
                newlist(i,:) = AssignList(j,:);
                satisfied(i) = 0;
            end
        end
    end
end

evolved = newevol;
timestamp = newts;
AssignList = newlist;
newconum = zeros(agentnum, tasknum);
for i = 1:agentnum
    for j = 1:agentnum
        newconum(i,AssignList(i,j)) = newconum(i,AssignList(i,j)) + 1;
    end
end
coAgentNums = newconum;

end

%% 全局效用
function [Ug] = GlobalUtility(model,plan)
Nt = model.num_targets;
Ut = zeros(Nt,1);
flag = 1;
for j=1:Nt
    part_missiles = find(plan == j);
    num_missiles = length(part_missiles);
    if num_missiles > model.target_require_num_list(j)
        flag = 0;
    end
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = sum(model.Time_to_go(part_missiles,j));
        J_sum = sum(model.Energy_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,model.Targets.value(j)/model.target_require_num_list(j)*num_missiles - cost);
    end
    
end
if flag == 0
    Ug = 0;
else
    Ug = sum(Ut);
end
end