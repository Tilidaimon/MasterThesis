function [FinalAssignList, globalutility, time] = HCGFunction2(e,model, Epis)

%% Parameter Initialisation
AgentNum = model.num_missiles;
TargetNum = model.num_targets;
Adjacent = model.Adjacent;

EvolveNum = zeros(1,AgentNum); % an integer variable to represent how many times partition has evolved
TimeStamp = zeros(1,AgentNum); % a uniform-random time stamp

% AssignList = ones(AgentNum); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(AgentNum, TargetNum+1); % Initial Partition
% coAgentNums(:,1) = AgentNum*ones(AgentNum,1); % Initial numbers distribution of cooperated agents
AssignList = randi(TargetNum, AgentNum, 1);

coAgentNums = zeros(TargetNum,1); % Initial Partition
for j=1:TargetNum
    coAgentNums(j) = sum(AssignList == j);
end
Satisfied = zeros(1,AgentNum); % whether or not the agent is satisfied with the partition


%% Decision-making process begins
T = Epis;
globalutility = zeros(1,T);
tic
for t = 1:T
    for i = 1:AgentNum
%         newAssignList = AssignList;
%         newcoAgentNum = coAgentNums;
        best_task = AssignList(i); % current task of agent i
        %         curUtility = AgentTaskUtility(TaskCat(bestTask), coAgentNums(i,bestTask), Distance(i,bestTask), Rmax(bestTask), DesiredNum(bestTask), Rmin(bestTask), epi);
        curUtility = AuxiliaryUtility(e, i, best_task, coAgentNums(best_task), model);
        bestUtility = curUtility;
        
        % update task of agent i
        % choose best task for agent i
        for j = 1:TargetNum
            if AssignList(i) == j
                continue
            else
                %                     tmpUti = AgentTaskUtility(j, TaskCat(j), coAgentNums(i,j)+1, Distance(i,j), Rmax(j), DesiredNum(j), Rmin(j), epi);
                tmpUti = AuxiliaryUtility(e, i, j, coAgentNums(j)+1, model);
                if tmpUti > bestUtility
                    bestUtility = tmpUti;
                    best_task = j;
                end
            end
        end
        % update task
        if bestUtility > curUtility
            curTask = AssignList(i);
            coAgentNums(curTask) = coAgentNums(curTask) - 1;
            AssignList(i) = best_task;
            coAgentNums(best_task) = coAgentNums(best_task) + 1;
            EvolveNum(i) = EvolveNum(i) + 1;
            TimeStamp(i) = rand;
            
        end
        %             Satisfied(i) = Satisfied(i) + 1;
    end
%     AssignList = newAssignList;
%     coAgentNums = newcoAgentNum;
    %     [EvolveNum, TimeStamp, coAgentNums, AssignList, Satisfied] = DMutex(Satisfied, Adjacent, EvolveNum, TimeStamp, coAgentNums, AssignList);
    
    globalutility(t) = GlobalUtility(model,AssignList);
    if sum(Satisfied==20) == AgentNum
        time = toc
        FinalAssignList = AssignList;
        maxutility = globalutility(t);
        globalutility(t+1:T) = maxutility*ones(1,T-t);
        break;
    end
end
time = toc
FinalAssignList = AssignList;
end


%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(e, agent, target, coAgentNum, model)

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