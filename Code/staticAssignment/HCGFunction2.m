function [FinalAssignList, maxutility, globalutility] = HCGFunction2(AgentNum, TaskNum, Adjacent, Task_Value, Desired_Num, Max_limit_num, Time_to_go, J_opt, Epis)

%% Parameter Initialisation
EvolveNum = zeros(1,AgentNum); % an integer variable to represent how many times partition has evolved
TimeStamp = zeros(1,AgentNum); % a uniform-random time stamp

AssignList = randi(TaskNum,AgentNum,AgentNum); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(AgentNum, TaskNum); % Initial Partition
for j=1:TaskNum
    for i=1:AgentNum
        coAgentNums(i,j) = length(find(AssignList(i,:)==j));
    end
end% Initial numbers distribution of cooperated agents
Satisfied = zeros(1,AgentNum); % whether or not the agent is satisfied with the partition


%% Decision-making process begins
T = Epis; 

globalutility = zeros(1,T);
for t = 1:T
    localutility = zeros(1,AgentNum);
    for i = 1:AgentNum
        
        best_task = AssignList(i,i); % current task of agent i
%         curUtility = AgentTaskUtility(TaskCat(bestTask), coAgentNums(i,bestTask), Distance(i,bestTask), Rmax(bestTask), DesiredNum(bestTask), Rmin(bestTask), epi);
        curUtility = AuxiliaryUtility(coAgentNums(i,best_task), Task_Value(best_task), Desired_Num(best_task), Max_limit_num(best_task), Time_to_go(i,best_task), J_opt(i,best_task));
        bestUtility = curUtility;
        localutility(i) = curUtility;
        globalutility(t) = sum(localutility);
        % update task of agent i
        if Satisfied(i) == 0
            % choose best task for agent i
            for j = 1:TaskNum
                if AssignList(i,i) == j
                    continue
                else
%                     tmpUti = AgentTaskUtility(j, TaskCat(j), coAgentNums(i,j)+1, Distance(i,j), Rmax(j), DesiredNum(j), Rmin(j), epi);
                    tmpUti = AuxiliaryUtility(coAgentNums(i,j)+1, Task_Value(j), Desired_Num(j), Max_limit_num(j), Time_to_go(i,j), J_opt(i,j));
                    if tmpUti > bestUtility
                        bestUtility = tmpUti;
                        best_task = j;
                    end
                end
            end
            % update task
            if bestUtility > curUtility
                curTask = AssignList(i,i);
                coAgentNums(i, curTask) = coAgentNums(i, curTask) - 1;
                AssignList(i,i) = best_task;
                coAgentNums(i, best_task) = coAgentNums(i, best_task) + 1;
                EvolveNum(i) = EvolveNum(i) + 1;
                TimeStamp(i) = rand;
                
            end
            Satisfied(i) = 1;
        end
    end
    
    [EvolveNum, coAgentNums, AssignList, Satisfied] = ...
        DMutex(Satisfied, Adjacent, EvolveNum, coAgentNums, AssignList, Task_Value, Desired_Num, Max_limit_num, Time_to_go, J_opt);
    if sum(Satisfied) == AgentNum
        FinalAssignList = AssignList(1,:)';
        maxutility = globalutility(t);
        globalutility(t+1:T) = maxutility*ones(1,T-t);
        break;
    end
end

% FinalAssignList = AssignList;
% maxutility = globalutility(T);


end


%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(coAgentNum, task_value, desired_num, max_limit_num, time_to_go, j_opt)

if coAgentNum > max_limit_num
    AU = 0;
else
    R = (task_value/desired_num) * exp(-coAgentNum/desired_num+1);
    cost = time_to_go + j_opt;
    AU = R - cost;
end
end

%% Distributed Mutual Exclusion Subroutine Function
function [evolved, coAgentNums, AssignList, satisfied] ...
    = DMutex(satisfied, adjacent, evolved, coAgentNums, AssignList, task_value, desired_num, max_limit_num, time_to_go, j_opt)

[agentnum, tasknum] = size(coAgentNums);
newevol = evolved;
% newts = timestamp;
newlist = AssignList;

for i = 1:agentnum
    satisfied(i) = 1;
    for j = 1:agentnum
        if i==j
            continue 
        elseif adjacent(i,j) == 1
            taski_for_i = AssignList(i,i);
            taskj_for_i = AssignList(j,i);
            Ui_for_i = AuxiliaryUtility(coAgentNums(i,taski_for_i), task_value(taski_for_i), ...
                desired_num(taski_for_i),max_limit_num(taski_for_i),time_to_go(taski_for_i),j_opt(taski_for_i));
            Uj_for_i = AuxiliaryUtility(coAgentNums(j,taskj_for_i), task_value(taskj_for_i), ...
                desired_num(taskj_for_i),max_limit_num(taskj_for_i),time_to_go(taskj_for_i),j_opt(taskj_for_i));
            if (evolved(j) > evolved(i))||((evolved(j)==evolved(i))&&(Uj_for_i>Ui_for_i))
                newevol(i) = evolved(j);
%                 newts(i) = timestamp(j);
                newlist(i,:) = AssignList(j,:);
                satisfied(i) = 0;
            end
        end
    end
end

evolved = newevol;
% timestamp = newts;
AssignList = newlist;
newconum = zeros(agentnum, tasknum);
for i = 1:agentnum
    for j = 1:agentnum
        newconum(i,AssignList(i,j)) = newconum(i,AssignList(i,j)) + 1;
    end
end
coAgentNums = newconum;

end