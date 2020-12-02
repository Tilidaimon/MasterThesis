function [FinalAssignList, maxutility, globalutility] = HCGFunctionSAP(AgentNum, TaskNum, Adjacent, Task_Value, Max_limit_num, Time_to_go, J_opt, Epis)

%% Parameter Initialisation
EvolveNum = zeros(1,AgentNum); % an integer variable to represent how many times partition has evolved
TimeStamp = zeros(1,AgentNum); % a uniform-random time stamp

AssignList = randi(TaskNum,AgentNum,1); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(TaskNum,1); % Initial Partition
for i = 1:TaskNum
    coAgentNums(i) = length(find(AssignList==i));
end
% coAgentNums(:,1) = AgentNum*ones(AgentNum,1); % Initial numbers distribution of cooperated agents


%% Decision-making process begins
T = Epis; 

globalutility = zeros(1,T);
for t = 1:T
    localutility = zeros(1,AgentNum);
    i = randi(AgentNum);
    best_task = AssignList(i); % current task of agent i
    curUtility = AuxiliaryUtility(coAgentNums(best_task), Task_Value(best_task), Max_limit_num(best_task), Time_to_go(i,best_task), J_opt(i,best_task));
    % update task of agent i
    % choose best task for agent i
    tmpUti = zeros(1,TaskNum);
    tau = 10/t^2;
    for j = 1:TaskNum
        if AssignList(i) == j
            tmpUti(j) = curUtility;
        elseif Adjacent(i,j) == 1
            %                     tmpUti = AgentTaskUtility(j, TaskCat(j), coAgentNums(i,j)+1, Distance(i,j), Rmax(j), DesiredNum(j), Rmin(j), epi);
            tmpUti(j) = AuxiliaryUtility(coAgentNums(j)+1, Task_Value(j), Max_limit_num(j), Time_to_go(i,j), J_opt(i,j));
        end
    end
    % update task
    AU = tmpUti/tau;
    newProb = exp(AU)/sum(exp(AU));
    flag = 0;
    for j = 1:TaskNum
        if isnan(newProb(j))
            flag = 1;
            break;
        end
    end
    if flag == 0
        ProbDist = newProb;
    else
        AU = AU/max(abs(AU))*t/2;
        ProbDist = exp(AU)/sum(exp(AU));
    end
    e = rand;
    cp = 0;
    for j = 1:TaskNum
        cp = cp + ProbDist(j);
        if e <= cp
            curTask = AssignList(i);
            coAgentNums(curTask) = coAgentNums(curTask) - 1;
            AssignList(i) = j;
            coAgentNums(j) = coAgentNums(j) + 1;
            break
        end
    end

    for i=1:AgentNum
        best_task = AssignList(i); % current task of agent i
        localutility(i) = AuxiliaryUtility(coAgentNums(best_task), Task_Value(best_task), Max_limit_num(best_task), Time_to_go(i,best_task), J_opt(i,best_task));
    end
    globalutility(t) = sum(localutility);
end

FinalAssignList = AssignList;
maxutility = globalutility(T);

end


%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(coAgentNum, task_value, max_limit_num, time_to_go, j_opt)

e = 3;
rmin = task_value;
if coAgentNum > max_limit_num
    AU = 0;
else
    %R = (task_value/desired_num) * exp(-coAgentNum/desired_num+1);
    R = rmin*log(coAgentNum+e-1)/coAgentNum;
    cost = time_to_go + j_opt;
    AU = R - cost;
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