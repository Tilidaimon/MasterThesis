function [FinalAssignList, globalutility, maxutility] = HedonicGameFunction(AgentNum, TaskNum, TaskCat, Distance, Adjacent, TaskValue, DesiredNum, epi, CompleteProb)

%% Parameter Initialisation
EvolveNum = zeros(1,AgentNum); % an integer variable to represent how many times partition has evolved
TimeStamp = zeros(1,AgentNum); % a uniform-random time stamp

AssignList = ones(AgentNum); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(AgentNum, TaskNum); % Initial Partition
coAgentNums(:,1) = AgentNum*ones(AgentNum,1); % Initial numbers distribution of cooperated agents
Satisfied = zeros(1,AgentNum); % whether or not the agent is satisfied with the partition

Rlimit = ones(1,TaskNum);

%% Decision-making process begins
T = 500; 
globalutility = zeros(1,T);
for t = 1:T
    localutility = zeros(1,AgentNum);
    for i = 1:AgentNum
        
        bestTask = AssignList(i,i); % current task of agent i
%         curUtility = AgentTaskUtility(TaskCat(bestTask), coAgentNums(i,bestTask), Distance(i,bestTask), Rmax(bestTask), DesiredNum(bestTask), Rmin(bestTask), epi);
        curUtility = AuxiliaryUtility(i, bestTask, Rlimit, TaskCat, coAgentNums(i,bestTask), Distance, TaskValue, DesiredNum, epi, CompleteProb);
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
                    tmpUti = AuxiliaryUtility(i, j, Rlimit, TaskCat, coAgentNums(i,j)+1, Distance, TaskValue, DesiredNum, epi, CompleteProb);
                    if tmpUti > bestUtility
                        bestUtility = tmpUti;
                        bestTask = j;
                    end
                end
            end
            % update task
            if bestUtility > curUtility
                curTask = AssignList(i,i);
                coAgentNums(i, curTask) = coAgentNums(i, curTask) - 1;
                AssignList(i,i) = bestTask;
                coAgentNums(i, bestTask) = coAgentNums(i, bestTask) + 1;
                EvolveNum(i) = EvolveNum(i) + 1;
                TimeStamp(i) = rand;
                
            end
            Satisfied(i) = 1;
        end
    end
    
    [EvolveNum, TimeStamp, coAgentNums, AssignList, Satisfied] = DMutex(Satisfied, Adjacent, EvolveNum, TimeStamp, coAgentNums, AssignList);
end

FinalAssignList = AssignList(1,:);
maxutility = globalutility(T);

end

%% Utility function
function [U] = AgentTaskUtility(taskCat, coAgentNum, cost, targetvalue, DesiredNum, epi, completeprob)

% if task == 1
%     U = 0; % Utility is zero if no task
% else
    % Reward according to task category
    if taskCat == 1
        R = (targetvalue/DesiredNum) * exp(-coAgentNum/DesiredNum+1);
    else
        R = 0.5 * targetvalue * log(coAgentNum + epi - 1)/(coAgentNum*log(epi));
    end
    
    U = (R - cost)*completeprob;
% end
end

%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(agent, task, Rlimit, taskCat, coAgentNum, distance, taskvalue, DesiredNum, epi, completeprob)

% beta = 50;
% maxU = 0;
% [~,tasknum] = size(coAgentNum);
% if coAgentNum < Rlimit
%     for i = 1:tasknum
%         u = AgentTaskUtility(taskCat(i), Rlimit(i)+1, distance(agent,i), taskvalue(task), DesiredNum(task),epi,completeprob(task));
%         if u>maxU
%             maxU = u;
%         end
%     end
%     AU = maxU + beta;

if coAgentNum > Rlimit(task)
    AU = 0;
else
    AU = AgentTaskUtility(taskCat(task), coAgentNum, distance(agent,task), taskvalue(task), DesiredNum(task),epi,completeprob(agent,task));
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