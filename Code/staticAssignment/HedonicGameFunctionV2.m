function [FinalAssignList, globalutility, maxutility] = HedonicGameFunctionV2(AgentNum, TaskNum, Adjacent, Utility)

%% Parameter Initialisation
EvolveNum = zeros(1,AgentNum); % an integer variable to represent how many times partition has evolved
TimeStamp = zeros(1,AgentNum); % a uniform-random time stamp

% AssignList = ones(AgentNum); % Assign List for each agent
AssignList = randi(TaskNum, AgentNum);
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
%         curUtility = AuxiliaryUtility(i, bestTask, Rlimit, TaskCat, coAgentNums(i,bestTask), Distance, TaskValue, DesiredNum, epi, CompleteProb);
        curUtility = AuxiliaryUtility(i, bestTask, Rlimit, coAgentNums(i,bestTask), Utility);
        bestUtility = curUtility;
        localutility(i) = curUtility;
        
        % update task of agent i
        if Satisfied(i) == 0
            % choose best task for agent i
            for j = 1:TaskNum
                if AssignList(i,i) == j
                    continue
                else
%                     tmpUti = AuxiliaryUtility(i, j, Rlimit, TaskCat, coAgentNums(i,j)+1, Distance, TaskValue, DesiredNum, epi, CompleteProb);
                    tmpUti = AuxiliaryUtility(i, j, Rlimit, coAgentNums(i,j)+1, Utility);
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
    
    for i = 1:AgentNum-1
        for j = (i+1):AgentNum
            if Adjacent(i,j) == 1
                [AssignList(i,:),AssignList(j,:),Satisfied] = ZeuthenNegotiation(i,j,AssignList(i,:),AssignList(j,:), TaskNum, Rlimit, Utility, Satisfied);
            end
        end
    end
    globalutility(t) = sum(localutility);
%     [EvolveNum, TimeStamp, coAgentNums, AssignList, Satisfied] = DMutex(Satisfied, Adjacent, EvolveNum, TimeStamp, coAgentNums, AssignList);
end

FinalAssignList = AssignList(1,:);
maxutility = globalutility(T);

end


%% Auxiliary individual utility function
function [AU] = AuxiliaryUtility(agent, task, Rlimit, coAgentNum, Utility)

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
    AU = Utility(agent,task)*exp(-coAgentNum^2);
else
    AU = Utility(agent,task);
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


%% Zeuthen Negotiation Function
function [newassigni,newassignj,satisfied] = ZeuthenNegotiation(i, j, assigni, assignj, TaskNum, Rlimit, Utility, satisfied)

[riski,utilityi1,utilityi2] = RiskCalculation(i,assigni,assignj,Rlimit,Utility);
[riskj,utilityj1,utilityj2] = RiskCalculation(j,assignj,assigni,Rlimit,Utility);

if utilityi2>=utilityi1
    newassigni = assignj;
    newassignj = assignj;
    satisfied(i) = 0;
end

if utilityj2>=utilityj1
    newassignj = assigni;
    newassigni = assigni;
    satisfied(j) = 0;
end

if (utilityi2<utilityi1) && (utilityj2<utilityj1)
    satisfied(i) = 0;
    satisfied(j) = 0;
    newassigni = assigni;
    newassignj = assignj;
    if riski>riskj % j风险小，j做让步
        minriskj = 1;
        for k = 1:TaskNum
            assigntmp = assignj;
            assigntmp(i) = k;
            risktmpi = RiskCalculation(i,assigni,assigntmp,Rlimit,Utility);
            risktmpj = RiskCalculation(j,assignj,assigntmp,Rlimit,Utility);
            if (risktmpi<riski) && (risktmpj<=minriskj)
                newassignj(i) = k;
                newassigni = newassignj;
                minriskj = risktmpj;
                riski = risktmpi;
            end
        end
    elseif riski<riskj % i风险小，i做让步
        minriski = 1;
        for k = 1:TaskNum
            assigntmp = assigni;
            assigntmp(j) = k;
            risktmpj = RiskCalculation(j,assignj,assigntmp,Rlimit,Utility);
            risktmpi = RiskCalculation(i,assigni,assigntmp,Rlimit,Utility);
            if (risktmpj<riskj) && (risktmpi<=minriski)
                newassigni(j) = k;
                newassignj = newassigni;
                minriski = risktmpi;
                riskj = risktmpj;
            end
        end
    else
        tmp = rand;
        if tmp<0.5 % 随机选择i做让步
            newassignj = assignj;
            newassigni = assigni;
            minriski = 1;
            for k = 1:TaskNum
                assigntmp = assigni;
                assigntmp(j) = k;
                risktmpj = RiskCalculation(j,assignj,assigntmp,Rlimit,Utility);
                risktmpi = RiskCalculation(i,assigni,assigntmp,Rlimit,Utility);
                if (risktmpj<riskj) && (risktmpi<=minriski)
                    newassigni(j) = k;
                    newassignj = newassigni;
                    minriski = risktmpi;
                    riskj = risktmpj;
                end
            end
        else % 随机选择j做让步
            minriskj = 1;
            for k = 1:TaskNum
                
                assigntmp = assigni;
                assigntmp(i) = k;
                risktmpj = RiskCalculation(j,assignj,assigntmp,Rlimit,Utility);
                risktmpi = RiskCalculation(i,assigni,assigntmp,Rlimit,Utility);
                if (risktmpi<riski) && (risktmpj<=minriskj)
                    newassignj(i) = k;
                    newassigni = newassignj;
                    minriskj = risktmpj;
                    riski = risktmpi;
                end
                
            end
        end
    end
    
end

end

%% Risk Calculation
function [risk,utility1,utility2] = RiskCalculation(i,assign1,assign2,Rlimit,Utility)

task1 = assign1(i);
task2 = assign2(i);
coagent1 = length(find(assign1==task1));
coagent2 = length(find(assign2==task2));

utility1 = AuxiliaryUtility(i, task1, Rlimit, coagent1, Utility);
utility2 = AuxiliaryUtility(i, task2, Rlimit, coagent2, Utility);

if utility1 == 0
    risk = 1;
else
    risk = (utility1 - utility2)/utility1;
end

end