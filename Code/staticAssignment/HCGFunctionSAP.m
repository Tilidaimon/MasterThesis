function [FinalAssignList, globalutility, time] = HCGFunctionSAP(epsilon, model, Epis)

%% Parameter Initialisation
AgentNum = model.num_missiles;
TargetNum = model.num_targets;
Adjacent = model.Adjacent;

AssignList = randi(TargetNum,AgentNum,1); % Assign List for each agent
% AssignList = randi(TaskNum, AgentNum);
% AssignList = randi(TaskNum,AgentNum,1)*ones(1,AgentNum);

coAgentNums = zeros(TargetNum,1); % Initial Partition
for i = 1:TargetNum
    coAgentNums(i) = length(find(AssignList==i));
end
% coAgentNums(:,1) = AgentNum*ones(AgentNum,1); % Initial numbers distribution of cooperated agents


%% Decision-making process begins
T = Epis; 
hcg_satisfied = false(AgentNum,1);
ex_satisfied = false(AgentNum,1);

globalutility = zeros(1,T);
tic
for t = 1:T
    i = randi(AgentNum);
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
    % update task
    AU = tmpUti/tau;
    deno = 0;
    newProb = zeros(TargetNum,1);
    for k=1:TargetNum
        if AU(k) ~= 0
            deno = deno + exp(AU(k));
        end
    end
    for k=1:TargetNum
        if AU(k)~=0
            newProb(k) = exp(AU(k))/deno;
        else
            newProb(k) = 0;
        end
    end
    flag = 0;
    for j = 1:TargetNum
        if isnan(newProb(j))
            flag = 1;
            break;
        end
    end
    if flag == 0
        ProbDist = newProb;
    else
        AU = AU/max(abs(AU))*t/2;
        deno = 0;
        for k=1:TargetNum
            if AU(k) ~= 0
                deno = deno + exp(AU(k));
            end
        end
        for k=1:TargetNum
            if AU(k)~=0
                newProb(k) = exp(AU(k))/deno;
            else
                newProb(k) = 0;
            end
        end
        ProbDist = newProb;
    end
    e = rand;
    cp = 0;
    for j = 1:TargetNum
        cp = cp + ProbDist(j);
        if e <= cp
            curTask = AssignList(i);
            coAgentNums(curTask) = coAgentNums(curTask) - 1;
            AssignList(i) = j;
            coAgentNums(j) = coAgentNums(j) + 1;
            hcg_satisfied(i) = 1;
            break
        end
    end

    neighbors = find(Adjacent(i,:)==1);
    num_neighbors = length(neighbors);
    if sum(hcg_satisfied(neighbors))==num_neighbors
        for p=1:num_neighbors
            if p==i
                ex_satisfied(p) = 1;
                continue
            end
            ex_satisfied(neighbors) = 1;
            AssignEx1 = AssignList;
            AssignEx2 = AssignList;
            task1 = AssignEx2(i);
            AssignEx2(i) = AssignEx2(p);
            AssignEx2(p) = task1;
            AgentU1 = AuxiliaryUtility(i,AssignEx2(i),coAgentNums(AssignEx2(i)),model,e)-AuxiliaryUtility(i,AssignEx1(i),coAgentNums(AssignEx1(i)),model,e);
            AgentU2 = AuxiliaryUtility(p,AssignEx2(p),coAgentNums(AssignEx2(p)),model,e)-AuxiliaryUtility(p,AssignEx1(p),coAgentNums(AssignEx1(p)),model,e);
            Net_utility = AgentU1+AgentU2;
            if Net_utility>0
                extask = AssignList(i);
                AssignList(i) = AssignList(p);
                AssignList(p) = extask;
                hcg_satisfied(i) = 0;
                hcg_satisfied(p) = 0;
                ex_satisfied(i) = 0;
                ex_satisfied(p) = 0;
            end
            
        end
    end
    
    globalutility(t) = GlobalUtility(model,AssignList);
    if(sum(hcg_satisfied==20)==model.num_missiles && sum(ex_satisfied)==model.num_missiles)
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

% e = 3;
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
