% ����Lagrange���ӷ���SAP�㷨
function [Ug,assign_plan] = SAPLagrangeFunction(model,Epis,initial_plan)
% ���룺 model:ģ��
%        Epis:��������
% ����� assign_plan�����䷽��

AgentNum = model.Missiles.alive_num;
TargetNum = model.Targets.alive_num;
ProbDist = ones(AgentNum, TargetNum);%������ѡ��������ʷֲ���洢
assign_plan = initial_plan;%randi(model.num_targets,model.num_missiles,1);%���������

if AgentNum == 1
    alive_target = model.Targets.alive_targets;
    Ug = max(model.Targets.value(alive_target));
    assign_plan = 1;
    return
end
if TargetNum == 1
    alive_target = model.Targets.alive_targets;
    Ug = max(model.Targets.value(alive_target));
    assign_plan = ones(AgentNum,1);
    return
end
%% ��������
T = Epis;%�ܵ�������
lambda = max(model.Targets.value);

sap_satisfied = zeros(AgentNum,1);
ex_satisfied = zeros(AgentNum,1);

% for j = 1:TargetNum
%     GlobalUtility(1) = GlobalUtility(1) + TargetUtility(1, j, model, assign_plan, lambda);
% end
attack_flag = MissileCanChangeTarget(model);
Ug = 0;
for k = 2:T
    %for i = 1:AgentNum
        
        i = randi(model.Missiles.alive_num);
        % �����i�������岻ѡ����
        Assign2 = assign_plan;
        Assign2(i) = 0;
        
        tau = 10/k^2;
        % ����k+1ʱ�̲�ͬ�����Ч��ֵ
        AgentUtility = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ�����Ч��ֵ
        
        for j = 1:TargetNum
            %if model.Missiles.target_set(i,j) == true
            if attack_flag(i,j)
                Assign1 = Assign2;
                Assign1(i) = j;
                % �����i��������ı߼�����
                AgentUtility(j) = TargetUtility(i,j,model,Assign1,lambda) - TargetUtility(i,j,model,Assign2,lambda);
            else
                AgentUtility(j) = 0;
            end
        end
        %������ʷֲ����˴����������Ϊԭ��ѡ����ʷֲ�
        AU = AgentUtility/tau;
        newProb = exp(AU)/sum(exp(AU));
        newProb(model.Missiles.target_set(i,:)==false) = 0;
        newProb = newProb/sum(newProb);
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
            Prob = exp(AU)/sum(exp(AU));
            Prob(model.Missiles.target_set(i,:)==0) = 0;
            ProbDist(i,:) = Prob/sum(Prob);
        end
        % �������
        e = rand;
        cp = 0;
        for j = 1:TargetNum
            cp = cp + ProbDist(i,j);
            if e <= cp
                if j==assign_plan(i)
                    sap_satisfied(i) = sap_satisfied(i) + 1;
                else
                    sap_satisfied(i) = 0;
                end
                assign_plan(i) = j;
                break
            end
        end
        
%         if (sum(sap_satisfied) == AgentNum)
%             break;
%         end
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
    
   Ug = GlobalUtility(model,assign_plan);
%     
    if (sum(sap_satisfied>20) == AgentNum)&&(sum(ex_satisfied)==AgentNum)
%         MaxGlobalUtility = GlobalUtility(k);
%         GlobalUtility(k+1:T) = MaxGlobalUtility*ones(1,T-k);
        break
    end
%     if k==T
%         MaxGlobalUtility = GlobalUtility(T);
%     end
end

%-lambda*(sum(RequireNum)-AgentNum);

end

%% ����Ч�ú���
function [U] = TargetUtility(agent, target, model, assign_plan, lambda)

neighbors = find(model.Adjacent(agent,:)==true);
part_missiles = find(assign_plan(neighbors) == target);
num_missiles = length(part_missiles);

if target==0
    global_target = 0;
    num_require = 0;
else
    global_target = model.Targets.alive_targets(target);
    num_require = max(0,model.target_require_num_list(global_target) - model.Targets.num_attacked(global_target));
end
if (num_missiles == 0) || (target==0)
    con = lambda*min(num_require-num_missiles,0);
    U = con;
else
    Time_to_go_max = max(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = model.Targets.value(global_target) - task_cost + con;
end
end

function [U] = Agent_Utility(agent,model,assign_plan,lambda)

Assign = assign_plan;
Assign(agent) = 0;
target = assign_plan(agent);
U = TargetUtility(agent,target,model,assign_plan,lambda)...
    -TargetUtility(agent,target,model,Assign,lambda);
end

%% ������Ч��
function [Ug] = GlobalUtility(model,plan)
Nt = model.Targets.alive_num;
Ut = zeros(Nt,1);
flag = 1;
for j=1:Nt
    part_missiles = find(plan == j);
    num_missiles = length(part_missiles);
    if num_missiles > model.target_require_num_list(j) - model.Targets.num_attacked(j)
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