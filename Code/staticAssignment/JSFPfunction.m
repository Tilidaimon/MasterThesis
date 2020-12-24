function [AssignTarget, Ug, time] = JSFPfunction(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;

AssignTarget = randi(TargetNum,AgentNum,1);%���������
%% ��������
T = Epis;%�ܵ�������
% Ug = zeros(1,T);
MeanUtility = zeros(AgentNum,TargetNum);%��ھ���
Ug = zeros(T,1);
Ug(1) = GlobalUtility(model,AssignTarget);
% �ڵ�k������ʱ������ÿ��������ĺ��ֵ����ǰk-1��û��ѡ������l�ĺ�ڳ̶ȣ���ڳ̶�������l������Ч����ߺ���

satisfied = zeros(AgentNum,1);
lambda = max(model.Targets.value);
tic
for k=1:T
    newAssignTarget = AssignTarget;
    for i = 1:AgentNum
        % �����i�������岻ѡ����
        Assign2 = AssignTarget;
        Assign2(i) = 0;
        
        % ����k+1ʱ�̲�ͬ�����Ч��ֵ�ͺ��ֵ
        AgentUtility = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ�����Ч��ֵ
        %Regret = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ����ĺ��ֵ
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,model,Assign1,lambda)...
                - TargetUtility(j,model,Assign2,lambda);% �����i��������ı߼�����
            %����k+1ʱ�̵ĺ��ֵ
            MeanUtility(i,j) = ((k-1)*MeanUtility(i,j) + 1*AgentUtility(j))/(k);
        end
        % ����������i��k+1ʱ��ִ�в�ͬ����ĸ��ʷֲ�
        %RM = Regret/sum(Regret);% ���ƥ�亯��
        
        
        % �������
        select_target = find(MeanUtility(i,:)==max(MeanUtility(i,:)));
        prob = rand;
        acp = 0.5;
        if length(select_target)>1
            select_index = randi(length(select_target));
            newassign = select_target(select_index);
        else
            newassign = select_target;
        end
        if prob>acp
            newAssignTarget(i) = newassign;
        end
        if newAssignTarget(i)==AssignTarget(i)
            satisfied(i) = satisfied(i) + 1;
        else
            satisfied(i) = 0;
        end
    end
    AssignTarget = newAssignTarget;
    %������Ч��
    Ug(k) = GlobalUtility(model,AssignTarget);
    if sum(satisfied>20)==AgentNum
        time = toc
        Ug(k+1:T) = Ug(k)*ones(T-k,1);
        break
    end
end
time = toc

end

%% ����Ч�ú���
function [U] = TargetUtility(target,model,AssignTarget,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = model.target_require_num_list(target);
if num_missiles == 0
%     con = lambda*(num_require-num_missiles);
%     task_cost = 0;
    U = 0;
else
    Time_to_go_max = sum(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = model.Targets.value(target)*num_missiles - task_cost + con;
end


end

%% ȫ��Ч��
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