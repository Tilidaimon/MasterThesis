function [AssignTarget, Ug, time] = JSFPfunction(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;


AssignTarget = randi(TargetNum,AgentNum,1);%���������
%% ��������
T = Epis;%�ܵ�������
% Ug = zeros(1,T);
MeanUtility = zeros(AgentNum,TargetNum);%��ھ���
Ug(1) = GlobalUtility(model,AssignTarget);
% �ڵ�k������ʱ������ÿ��������ĺ��ֵ����ǰk-1��û��ѡ������l�ĺ�ڳ̶ȣ���ڳ̶�������l������Ч����ߺ���

satisfied = false(AgentNum,1);
lambda = 100;
tic
for k = 2:T
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
        if length(select_target)>1
            select_index = randi(length(select_target));
            AssignTarget(i) = select_target(select_index);
        else
            AssignTarget(i) = select_target;
        end
        if j==AssignTarget(i)
            satisfied(i) = 1;
        else
            satisfied(i) = 0;
        end
    end
    
    %������Ч��
    Ug(k) = GlobalUtility(model,AssignTarget);
%     if sum(satisfied)==AgentNum
%         break
%     end
end
time = toc

end

%����Ч�ú���
function [U] = TargetUtility(target,model,AssignTarget,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = model.target_require_num_list(target);
if num_missiles == 0
    con = lambda*min(num_require-num_missiles,0);
    U = con;
else
    Time_to_go_max = max(model.Time_to_go(part_missiles,target));
    %J_max_task = max(J_opt(part_missiles,j));
    task_cost = Time_to_go_max+sum(model.Energy_opt(part_missiles,target));
    con = lambda*min(num_require-num_missiles,0);
    U = model.Targets.value(target) - task_cost + con;
end

end

%% ȫ��Ч��
function [Ug] = GlobalUtility(model,plan)
Nt = model.num_targets;
Ut = zeros(Nt,1);

for j=1:Nt
    part_missiles = find(plan == j);
    num_missiles = length(part_missiles);
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = max(model.Time_to_go(part_missiles,j));
        J_sum = sum(model.Energy_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,model.Targets.value(j) - cost);
    end
    
end
Ug = sum(Ut);
end