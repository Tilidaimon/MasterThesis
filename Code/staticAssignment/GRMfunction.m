function [AssignTarget, Ug, time] = GRMfunction(model,Epis)

AgentNum = model.num_missiles;
TargetNum = model.num_targets;


ProbDist = rand(AgentNum, TargetNum);%������ѡ��������ʷֲ���洢
AssignTarget = randi(TargetNum,AgentNum,1);%���������
%% ��������
T = Epis;%�ܵ�������
% Ug = zeros(1,T);
RegretMatrix = zeros(AgentNum,TargetNum);%��ھ���
Ug = zeros(T,1);
Ug(1) = GlobalUtility(model,AssignTarget);
% �ڵ�k������ʱ������ÿ��������ĺ��ֵ����ǰk-1��û��ѡ������l�ĺ�ڳ̶ȣ���ڳ̶�������l������Ч����ߺ���

rho = 0.1;% ��������
alpha = 0.8;% ��������
lambda = max(model.Targets.value);
satisfied = zeros(AgentNum,1);

tic
for k = 2:T
    newAssignTarget = AssignTarget;
    for i = 1:AgentNum
        % �����i�������岻ѡ����
        Assign2 = AssignTarget;
        Assign2(i) = 0;
        Assign1 = Assign2;
        % ����k+1ʱ�̲�ͬ�����Ч��ֵ�ͺ��ֵ
        AgentUtility = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ�����Ч��ֵ
        Regret = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ����ĺ��ֵ
        for j = 1:TargetNum
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,model,Assign1,lambda)...
                - TargetUtility(j,model,Assign2,lambda);% �����i��������ı߼�����
           
        end
        for j = 1:TargetNum
             %����k+1ʱ�̵ĺ��ֵ
            Regret(j) = (1-rho)*RegretMatrix(i,j) + rho*(AgentUtility(j)-AgentUtility(AssignTarget(i)));
        end
        RegretMatrix(i,:) = Regret;
        % ����������i��k+1ʱ��ִ�в�ͬ����ĸ��ʷֲ�
        Regret(Regret<0) = 0;% ���ֵС��0��ֵȡ0
        %RM = Regret/sum(Regret);% ���ƥ�亯��
        RM = zeros(1,TargetNum);
        tau = 10/k^2;
        deno = 0;
        flag = 0;
        for s = 1:TargetNum
            if Regret(s) ~= 0
                deno = deno + exp(Regret(s)/tau);
                
            end
        end
        if isinf(deno)
            flag = 1;
        end
        if flag==1
            
            Rm = Regret/tau;
            maxR = max(abs(Rm));
            modiRM = Rm/maxR*k/2;
            deno = 0;
            for s = 1:TargetNum
                if Regret(s) ~= 0
                    deno = deno + exp(modiRM(s));
                end
            end
            for s = 1:TargetNum
                if Regret(s)>0
                    RM(s) = exp(modiRM(s))/deno;
                else
                    RM(s) = 0;
                end
            end
        else
            for s = 1:TargetNum
                if Regret(s)>0
                    RM(s) = exp(Regret(s)/tau)/deno;
                else
                    RM(s) = 0;
                end
            end
        end
        vertex = zeros(1,TargetNum);
        vertex(AssignTarget(i)) = 1;
        ProbDist(i,:) = alpha*RM + (1-alpha)*vertex;% ��һ��ֵ��Ϊk+1ʱ��ѡ��������ʷֲ�
        
        % �������
        e = rand;
        cp = 0;
        for j = 1:TargetNum
            cp = cp + ProbDist(i,j);
            if e <= cp
                if j==AssignTarget(i)
                    satisfied(i) = satisfied(i) + 1;
                else
                    satisfied(i) = 0;
                end
                newAssignTarget(i) = j;
                
                break
            end
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

%����Ч�ú���
function [U] = TargetUtility(target,model,AssignTarget,lambda)

part_missiles = find(AssignTarget == target);
num_missiles = length(part_missiles);
num_require = model.target_require_num_list(target);
if num_missiles == 0
%     con = lambda*min(num_require-num_missiles,0);
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
        Ut(j) = max(0,model.Targets.value(j)/model.target_require_num_list(j)*num_missiles - cost);
    end
    
end
if flag == 0
    Ug = 0;
else
    Ug = sum(Ut);
end
end