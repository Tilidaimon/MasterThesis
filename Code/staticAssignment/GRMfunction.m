function [AssignTarget, MaxGlobalUtility, GlobalUtility] = GRMfunction(TargetValue,CompleteProb)

[AgentNum,TargetNum] = size(CompleteProb);


ProbDist = rand(AgentNum, TargetNum);%������ѡ��������ʷֲ���洢
AssignTarget = randperm(AgentNum);%���������
%% ��������
T = 500;%�ܵ�������
GlobalUtility = zeros(1,T);
RegretMatrix = zeros(AgentNum,TargetNum);%��ھ���
for j = 1:TargetNum
        GlobalUtility(1) = GlobalUtility(1) + TargetUtility(j, TargetValue, AssignTarget, CompleteProb);
end
% �ڵ�k������ʱ������ÿ��������ĺ��ֵ����ǰk-1��û��ѡ������l�ĺ�ڳ̶ȣ���ڳ̶�������l������Ч����ߺ���

rho = 0.1;% ��������
alpha = 0.5;% ��������

for k = 2:T
    for i = 1:AgentNum
        % �����i�������岻ѡ����
        if i==1
            Assign2 = [0,AssignTarget(2:AgentNum)];
        elseif i == AgentNum
            Assign2 = [AssignTarget(1:AgentNum-1),0];
        else
            Assign2 = [AssignTarget(1:i-1),0,AssignTarget(i+1:AgentNum)];
        end
        
        % ����k+1ʱ�̲�ͬ�����Ч��ֵ�ͺ��ֵ
        AgentUtility = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ�����Ч��ֵ
        Regret = zeros(1,TargetNum);%������i��k+1ʱ��ִ�в�ͬ����ĺ��ֵ
        for j = 1:TargetNum
            Assign1 = Assign2;
            Assign1(i) = j; 
            AgentUtility(j) = TargetUtility(j,TargetValue,Assign1,CompleteProb)...
                - TargetUtility(j,TargetValue,Assign2,CompleteProb);% �����i��������ı߼�����
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
        for s = 1:TargetNum
            if Regret(s) ~= 0
                deno = deno + exp(Regret(s)/tau);
            end
        end
        for s = 1:TargetNum
            if Regret(s)>0
                RM(s) = exp(Regret(s)/tau)/deno;
            else
                RM(s) = 0;
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
                AssignTarget(i) = j;
                break
            end
        end
    end
    
    %������Ч��
    for j = 1:TargetNum
        GlobalUtility(k) = GlobalUtility(k) + TargetUtility(j, TargetValue, AssignTarget, CompleteProb);
    end
end

MaxGlobalUtility = GlobalUtility(T);

end

%����Ч�ú���
function [U] = TargetUtility(target,targetValue,AssignTarget,completeProb)

Nmax = 1;
totalProb = 1;
% Agents = find(AssignTarget==target);
% nAgent = length(Agents);
[~,TargetNum] = size(AssignTarget);
nAgent = zeros(1,TargetNum);
for i = 1:TargetNum
    nAgent(i) = length(find(AssignTarget==i));
end
if nAgent(target)>0
    Agents = find(AssignTarget == target);
    for i = 1:nAgent(target)
        totalProb = totalProb * (1-completeProb(Agents(i),target));
    end
else
    totalProb = 1;
end
con = sum(min(Nmax-nAgent,0));
U = targetValue(target)*(1-totalProb)+con;

end