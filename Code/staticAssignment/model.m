clear
clc

%% 导弹与目标位置初始化
Nm = 11;
Nt = 10;
Missile.p = 10*rand(Nm,2);
Missile.v = 0.4*ones(Nm,1);
Missile.angle = 0.1*pi*rand(Nm,1);
Missile.acc = zeros(Nm,1);
Missile.adjacent = eye(Nm)+[zeros(1,Nm-1),1;eye(Nm-1),zeros(Nm-1,1)];

Target.p = 20*rand(Nt,2);
Target.p(:,1) = 10 + 5*rand(Nt,1);
Target.v = 0.3*ones(Nt,1);
Target.angle = 10/180*pi*rand(Nt,1);
Target.acc = (-2+4*rand(Nt,1))*0.0098;

important_target_index = randi(Nt);
Target.value = 200*ones(Nt,1);
% Target.value(important_target_index) = 300;

Target.requireNum = ones(Nt,1);
Target.requireNum(important_target_index) = 2;

%% 剩余攻击时间估计
N = 4;% 比例系数

%计算弹目距离
Distance = zeros(Nm,Nt);
for i=1:Nm
    for j=1:Nt
        Distance(i,j) = norm(Missile.p(i,:)-Target.p(j,:));
    end
end


T_go = zeros(Nm,Nt);
for i=1:Nm
    for j=1:Nt
        d = Distance(i,j);
        aM = [cos(Missile.angle(i)),sin(Missile.angle(i))]; % 导弹方向向量
        b = [Target.p(j,1)-Missile.p(i,1),Target.p(j,2)-Missile.p(i,2)]; % 视线角向量
        alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
        aT = [cos(Target.angle(j)),sin(Target.angle(j))]; % 目标方向向量
        beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
        re_v = Missile.v(i) - Target.v(j)*cos(beta);
        T_go(i,j) = d/re_v*(1+0.5*alpha^2/(2*N-1));
    end
end


%% 计算最优能量指标
A_c = zeros(Nm,Nt); % 制导指令
A = zeros(Nm,Nt);
J_opt = zeros(Nm,Nt);
tau = 0.4; % 导弹等效时间常数
for i=1:Nm
    for j=1:Nt
        t_go = T_go(i,j);
        T = t_go/tau;
        Nopl = T^2*(exp(-T)-1+T)/(-0.5*exp(-2*T)-2*T*exp(-T)+T^3/3-T^2+T+0.5);
        R = norm([Target.p(j,1)-Missile.p(i,1),Target.p(j,2)-Missile.p(i,2)]);
        aM = [cos(Missile.angle(i)),sin(Missile.angle(i))]; % 导弹方向向量
        b = [Target.p(j,1)-Missile.p(i,1),Target.p(j,2)-Missile.p(i,2)]; % 视线角向量
        alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % 方位角[0,pi]
        aT = [cos(Target.angle(j)),sin(Target.angle(j))]; % 目标方向向量
        beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % 进入角[0,pi]
        dotR = Target.v(j)*cos(beta) - Missile.v(i)*cos(alpha);
        A_c(i,j) = Nopl/t_go^2 * (R+dotR*t_go - tau^2*(exp(-T)-1+T)*Missile.acc(i)-...
            t_go^2*Target.acc(j)/2);
        
        Npng = 4;
        dotlambda = (Target.v(j)*sin(beta) - Missile.v(i)*sin(alpha))/R;
        A(i,j) = Npng*dotlambda*Missile.v(i);
        P = [1;t_go;tau^2*(exp(-T)-1+T)];
        X_state = [R;dotR;A_c(i,j)];
        J_opt(i,j) = 0.5*X_state'*(P*P')*X_state;
    end
end

%% 分配
Epis = 500;
tic
for r = 1:100
[SAP_lagrange2_assign, SAP_lagrange2_global, SAP_lagrange2_global_series] = SAP_lagrange_function2(Nm,Nt,Missile.adjacent,Target.value,Target.requireNum,T_go,J_opt,Epis);
SAP_lagrange2_time = toc;
tic
[HCG_assign, HCG_global, HCG_global_series] = HCGFunction(Nm,Nt,Missile.adjacent,Target.value,Target.requireNum,T_go,J_opt,Epis);
HCG_time = toc;
tic
[HCG2_assign, HCG2_global, HCG2_global_series] = HCGFunction2(Nm,Nt,Missile.adjacent,Target.value,Target.requireNum,T_go,J_opt,Epis);
HCG_time = toc;
tic
[HCGSAP_assign, HCGSAP_global, HCGSAP_global_series] = HCGFunctionSAP(Nm,Nt,Missile.adjacent,Target.value,Target.requireNum,T_go,J_opt,Epis);
HCG2_time = toc;
% [HCGZ_assign, HCGZ_global, HCGZ_global_series] = HCGFunction_with_Zeuthen(Nm,Nt,Missile.adjacent,Target.value,Target.requireNum,Target.requireNum,T_go,J_opt,Epis);

Ug_SAP_lagrange2(r) = GlobalUtility(SAP_lagrange2_assign,Target.value,T_go,J_opt);
Ug_HCG(r) = GlobalUtility(HCG_assign,Target.value,T_go,J_opt);
% Ug_HCG2(r) = GlobalUtility(HCG2_assign,Target.value,T_go,J_opt);
Ug_HCGSAP(r) = GlobalUtility(HCGSAP_assign,Target.value,T_go,J_opt);

% episode = 1:Epis;
% plot(episode,SAP_lagrange2_global_series,episode,HCG_global_series,episode,HCG2_global_series,episode,HCGSAP_global_series)%,episode,HCGZ_global_series);
% legend('SAP lagrange','HCG','HCG2','HCGSAP')%,'HCGZ')

% % Ug_SAP(r) = GlobalUtility(SAP_assign,Target.value,T_go,J_opt);
% Ug_SAP_lagrange(r) = GlobalUtility(SAP_lagrange_assign,Target.value,T_go,J_opt);
% Ug_SAP_lagrange2(r) = GlobalUtility(SAP_lagrange2_assign,Target.value,T_go,J_opt);
% Ug_HCG(r) = GlobalUtility(HCG_assign,Target.value,T_go,J_opt);
% Ug_HCG2(r) = GlobalUtility(HCG2_assign,Target.value,T_go,J_opt);
% 
% % [Time_SAP(r),J_SAP(r)] = Time_and_Energy(SAP_assign,T_go,J_opt);
% [Time_SAP_lagrange(r),J_SAP_lagrange(r)] = Time_and_Energy(SAP_lagrange_assign,T_go,J_opt);
% [Time_SAP_lagrange2(r),J_SAP_lagrange2(r)] = Time_and_Energy(SAP_lagrange2_assign,T_go,J_opt);
% [Time_HCG(r),J_HCG(r)] = Time_and_Energy(HCG_assign,T_go,J_opt);
% [Time_HCG2(r),J_HCG2(r)] = Time_and_Energy(HCG2_assign,T_go,J_opt);
% % 
end
% 
figure(1)
boxplot([Ug_SAP_lagrange2;Ug_HCG;Ug_HCGSAP]');
tick = {'SAPlagrange2','HCG','HCGSAP'};
set(gca,'XTickLabel',tick);
% 
% figure(2)
% boxplot([Time_SAP_lagrange;Time_SAP_lagrange2;Time_HCG;Time_HCG2]');
% tick = {'SAPlagrange','SAPlagrange2','HCG','HCG2'};
% set(gca,'XTickLabel',tick);
% 
% figure(3)
% boxplot([J_SAP_lagrange;J_SAP_lagrange2;J_HCG;J_HCG2]');
% tick = {'SAPlagrange','SAPlagrange2','HCG','HCG2'};
% set(gca,'XTickLabel',tick);
%% 分配结果效用
function [Ug] = GlobalUtility(Assign,Value,Time_to_go,J_opt)
[Nt,~] = size(Value);
Ug = 0;
Ut = zeros(Nt,1);

for j=1:Nt
    part_missiles = find(Assign == j);
    num_missiles = length(part_missiles);
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = max(Time_to_go(part_missiles,j));
        J_sum = sum(J_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,Value(j) - cost);
    end
    
end
Ug = sum(Ut);
end

function [Time,J] = Time_and_Energy(Assign,T_go,J_opt)

[Nm,~] = size(Assign);
for i=1:Nm
    time_assign(i) = T_go(i,Assign(i));
    J_assign(i) = J_opt(i,Assign(i));
end

Time = max(time_assign);
J = sum(J_assign);
end