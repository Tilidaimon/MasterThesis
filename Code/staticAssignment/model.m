clear
clc

%% 导弹与目标位置初始化
Nm = 11;
Nt = 10;
Missile.p = 10*rand(Nm,2);
Missile.v = 0.4*ones(Nm,1);
Missile.angle = 0.1*pi*rand(Nm,1);
Missile.acc = zeros(Nm,1);

Target.p = 20*rand(Nt,2);
Target.p(:,1) = 10 + 5*rand(Nt,1);
Target.v = 0.3*ones(Nt,1);
Target.angle = 10/180*pi*rand(Nt,1);
Target.acc = (-2+4*rand(Nt,1))*0.0098;

important_target_index = randi(Nt);
Target.value = 200*ones(Nt,1);
Target.value(important_target_index) = 300;

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
[SAP_assign, SAP_global, SAP_global_series] = SAPfunction(Target.value,Target.requireNum,T_go,J_opt,Epis);
[SAP_lagrange_assign, SAP_lagrange_global, SAP_lagrange_global_series] = SAP_lagrange_function(Target.value,Target.requireNum,T_go,J_opt,Epis);

episode = 1:Epis;
plot(episode,SAP_global_series,episode,SAP_lagrange_global_series);
legend('SAP with newU','SAP lagrange')
