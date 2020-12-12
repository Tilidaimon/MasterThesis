% ����������Ŀ�깥����̬ģ��
% ��������� Chapter 5
% ���ģ�� Stochastic Game Model

classdef DynamicMissileAndTarget
    % ������Ŀ��ģ�� Model of Missiles and Target
    % ��Ҫָ�������� properties
    properties
        num_missiles; % �������� numbers of missiles
        num_targets; % Ŀ������ numbers of targets
        target_require_num_list; % Ŀ�����赼������ required numbers of missiles to attack targets
        order_targets;
        
        % ������ṹ�壬����λ�á��ٶȡ����ٶȡ������
        Missiles;
        Targets;
        
    end
    
    % ����
    properties
        dT = 0.02; % ����ʱ��0.02s
        tau_m = 0.04; % ������Чʱ�䳣��
        tau_t = 0.06; % Ŀ���Чʱ�䳣��
            
        maxVOfMissile  = 1.5; % ��������ٶ�
        maxVOfTarget = 0.3; % Ŀ������ٶ�
        GOfMissile = 0.5; % ���������ٶ�Gֵ
        GOfTarget = 0.03; % Ŀ����ٶ�Gֵ
        
        max_inexcapable_radius = 15; % ��󲻿��������뾶
        communication_radius = 30; % ͨ�Ű뾶
        radar_radius = 60; % �����״�̽��뾶
    end
    
    % ����õ��ĳ���
    properties
        Tmax = 10000; % �ο���󹥻�ʱ�䣬���ʱ���ɵ�һ�ι����õ�
        Time_to_go =[]; % ����Ԥ��ʣ�๥��ʱ��
        extend_time_togo = []; % ��չĿ��Ԥ��ʣ�๥��ʱ��
        Energy_opt = []; % ���ſ���������
        extend_energy = []; % ��չĿ��Ԥ������
        extend_value = [];
        Adjacent = []; % �ڽӾ���
        distance_missiles_targets = []; % ������Ŀ��������
        oig_acc = []; % ���ſ����Ƶ��ɾ���
        reassign_flag = 0; % �Ƿ������ط���
        maxtime = 0;
    end
    
    methods
        %% ���캯��
        function model = DynamicMissileAndTarget(num_missiles, num_targets)
            model.num_missiles = num_missiles;
            model.num_targets = num_targets;
        end
        
        %% �������̬�Ʋ���
        function model = RandomSetSituation(model)
            model.Missiles.p = 10 * rand(model.num_missiles, 2);
%             model.Missiles.p(model.num_missiles/2+1:model.num_missiles,2) = 40 + 10*rand(model.num_missiles/2, 1);
            %model.Missiles.p(:,2) = 20 * rand(model.num_missiles, 1);
            model.Missiles.v = 1 * ones(model.num_missiles, 1);
            model.Missiles.angle = 0.1 * pi * (rand(model.num_missiles, 1)-0.5);
            model.Missiles.alive = true(model.num_missiles,1);
            model.Missiles.alive_missiles = 1:model.num_missiles;
            model.Missiles.alive_num = model.num_missiles;
            model.Missiles.diff_R = -1*ones(model.num_missiles,1);
            
            
            model.Targets.p = 20*rand(model.num_targets, 2);
            model.Targets.p(:,1) = 30 + 5*rand(model.num_targets, 1);
            model.Targets.v = 0.6*ones(model.num_targets, 1);
            model.Targets.angle = 10/180*pi*rand(model.num_targets, 1);
            model.Targets.acc = zeros(model.num_targets, 1);
            model.Targets.alive = true(model.num_targets,1);
            model.Targets.num_attacked = zeros(model.num_targets,1);
            model.Targets.alive_targets = 1:model.num_targets;
            model.Targets.alive_num = model.num_targets;
            
            % ��Ҫ����ö����������Ŀ�����
            index_important_targets = randperm(model.num_targets, model.num_missiles-model.num_targets);
            model.target_require_num_list = ones(model.num_targets,1);
            model.target_require_num_list(index_important_targets) = model.target_require_num_list(index_important_targets) +1;
            
            % ����Ŀ����������
            k=1;
            model.order_targets = zeros(model.num_targets,1);
            for i=1:model.num_targets
                for j=1:model.target_require_num_list(i)
                    model.order_targets(k) = i;
                    k = k+1;
                end
            end
            
            % Ŀ���ֵ
            model.Targets.value = 200*ones(model.num_targets,1);
            model.Targets.value(index_important_targets) = 2*model.Targets.value(index_important_targets);
            model.extend_value = model.Targets.value(model.order_targets);
            
            % ��������ڽӾ���
            model.Missiles.distance_missiles = DistanceMatrixMissiles(model);
            model.Adjacent = GetAdjacentMatrix(model);
            
            % ���㵼����Ŀ��������
            model.distance_missiles_targets = DistanceMatrixMissileAndTarget(model);
            model.Missiles.target_set = DetectTargetSet(model);
            
            % ����Ŀ����Ϣ
            model = TargetShared(model);
            
            % ����ʣ�๥��ʱ�����
            model.Time_to_go = AttackTimeToGo(model);
            model.extend_time_togo = model.Time_to_go(:,model.order_targets);
            model.maxtime = max(max(model.Time_to_go));
            
            % ���������Ƶ��ɼ�Ԥ����������
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model,0);
            model.extend_energy = model.Energy_opt(:,model.order_targets);
            
            
        end
        
        %% �û�ָ��̬�Ʋ���
        function model = SetSituation(model, Missiles, Targets, index_important_targets, target_value)
            model.Missiles = Missiles;
            model.Targets = Targets;
            
            model.target_require_num_list = ones(model.num_targets,1);
            model.target_require_num_list(index_important_targets) = model.target_require_num_list(index_important_targets) +1;
            
            % Ŀ���ֵ
            model.Targets.value = target_value;
            
            % ��������ڽӾ���
            model.Missiles.distance_missiles = DistanceMatrixMissiles(model);
            model.Adjacent = GetAdjacentMatrix(model);
            
            % ���㵼����Ŀ��������
            model.distance_missiles_targets = DistanceMatrixMissileAndTarget(model);
            model.Missiles.target_set = DetectTargetSet(model);
            
            % ����Ŀ����Ϣ
            model = TargetShared(model);
            
            % ����ʣ�๥��ʱ�����
            model.Time_to_go = AttackTimeToGo(model);
            
            % ���������Ƶ��ɼ�Ԥ����������
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model);
            model.maxtime = max(max(model.Time_to_go));
        end
        
        %% ���ݿ�ѡĿ�꼯��ó�ʼ��
        function initial_plan = GetInitialPlan(model)
            initial_plan = zeros(model.Missiles.alive_num,1);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                select_targets = find(model.Missiles.target_set(m,:)==1);
                select_index = randi(length(select_targets));
                initial_plan(i) = select_targets(select_index);
            end
        end
        
        %% �����δ��ɹ����ĵ�������
        function model = GetAliveMissiles(model)
            index = 1:model.num_missiles;
            model.Missiles.alive_missiles = index(model.Missiles.alive==1);
            model.Missiles.alive_num = sum(model.Missiles.alive);
        end
        
        %% �����δ�������Ŀ�꼯��
        function model = GetAliveTarget(model)
            index = 1:model.num_targets;
            model.Targets.alive_targets = index(model.Targets.alive==1);
            model.Targets.alive_num = sum(model.Targets.alive);
        end
        
        %% ȡ�����浼����Ŀ��ķ����Ӽ�
        function sub_plan = GetAliveAssignPlan(model,assign_plan)
            sub_plan = zeros(model.Missiles.alive_num,1);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                if t == 0
                    sub_plan(i) = 0;
                else
                    sub_plan(i) = find(model.Targets.alive_targets==t);
                end
            end
        end
        
        %% �ɷ����Ӽ����Ƶõ�ȫ�����
        function new_assign_plan = GetAllAssignPlan(model,assign_plan,sub_plan)
            new_assign_plan = assign_plan;
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                new_assign_plan(m) = model.Targets.alive_targets(sub_plan(i));
            end
        end
        
        %% ���㵼���ڽӾ���
        function adjacent = GetAdjacentMatrix(model)
            n = model.Missiles.alive_num;
            adjacent = false(n);
            for i=1:n
                for j=1:n
                    if model.Missiles.distance_missiles(i,j) < model.communication_radius
                        adjacent(i,j) = true;
                    end
                end
            end
        end
        
        %% ���㵼��̽�⵽��Ŀ�꼯��
        function target_set = DetectTargetSet(model)
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            target_set = false(Nm, Nt);
            for i=1:Nm
                for j=1:Nt
                    if model.distance_missiles_targets(i,j) < model.radar_radius
                        target_set(i,j) = true;
                    end
                end
            end
        end
        
        %% ��������Ŀ�꼯�� Missiles in same group share targets
        function model = TargetShared(model)
            %target_list = false(model.num_missiles,model.num_targets);
            for i=1:model.Missiles.alive_num
                for j=1:model.Missiles.alive_num
                    if model.Adjacent(i,j) == true
                        model.Missiles.target_set(j,:) = model.Missiles.target_set(j,:) | model.Missiles.target_set(i,:);
                    end
                end
                %target_list(model.Adjacent(i,:)==true,:) = ones(num_neighbors,1) * target;
            end
        end
        
        %% ����ʣ�๥��ʱ�����
        function TimeToGo = AttackTimeToGo(model)
            N = 3; % ����ϵ��
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            TimeToGo = zeros(Nm, Nt);
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                for j=1:Nt
                    t = model.Targets.alive_targets(j);
                    d = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(t))]; % ������������
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % ���߽�����
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % ��λ��[0,pi]
                    aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % Ŀ�귽������
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % �����[0,pi]
%                     relative_theta = (model.Missiles.v(m)*sin(alpha) - model.Targets.v(t)*sin(beta))/d;
                    relative_v = model.Missiles.v(m) - model.Targets.v(t)*cos(beta);
                    TimeToGo(i,j) = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                end
            end
        end
        
        %% ���������Ƶ��ɼ���Ԥ����������
        function [acc, energy] = OptimalInterceptGuidance(model,step)
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            acc = zeros(Nm, Nt);
            energy = zeros(Nm, Nt);
            t_go = model.maxtime - step*model.dT;
            T_m = t_go/model.tau_m;
            T_t = t_go/model.tau_t;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                for j=1:Nt
                    t = model.Targets.alive_targets(j);
                    
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
                    R = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(t))]; % ������������
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % ���߽�����
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % ��λ��[0,pi]
                    aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % Ŀ�귽������
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % �����[0,pi]
%                     dotR = model.Targets.v(t)*cos(beta) - model.Missiles.v(m)*cos(alpha); % ��Ŀ����仯��
                    
                    
                    % ����ǰһʱ��λ��
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % Ŀ��ǰһʱ��λ��
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % ǰһʱ�̵�Ŀ����
                    dotR = (R-pre_R)/model.dT; % ��Ŀ����仯��
                    
                    acc(i,j) = Noig/t_go^2 * (R + dotR*t_go + model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i) - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(j));
%                     acc(i,j) = sign(acc(i,j)) * min(abs(acc(i,j)),model.GOfMissile);
                    
                    P = [1;t_go;model.tau_m^2*(exp(-T_m)-1+T_m)];
                    X_state = [R;dotR;model.Missiles.acc(i)];
                    energy(i,j) = 0.5*X_state'*(P*P')*X_state;
                end
            end
        end
        
        
        
        %% �������������Ƶ���OIG�ƶ�
        function [model,assign_plan] = MissilesMoveByOIG(model,assign_plan,step)
            N=3;
            Nm = model.Missiles.alive_num;
            time_to_go = model.maxtime - step*model.dT;
            T_m = time_to_go/model.tau_m;
            T_t = time_to_go/model.tau_t;
%             Nt = model.Targets.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                %target = model.Targets.alive_targets(t);
                if t ~= 0 % ���������䵽Ŀ��
%                     d = model.distance_missiles_targets(i,model.Targets.alive_targets(t));
                    aM = [cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))]; % ������������
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % ���߽�����
                    R = norm(b);
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % ��λ��[0,pi]
%                     aT = [cos(model.Targets.angle(t)),sin(model.Targets.angle(t))]; % Ŀ�귽������
%                     beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % �����[0,pi]
%                     relative_v = model.Missiles.v(m) - model.Targets.v(t)*cos(beta);
                    

                    % ����ǰһʱ��λ��
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % Ŀ��ǰһʱ��λ��
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % ǰһʱ�̵�Ŀ����
                    diff_R = (R-pre_R)/model.dT; % ��Ŀ����仯��
                    
%                     time_to_go = R/(-diff_R)*(1+0.5*alpha^2/(2*N-1));
                    
                    
                    Noig = T_m^2*(exp(-T_m)-1+T_m)/(-0.5*exp(-2*T_m)-2*T_m*exp(-T_m)+T_m^3/3-T_m^2+T_m+0.5);
%                     TimeToGo = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                    
                    %                     dotR = model.Targets.v(target)*cos(beta) - model.Missiles.v(i)*cos(alpha); % ��Ŀ����仯��
                    acc = Noig/time_to_go^2 * (R + diff_R*time_to_go - model.tau_m^2*(exp(-T_m)-1+T_m)*model.Missiles.acc(i)...
                        - model.tau_t^2*(exp(-T_t)-1+T_t)*model.Targets.acc(t));
                    model.Missiles.acc(i) = sign(acc)*min(model.GOfMissile,abs(acc));
                    if R<0.05 || diff_R>0% �������30m
                        model.Targets.num_attacked(t) = model.Targets.num_attacked(t) + 1;
                        model.Missiles.alive(m) = 0;
                        model.Missiles.v(m) = 0; % ����Ŀ��ĵ���ֹͣ�ƶ�
                        model.Missiles.target_set(i,:) = false;
                        
                        
                        if model.Targets.num_attacked(t) == model.target_require_num_list(t)
                            
                            model.Targets.v(t) = 0; % Ŀ�걻���У�ֹͣ�ƶ�
                            model.Targets.alive(t) = 0;
                            model.Missiles.target_set(:,t) = false; % �����е�Ŀ�겻�ٳ�Ϊ��ѡĿ��
                            if sum(assign_plan==t)>1
                                co_missiles = assign_plan == t;
                                co_m = co_missiles & model.Missiles.alive;
                                assign_plan(co_m==1) = 0;
                                model.reassign_flag = 1;
                            end
                            
                                
                        end
                    else
                        
                        accMX = model.Missiles.acc(m)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(m)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                        model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                        model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                        model.Missiles.p(m,:) = model.Missiles.p(m,:) + model.Missiles.v(m)*model.dT*...
                            [cos(model.Missiles.angle(m)), sin(model.Missiles.angle(m))];
%                         model.Missiles.angle(m) = r;
                    end
                end
            end
        end
        
        %% ������PNG��ʽ�ƶ�
        function [model,assign_plan] = MissilesMoveByPNG(model,assign_plan)
            N=4;
            Nm = model.Missiles.alive_num;
            for i=1:Nm
                m = model.Missiles.alive_missiles(i);
                t = assign_plan(m);
                %t_i = find(model.Targets.alive_targets==t);
                if t ~= 0 % ���������䵽Ŀ��
                    %R = model.distance_missiles_targets(i,t_i);
                    b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % ���߽�����
                    R = norm(b);
                    v = [-model.Missiles.v(m)*cos(model.Missiles.angle(m)) + model.Targets.v(t)*cos(model.Targets.angle(t)),...
                        -model.Missiles.v(m)*sin(model.Missiles.angle(m)) + model.Targets.v(t)*sin(model.Targets.angle(t))];
                    
                    % ����ǰһʱ��λ��
                    pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                    % Ŀ��ǰһʱ��λ��
                    pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
%                   
                    pre_R = norm(pre_missiles_p - pre_target_p); % ǰһʱ�̵�Ŀ����
                    diff_R = (R-pre_R)/model.dT; % ��Ŀ����仯��
                    pre_diff_R = model.Missiles.diff_R(m);
                    model.Missiles.diff_R(m) = diff_R;
%                   
%                     lambda = atan2(b(2),b(1));
%                     pre_lambda = atan2(pre_target_p(2)-pre_missiles_p(2),pre_target_p(1)-pre_missiles_p(1));
%                     Qtm = (lambda-pre_lambda)/model.dT;
                    
                    
                    Vtm = -(b(1)*v(1)+b(2)*v(2))/R;
                    Qtm = (b(1)*v(2)-b(2)*v(1))/(R^2);
                    if Qtm<-pi
                        Qtm = Qtm + 2*pi;
                    elseif Qtm>pi
                        Qtm = Qtm - 2*pi;
                    end
                    acc = N*Vtm*Qtm;
                    model.Missiles.acc(m) = sign(acc)*min(abs(acc),model.GOfMissile);
%                     r = model.Missiles.angle(m) + model.Missiles.acc(m)*model.dT;
                    if R<0.05 || diff_R>0% �������30m
                        model.Targets.num_attacked(t) = model.Targets.num_attacked(t) + 1;
                        model.Missiles.alive(m) = 0;
                        model.Missiles.v(m) = 0; % ����Ŀ��ĵ���ֹͣ�ƶ�
                        model.Missiles.target_set(i,:) = false;
                        
                        
                        if model.Targets.num_attacked(t) == model.target_require_num_list(t)
                            
                            model.Targets.v(t) = 0; % Ŀ�걻���У�ֹͣ�ƶ�
                            model.Targets.alive(t) = 0;
                            model.Missiles.target_set(:,t) = false; % �����е�Ŀ�겻�ٳ�Ϊ��ѡĿ��
                            if sum(assign_plan==t)>1
                                co_missiles = assign_plan == t;
                                co_m = co_missiles & model.Missiles.alive;
                                assign_plan(co_m==1) = 0;
                                model.reassign_flag = 1;
                            end
                            
                                
                        end
                    else
                        
                        accMX = model.Missiles.acc(m)*(-sin(model.Missiles.angle(m)));
                        accMY = model.Missiles.acc(m)*cos(model.Missiles.angle(m));
                        VX = model.Missiles.v(m)*cos(model.Missiles.angle(m)) + accMX*model.dT;
                        VY = model.Missiles.v(m)*sin(model.Missiles.angle(m)) + accMY*model.dT;
                        model.Missiles.angle(m) = atan2(VY,VX);
                        model.Missiles.angle(model.Missiles.angle>pi) = model.Missiles.angle(model.Missiles.angle>pi) - 2*pi;
                        model.Missiles.angle(model.Missiles.angle<-pi) = model.Missiles.angle(model.Missiles.angle<-pi) + 2*pi;
                        model.Missiles.p(m,:) = model.Missiles.p(m,:) + model.Missiles.v(m)*model.dT*...
                            [cos(model.Missiles.angle(m)), sin(model.Missiles.angle(m))];
%                         model.Missiles.angle(m) = r;
                    end
                end
               

            end
        end

        
        %% �жϵ����Ƿ���Ըı�Ŀ��
        function flag = MissileCanChangeTarget(model)
            flag = false(model.Missiles.alive_num,model.Targets.alive_num);
            N=4;
            Nm = model.Missiles.alive_num;
            Nt = model.Targets.alive_num;
            for i=1:Nm
                for j=1:Nt
                    m = model.Missiles.alive_missiles(i);
                    t = model.Targets.alive_targets(j);
                    %t_i = find(model.Targets.alive_targets==t);
                    if t ~= 0 % ���������䵽Ŀ��
                        %R = model.distance_missiles_targets(i,t_i);
                        b = [model.Targets.p(t,1)-model.Missiles.p(m,1), model.Targets.p(t,2)-model.Missiles.p(m,2)]; % ���߽�����
                        R = norm(b);
                        v = [-model.Missiles.v(m)*cos(model.Missiles.angle(m)) + model.Targets.v(t)*cos(model.Targets.angle(t)),...
                            -model.Missiles.v(m)*sin(model.Missiles.angle(m)) + model.Targets.v(t)*sin(model.Targets.angle(t))];
                        
                        % ����ǰһʱ��λ��
                        pre_missiles_p = model.Missiles.p(m,:) - model.Missiles.v(m)*model.dT*[cos(model.Missiles.angle(m)),sin(model.Missiles.angle(m))];
                        % Ŀ��ǰһʱ��λ��
                        pre_target_p = model.Targets.p(t,:) - model.Targets.v(t)*model.dT*[cos(model.Targets.angle(t)),sin(model.Targets.angle(t))];
                        %
                        pre_R = norm(pre_missiles_p - pre_target_p); % ǰһʱ�̵�Ŀ����
                        Vtm = (R-pre_R)/model.dT; % ��Ŀ����仯��
                        %pre_diff_R = model.Missiles.diff_R(m);
                        
                        lambda = atan2(b(2),b(1));
                        pre_lambda = atan2(pre_target_p(2)-pre_missiles_p(2),pre_target_p(1)-pre_missiles_p(1));
                        Qtm = (lambda-pre_lambda)/model.dT;
                        
                        acc = N*Vtm*Qtm;
                        
                        if abs(acc) <= model.GOfMissile
                            flag(i,j) = true & model.Missiles.target_set(i,j);
                        end
                    end
                    
                    %                     Vtm = -(b(1)*v(1)+b(2)*v(2))/R;
                    %                     Qtm = (b(1)*v(2)-b(2)*v(1))/(R^2);
                    
                end
            end
        end
        
        %% Ŀ���ƶ�
        function model = TargetMove(model)
            model.Targets.p = model.Targets.p + ...
                model.dT*[model.Targets.v .* cos(model.Targets.angle), model.Targets.v .* sin(model.Targets.angle)];
        end
        
%         %% ����ÿ�ε����ķ����ʼ��
%         function initial_plan = GetInitialAssignment(model,pre_plan)
%             % pre_planΪǰһʱ�̵ķ���⣬��һ��Ϊ�������
%             for i=1:model.num_missiles
%                 
%         end
        
        %% ���㵼��Ŀ�������� Get the distance matrix
        function dismatrix = DistanceMatrixMissileAndTarget(model)
            dismatrix = zeros(model.Missiles.alive_num, model.Targets.alive_num);
            for i=1:model.Missiles.alive_num
                m = model.Missiles.alive_missiles(i);
                for j=1:model.Targets.alive_num
                    t = model.Targets.alive_targets(j);
                    dismatrix(i,j) = norm(model.Missiles.p(m,:)-model.Targets.p(t,:));
                end
            end
            
        end
        
         %% ���㵼��Ŀ�������� Get the distance matrix
        function dismatrix = DistanceMatrixMissiles(model)
            dismatrix = zeros(model.Missiles.alive_num, model.Targets.alive_num);
            for i=1:model.Missiles.alive_num
                m1 = model.Missiles.alive_missiles(i);
                for j=1:model.Missiles.alive_num
                    m2 = model.Missiles.alive_missiles(j);
                    if i==j
                        dismatrix(i,i) = 0;
                    elseif i<j
                        dismatrix(i,j) = norm(model.Missiles.p(m1,:)-model.Missiles.p(m2,:));
                    else
                        dismatrix(i,j) = dismatrix(j,i);
                    end
                end
            end
            
        end
        
        %% ����
        function plan = DecodeFromExtendPlan(model,extend_plan)
            plan = zeros(model.num_missiles,1);
            for i=1:model.num_missiles
                plan(i) = model.order_targets(extend_plan(i));
            end
        end
        
    end
end



