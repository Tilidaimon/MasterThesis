% ����������Ŀ�깥����̬ģ��
% ��������� Chapter 5
% ���ģ�� Stochastic Game Model

classdef MissileAndTarget
    % ������Ŀ��ģ�� Model of Missiles and Target
    % ��Ҫָ�������� properties
    properties
        num_missiles; % �������� numbers of missiles
        num_targets; % Ŀ������ numbers of targets
        target_require_num_list; % Ŀ�����赼������ required numbers of missiles to attack targets
        assign_plan; % ������
        
        % ������ṹ�壬����λ�á��ٶȡ����ٶȡ������
        Missiles;
        Targets;
    end
    
    % ����
    properties
        dT = 0.02; % ����ʱ��0.02s
        maxVOfMissile  = 1.5; % ��������ٶ�
        maxVOfTarget = 0.3; % Ŀ������ٶ�
        GOfMissile = 0.06 * 9.8; % �����ļ��ٶ�Gֵ
        GOfTarget = 0.03 * 9.8; % �������ٶ�Gֵ
        
        max_inexcapable_radius = 25; % ��󲻿��������뾶
        communication_radius = 30; % ͨ�Ű뾶
        radar_radius = 50; % �����״�̽��뾶
    end
    
    % ����õ��ĳ���
    properties
        Tmax = 10000; % �ο���󹥻�ʱ�䣬���ʱ���ɵ�һ�ι����õ�
        Time_to_go =[]; % ����Ԥ��ʣ�๥��ʱ��
        Energy_opt = []; % ���ſ���������
        Adjacent = []; % �ڽӾ���
        distance_missiles_targets = []; % ������Ŀ��������
        oig_acc = []; % ���ſ����Ƶ��ɾ���
    end
    
    methods
        %% ���캯��
        function model = MissileAndTarget(num_missiles, num_targets)
            model.num_missiles = num_missiles;
            model.num_targets = num_targets;
        end
        
        %% �������̬�Ʋ���
        function model = RandomSet(model)
            model.Missiles.p = 100 * rand(model.num_missiles, 2);
            model.Missiles.v = 0.4 * rand(model.num_missiles, 1);
            model.Missiles.angle = 0.1 * pi * rand(model.num_missiles, 1);
            
            
            model.Targets.p = 20*rand(model.num_targets, 2);
            model.Targets.p(:,1) = 60 + 5*rand(model.num_targets, 1);
            model.Targets.v = 0.3*ones(model.num_targets, 1);
            model.Targets.angle = 10/180*pi*rand(model.num_targets, 1);
            
            index_important_targets = randperm(model.)
            model.Targets.value = 
            
            % ��������ڽӾ���
            model.Missiles.distance_missiles = DistanceMatrix(model.Missiles.p, model.Missiles.p);
            model.Adjacent = false(model.num_missiles);
            for i=1:model.num_missiles
                for j=1:model.num_missiles
                    if model.Missiles.distance_missiles(i,j) < model.communication_radius
                        model.Adjacent(i,j) = true;
                    end
                end
            end
            
            % ���㵼����Ŀ��������
            model.distance_missiles_targets = DistanceMatrix(model.Missiles.p, model.Targets.p);
            model.Missiles.target_set = false(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    if model.distance_missiles_targets(i,j) < model.radar_radius
                        model.Missiles.target_set(i,j) = true;
                    end
                end
            end
            
            % ����Ŀ����Ϣ
            model.Missiles.target_set = TargetShared(model.Missiles.target_set, model.Adjacent);
            
            % ����ʣ�๥��ʱ�����
            model.Time_to_go = AttackTimeToGo(model);
            
            % ���������Ƶ��ɼ�Ԥ����������
            model.Missiles.acc = zeros(model.num_missiles, 1);
            [model.oig_acc, model.Energy_opt] = OptimalInterceptGuidance(model);
            
        end
        
        %% ��������Ŀ�꼯�� Missiles in same group share targets
        function target_list = TargetShared(model)
            target_list = false(model.num_missiles);
            for i=1:numMissiles
                target = model.Missiles.target_set(i,:);
                num_neighbors = 0;
                for j=1:numMissiles
                    if model.Adjacent(i,j) == true
                        target = target | model.Missiles.target_set(j,:);
                        num_neighbors = num_neighbors+1;
                    end
                end
                target_list(model.Adjacent(i,:)==true,:) = ones(num_neighbors,1) * target;
            end
        end
        
        %% ����ʣ�๥��ʱ�����
        function TimeToGo = AttackTimeToGo(model)
            N = 4; % ����ϵ��
            TimeToGo = zeros(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    d = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(i)),sin(model.Missiles.angle(i))]; % ������������
                    b = [model.Targets.p(j,1)-model.Missiles.p(i,1), model.Targets.p(j,2)-model.Missiles.p(i,2)]; % ���߽�����
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % ��λ��[0,pi]
                    aT = [cos(model.Targets.angle(j)),sin(model.Targets.angle(j))]; % Ŀ�귽������
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % �����[0,pi]
                    relative_v = model.Missiles.v(i) - model.Targets.v(j)*cos(beta);
                    model.Time_to_go(i,j) = d/relative_v*(1+0.5*alpha^2/(2*N-1));
                end
            end
        end
        
        %% ���������Ƶ��ɼ���Ԥ����������
        function [acc, energy] = OptimalInterceptGuidance(model)
            tau = 0.4; % ������Чʱ�䳣��
            acc = zeros(model.num_missiles, model.num_targets);
            energy = zeros(model.num_missiles, model.num_targets);
            for i=1:model.num_missiles
                for j=1:model.num_targets
                    t_go = model.Time_to_go(i,j);
                    T = t_go/tau;
                    Noig = T^2*(exp(-T)-1+T)/(-0.5*exp(-2*T)-2*T*exp(-T)+T^3/3-T^2+T+0.5);
                    R = model.distance_missiles_targets(i,j);
                    aM = [cos(model.Missiles.angle(i)),sin(model.Missiles.angle(i))]; % ������������
                    b = [model.Targets.p(j,1)-model.Missiles.p(i,1), model.Targets.p(j,2)-model.Missiles.p(i,2)]; % ���߽�����
                    alpha = acos(max(-1,min(1,dot(aM,b)/(norm(aM)*norm(b)))));  % ��λ��[0,pi]
                    aT = [cos(model.Targets.angle(j)),sin(model.Targets.angle(j))]; % Ŀ�귽������
                    beta = acos(max(-1,min(1,dot(aT,b)/(norm(aT)*norm(b))))); % �����[0,pi]
                    dotR = model.Targets.v(j)*cos(beta) - model.Missiles.v(i)*cos(alpha); % ��Ŀ����仯��
                    acc(i,j) = Noig/t_go^2 * (R+dotR*t_go - tau^2*(exp(-T)-1+T)*model.Missiles.acc(i)-t_go^2*Target.acc(j)/2);
                    
                    P = [1;t_go;tau^2*(exp(-T)-1+T)];
                    X_state = [R;dotR;A_c(i,j)];
                    energy(i,j) = 0.5*X_state'*(P*P')*X_state;
                end
            end
        end
        
        
        
    end
end

%% ���������� Get the distance matrix
function [dismatrix] = DistanceMatrix(pos1, pos2)

[numMissiles,~] = size(pos1);
[numTargets,~] = size(pos2);

dismatrix = zeros(numMissiles, numTargets);
for i=1:pos1
    for j=1:pos2
        dismatrix(i,j) = norm(pos1(i,:)-pos2(j,:));
    end
end

end

