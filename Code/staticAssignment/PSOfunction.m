function [zbest,yy,time] = PSOfunction(model)

n = length(model.order_targets);

%% III. 参数初始化
c1 = 1.49445;
c2 = 1.49445;
ws = 0.9;
we = 0.4;

maxgen = 50;   % 进化次数  
sizepop = 100;   %种群规模

% V上下限的选取？
Vmax =round(n/2);
Vmin = round(-n/2);
popmax = n;
popmin = 1;
chao_num=10;


%for num=1:100
%% IV. 产生初始粒子和速度
for i = 1:sizepop
    % 随机产生一个种群
    pop(i,:) = randperm(n);    %初始种群
    V(i,:) = round(rand(1,n)*3);  %初始化速度
    % 计算适应度
    fitness(i) = GlobalUtility(model,pop(i,:));   %染色体的适应度
end

%% V. 个体极值和群体极值
[bestfitness bestindex] = max(fitness);
zbest = pop(bestindex,:);   %全局最佳
gbest = pop;    %个体最佳
fitnessgbest = fitness;   %个体最佳适应度值
fitnesszbest = bestfitness;   %全局最佳适应度值

%% VI. 迭代寻优
tic
for i = 1:maxgen
    w(i) = ws - (ws-we)*(i/maxgen);   %线性惯性权重
    for j = 1:sizepop
        % 速度更新
        %V(j,:) = round(w(i)*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));   %离散化取整
        V(j,:) = round(V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        V(j,find(V(j,:)>Vmax)) = Vmax;
        V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % 种群更新
        pop(j,:) = pop(j,:) + V(j,:);
        pop(j,find(pop(j,:)>popmax)) = popmax;
        pop(j,find(pop(j,:)<popmin)) = popmin;
        
        %% 找到被重复分配的任务，重新再分配
        [B, I] = unique(pop(j,:), 'first');
        repeat=setdiff(1:numel(pop(j,:)), I);
        [m_repeat,n_repeat]=size(repeat);
        set=1:n;
        temp_index=~ismember(set,pop(j,:));
        temp=set(temp_index);
        ass=randperm(n_repeat);
        for k=1:n_repeat
            pop(j,repeat(ass(k)))=temp(k);
        end
        
        % 适应度值更新
        fitness(j) = GlobalUtility(model,pop(j,:)); 
    end
    
    for j = 1:sizepop  
        % 个体最优更新
        if fitness(j) > fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % 群体最优更新
        if fitness(j) > fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    
%     %% 混沌优化    
%     if (i>0.9*maxgen)    
%         zbest_temp=(zbest-1)/(n-1);
%         chaotic_seq(1,:)=zbest_temp;
%         for c=2:chao_num
%             for j=1:n
%                 chaotic_seq(c,j)=4*chaotic_seq(c-1,j)*(1-chaotic_seq(c-1,j));
%             end
%         end
%         
%         chaotic=round(1+(n-1)* chaotic_seq);
%         for j=1:chao_num
%             [B, I] = unique(chaotic(j,:), 'first');
%             repeat=setdiff(1:numel(chaotic(j,:)), I);
%             [m_repeat,n_repeat]=size(repeat);
%             set=1:n;
%             temp_index=~ismember(set,chaotic(j,:));
%             temp=set(temp_index);
%             ass=randperm(n_repeat);
%             for k=1:n_repeat
%                 chaotic(j,repeat(ass(k)))=temp(k);
%             end
%             chao_fitness(j)=fun(chaotic(j,:),Val);
%         end
%         [chaobestfitness chaobestindex] = max(chao_fitness);
%         zbest = chaotic(chaobestindex,:);   %全局最佳
%         fitnesszbest = chaobestfitness;   %全局最佳适应度值
%     end

%% 混沌优化
% if i>0.8*maxgen
%       chaotic_seq(1,:)=rand(1,n);
%       for c=2:chao_num
%             for j=1:n
%                 chaotic_seq(c,j)=4*chaotic_seq(c-1,j)*(1-chaotic_seq(c-1,j));
%             end
%       end
%       
%       for j=1:chao_num
%           chaotic(j,:)=round(zbest+n*(2*chaotic_seq(j,:)-1)/2);
%           chaotic(j,find(chaotic(j,:)>popmax)) = popmax;
%           chaotic(j,find(chaotic(j,:)<popmin)) = popmin;
%           
%           [B, I] = unique(chaotic(j,:), 'first');
%           repeat=setdiff(1:numel(chaotic(j,:)), I);
%           [m_repeat,n_repeat]=size(repeat);
%           set=1:n;
%           temp_index=~ismember(set,chaotic(j,:));
%           temp=set(temp_index);
%           ass=randperm(n_repeat);
%           for k=1:n_repeat
%               chaotic(j,repeat(ass(k)))=temp(k);
%           end
%           chao_fitness(j)=fun(chaotic(j,:),Val);
%       end
%       [chaobestfitness chaobestindex] = max(chao_fitness);
%       if chaobestfitness>fitnesszbest
%           zbest = chaotic(chaobestindex,:);   %全局最佳
%           fitnesszbest = chaobestfitness;   %全局最佳适应度值
%       end
%       end
    
    yy(i) = GlobalUtility(model,zbest);            
end
time = toc
end

function [Ug] = GlobalUtility(model,plan)
Nm = model.num_missiles;
Nt = model.num_targets;
Ut = zeros(Nt,1);

assign = zeros(model.num_missiles,1);
for i=1:Nm
    assign(i) = model.order_targets(plan(i));
end

for j=1:Nt
    part_missiles = find(assign == j);
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


