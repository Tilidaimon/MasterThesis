function [zbest,fitnesszbest,yy] = PSOfunction(Val)

[m,n]=size(Val);

%% III. ������ʼ��
c1 = 1.49445;
c2 = 1.49445;
ws = 0.9;
we = 0.4;

maxgen = 50;   % ��������  
sizepop = 100;   %��Ⱥ��ģ

% V�����޵�ѡȡ��
Vmax =n/2;
Vmin = -n/2;
popmax = n;
popmin = 1;
chao_num=10;


%for num=1:100
%% IV. ������ʼ���Ӻ��ٶ�
for i = 1:sizepop
    % �������һ����Ⱥ
    pop(i,:) = randperm(n);    %��ʼ��Ⱥ
    V(i,:) = round(rand(1,n)*3);  %��ʼ���ٶ�
    % ������Ӧ��
    fitness(i) = fun(pop(i,:),Val);   %Ⱦɫ�����Ӧ��
end

%% V. ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex] = max(fitness);
zbest = pop(bestindex,:);   %ȫ�����
gbest = pop;    %�������
fitnessgbest = fitness;   %���������Ӧ��ֵ
fitnesszbest = bestfitness;   %ȫ�������Ӧ��ֵ

%% VI. ����Ѱ��
for i = 1:maxgen
    w(i) = ws - (ws-we)*(i/maxgen);   %���Թ���Ȩ��
    for j = 1:sizepop
        % �ٶȸ���
        %V(j,:) = round(w(i)*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));   %��ɢ��ȡ��
        V(j,:) = round(V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        V(j,find(V(j,:)>Vmax)) = Vmax;
        V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + V(j,:);
        pop(j,find(pop(j,:)>popmax)) = popmax;
        pop(j,find(pop(j,:)<popmin)) = popmin;
        
        %% �ҵ����ظ���������������ٷ���
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
        
        % ��Ӧ��ֵ����
        fitness(j) = fun(pop(j,:),Val); 
    end
    
    for j = 1:sizepop  
        % �������Ÿ���
        if fitness(j) > fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % Ⱥ�����Ÿ���
        if fitness(j) > fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    
%     %% �����Ż�    
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
%         zbest = chaotic(chaobestindex,:);   %ȫ�����
%         fitnesszbest = chaobestfitness;   %ȫ�������Ӧ��ֵ
%     end

%% �����Ż�
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
%           zbest = chaotic(chaobestindex,:);   %ȫ�����
%           fitnesszbest = chaobestfitness;   %ȫ�������Ӧ��ֵ
%       end
%       end
    
    yy(i) = fitnesszbest;            
end

end