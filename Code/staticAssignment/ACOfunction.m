function [Plan,max_P,Pseries] = ACOfunction(Val)
%% ��ʼ������
temp_Val=max(max(Val))+1-Val;
[p,q]=size(temp_Val);                     % ���ƾ������
n=max(p,q);
m = 100;                              % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                            % ����������Ҫ�̶�����
rho = 0.1;                           % ��Ϣ�ػӷ�����
Q = 1;                               % ��ϵ��
Eta = 1./temp_Val;                          % ��������
Tau = ones(n,n)./(n*sum(sum(temp_Val))');   % ��ʼ����Ϣ�ؾ���
Table = cell(m,n);                   % ·����¼��
iter = 1;                            % ����������ֵ
iter_max = 50;                      % ���������� 
Route_best = cell(iter_max,n);      % �������·��       
P_best = zeros(iter_max,1);     % �������·���ĳ���  
P_ave = zeros(iter_max,1);      % ����·����ƽ������  


pthreshold=1;

while iter <= iter_max
     pthreshold=pthreshold/log(1+iter);
     % ��������������ϵ�������
      start = cell(m,1);
      for i = 1:m
          temp1 = randperm(p);
          temp2 = randperm(q);
          start{i,1} = [temp1(1),temp2(1)];  
          Table{i,1}=start{i,1};
      end
 
      k=1;
      for i=1:n
          for j=1:n
              allocation_tabu_index{1,k}=[i,j];
              k=k+1;
          end
      end
      
      allocation_number=1:n*n;
      % �������·��ѡ��
      for i = 1:m
          % ����������ѡ��
         for j = 2:n
            %% ���ѷ��ʵĽڵ㼰������ж�������ɱ�
             tabu=[];
             for p=1:j-1
                 tabu=[tabu Table{i,p}'];
             end        
             tabu_index=[];
             for p=1:j-1
                a=(tabu(1,p)-1)*n+1:(tabu(1,p)-1)*n+n;
                tabu_index=[tabu_index a];
             end
             
             for p=1:j-1
                 a0=1:n;
                 a=tabu(2,p)+(a0-1)*n;
                 tabu_index=[tabu_index a];
             end
             tabu_index;
             allow_tabu_index = ~ismember(allocation_number,tabu_index);
             allow = allocation_number(allow_tabu_index);  % �����ʵ���������������
             P = allow;
             last_index=(tabu(1,j-1)-1)*n+tabu(2,j-1);
             % ������м�ת�Ƹ���
             for k = 1:length(allow)
                row=fix((allow(k)-1)/n)+1;
                column=allow(k)-(row-1)*n;
                 P(k) = Tau(row,column)^alpha * Eta(row,column)^beta;
             end
%              P = P/sum(P);
%             %%���̶ķ�ѡ����һ�����ʳ���
%             Pc = cumsum(P);     
%             target_tabu_index = find(Pc >= rand); 
%             target = allow(target_tabu_index(1));
%             target_row=fix((target-1)/n)+1;
%             target_column=target-(target_row-1)*n;
%             Table{i,j} = [target_row,target_column];
            
            if rand<pthreshold
                P = P/sum(P);
                %���̶ķ�ѡ����һ�����ʳ���
                Pc = cumsum(P);
                target_tabu_index = find(Pc >= rand);
                target = allow(target_tabu_index(1));
                target_row=fix((target-1)/n)+1;
                target_column=target-(target_row-1)*n;
                Table{i,j} = [target_row,target_column];
            else
                [~,Pindex]=max(P);
                target = allow(Pindex);
                target_row=fix((target-1)/n)+1;
                target_column=target-(target_row-1)*n;
                Table{i,j} = [target_row,target_column];
            end
            
         end
      end
      % ����������ϵ��������ȫ�̵�����ֵ
      totalP = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          matrixRoute=cell2mat(Route);
          for j = 1:n
              totalP(i) = totalP(i)+temp_Val(matrixRoute(j*2-1),matrixRoute(j*2));
          end
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_P,min_tabu_index] = min(totalP);
          P_best(iter) = min_P;  
          P_ave(iter) = mean(totalP);
          Route_best(iter,:) = Table(min_tabu_index,:);
      else
          [min_P,min_tabu_index] = min(totalP);
          P_best(iter) = min(P_best(iter - 1),min_P);
          P_ave(iter) = mean(totalP);
          if P_best(iter) == min_P
              Route_best(iter,:) = Table(min_tabu_index,:);
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      % ������Ϣ��
      Delta_Tau = zeros(n,n);
      % ������ϼ���
      for i = 1:m
          % ����������
          for j = 1:n
              temp=cell2mat(Table(i,j));
              Delta_Tau(temp(1,1),temp(1,2)) = Delta_Tau(temp(1,1),temp(1,2)) + Q/min_P;
          end
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = cell(m,n);
  
end
%% VI. �����ʾ
[Shortest_P,Shortest_index] = min(P_best);
max_P=-(Shortest_P-(n*(max(max(Val))+1)));
Shortest_Route = Route_best(Shortest_index,:);

for i=1:n
    temp=cell2mat(Shortest_Route(1,i));
    Plan(temp(1))=temp(2);
end

Pseries = -(P_best-n*(max(max(Val))+1));

end
