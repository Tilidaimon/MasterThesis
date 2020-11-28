function [Plan,max_P,Pseries] = ACOfunction(Val)
%% 初始化参数
temp_Val=max(max(Val))+1-Val;
[p,q]=size(temp_Val);                     % 优势矩阵参数
n=max(p,q);
m = 100;                              % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
rho = 0.1;                           % 信息素挥发因子
Q = 1;                               % 常系数
Eta = 1./temp_Val;                          % 启发函数
Tau = ones(n,n)./(n*sum(sum(temp_Val))');   % 初始化信息素矩阵
Table = cell(m,n);                   % 路径记录表
iter = 1;                            % 迭代次数初值
iter_max = 50;                      % 最大迭代次数 
Route_best = cell(iter_max,n);      % 各代最佳路径       
P_best = zeros(iter_max,1);     % 各代最佳路径的长度  
P_ave = zeros(iter_max,1);      % 各代路径的平均长度  


pthreshold=1;

while iter <= iter_max
     pthreshold=pthreshold/log(1+iter);
     % 随机产生各个蚂蚁的起点城市
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
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个分配情况选择
         for j = 2:n
            %% 将已访问的节点及其的行列都加入禁忌表
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
             allow = allocation_number(allow_tabu_index);  % 待访问的任务分配情况集合
             P = allow;
             last_index=(tabu(1,j-1)-1)*n+tabu(2,j-1);
             % 计算城市间转移概率
             for k = 1:length(allow)
                row=fix((allow(k)-1)/n)+1;
                column=allow(k)-(row-1)*n;
                 P(k) = Tau(row,column)^alpha * Eta(row,column)^beta;
             end
%              P = P/sum(P);
%             %%轮盘赌法选择下一个访问城市
%             Pc = cumsum(P);     
%             target_tabu_index = find(Pc >= rand); 
%             target = allow(target_tabu_index(1));
%             target_row=fix((target-1)/n)+1;
%             target_column=target-(target_row-1)*n;
%             Table{i,j} = [target_row,target_column];
            
            if rand<pthreshold
                P = P/sum(P);
                %轮盘赌法选择下一个访问城市
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
      % 计算各个蚂蚁的任务分配全程的优势值
      totalP = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          matrixRoute=cell2mat(Route);
          for j = 1:n
              totalP(i) = totalP(i)+temp_Val(matrixRoute(j*2-1),matrixRoute(j*2));
          end
      end
      % 计算最短路径距离及平均距离
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
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个任务计算
          for j = 1:n
              temp=cell2mat(Table(i,j));
              Delta_Tau(temp(1,1),temp(1,2)) = Delta_Tau(temp(1,1),temp(1,2)) + Q/min_P;
          end
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = cell(m,n);
  
end
%% VI. 结果显示
[Shortest_P,Shortest_index] = min(P_best);
max_P=-(Shortest_P-(n*(max(max(Val))+1)));
Shortest_Route = Route_best(Shortest_index,:);

for i=1:n
    temp=cell2mat(Shortest_Route(1,i));
    Plan(temp(1))=temp(2);
end

Pseries = -(P_best-n*(max(max(Val))+1));

end
