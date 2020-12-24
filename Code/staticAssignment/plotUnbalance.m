
algo_cell = {'HCGSAP';'HCGDMEA';'GRMFMI';'SAP';'JSFP'};
[l,~]=size(algo_cell);
filename = 'TestUnbalance/';


C =linspecer(l,'qualitative'); 

% 相同目标数量，不同导弹数量
Nm = 20;
Nt = 10;
adiN = Nm-Nt;


for i=1:l
    eval([algo_cell{i},'_max=zeros(adiN,100);'])
    eval([algo_cell{i},'_time=zeros(adiN,100);'])
    k=1;
    for k = 1:adiN
        load([filename,algo_cell{i},'_',num2str(k),'_max']);
        eval([algo_cell{i},'_max(k,:)=',algo_cell{i},'_',num2str(k),'_max;'])
        load([filename,algo_cell{i},'_',num2str(k),'_time']);
        eval([algo_cell{i},'_time(k,:)=',algo_cell{i},'_',num2str(k),'_time;'])
    end
end

for k=1:adiN
    GRMFMI_time(k,:) = GRMFMI_time(k,:)/(k*10);
%     JSFP_time(k,:) = JSFP_time(k,:)/(k*10);
    HCGDMEA_time(k,:) = HCGDMEA_time(k,:)/(k*10);
%     JSFP_max(k,JSFP_max(k,:)==0) = median(JSFP_max(k,JSFP_max(k,:)~=0));
end

%% 结果
figure(3)


axes('NextPlot','replacechildren','ColorOrder',C); 
Total_max = zeros(100,adiN*l);
Total_time = zeros(100,adiN*l);
for j=1:adiN
    for i=1:l
        Total_max(:,(j-1)*l+i) = eval([algo_cell{i},'_max(',num2str(j),',:)'';']);
        Total_time(:,(j-1)*l+i) = eval([algo_cell{i},'_time(',num2str(j),',:)'';']);
    end
end
% boxplot([JSFP_max;GRM_max;SAP_max;PSO_max;GA_max;SAGA_max;ACO_max]');
% tick = {'JSFP','GRM','SAP','PSO','GA','SAGA','ACOSAGA'};
boxplot(Total_max)

set(gca,'XTick',3:l:l*(adiN-1)+3,'XTickLabel',2:1+adiN);
set(gca,'YGrid','on');
h = findobj(gca,'Tag','Box');
for k=1:adiN
    for j=1:l
        
        box_color((k-1)*l+j)=patch(get(h((k-1)*l+j),'XData'),get(h((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end

hold on
plot(1:l:l*(adiN-1)+1,mean(HCGSAP_max,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    plot(j:l:l*(adiN-1)+j,eval(['mean(',algo_cell{j},'_max,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end
algo_cell1 = {'HCGSAP';'HCGDMCA';'GRMFMI';'PGSAP';'JSFP'};
xlabel('b_{max}');
ylabel('全局效用收敛值（归一化后）');
legend(box_color(l:-1:1),algo_cell1)
%% 时间
figure(2)
boxplot(Total_time)

set(gca,'XTick',3:l:l*(adiN-1)+3,'XTickLabel',2:adiN+1);
set(gca,'YGrid','on');
h = findobj(gca,'Tag','Box');
for k=1:adiN
    for j=1:l
        
        box_color((k-1)*l+j)=patch(get(h((k-1)*l+j),'XData'),get(h((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end

% hold on
% plot(1:5,median(JSFP_max,2)','.-','color',C(1,:),'MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(1:l:l*(adiN-1)+1,mean(HCGSAP_time,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    plot(j:l:l*(adiN-1)+j,eval(['mean(',algo_cell{j},'_time,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end
algo_cell2 = {'HCGSAP';'HCGDMCA';'GRMFMI';'PGSAP';'JSFP'};
xlabel('b_{max}');
ylabel('算法运行时间/s');
legend(box_color(l:-1:1),algo_cell2)