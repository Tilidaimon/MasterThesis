algo_cell = {'HCGSAP';'HCGDMEA';'JSFP';'GRMFMI';'SAP'};
[l,~]=size(algo_cell);
filename = 'TimeTest3/';


C =linspecer(l,'qualitative'); 

% 相同目标数量，不同导弹数量
MaxNm = 60;
MinNm = 20;
Nt = 10;
Nml = (MaxNm-MinNm)/Nt+1;
for i=1:l
    eval([algo_cell{i},'_max=zeros(Nml,100);'])
    eval([algo_cell{i},'_time=zeros(Nml,100);'])
    k=1;
    for Nm = MinNm:10:MaxNm
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_max']);
        eval([algo_cell{i},'_max(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max;'])
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_time']);
        eval([algo_cell{i},'_time(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_time;'])
        k = k+1;
    end
end
GRMFMI_max(5,GRMFMI_max(5,:)==0)=median(GRMFMI_max(5,:));
for k=1:Nml
    GRMFMI_time(k,:) = GRMFMI_time(k,:)/(k*10);
    JSFP_time(k,:) = JSFP_time(k,:)/(k*10);
    HCGDMEA_time(k,:) = HCGDMEA_time(k,:)/(k*10);
end

%% 结果
figure(3)


axes('NextPlot','replacechildren','ColorOrder',C); 
Total_max = zeros(100,Nml*l);
Total_time = zeros(100,Nml*l);
for j=1:Nml
    for i=1:l
        Total_max(:,(j-1)*l+i) = eval([algo_cell{i},'_max(',num2str(j),',:)'';']);
        Total_time(:,(j-1)*l+i) = eval([algo_cell{i},'_time(',num2str(j),',:)'';']);
    end
end
% boxplot([JSFP_max;GRM_max;SAP_max;PSO_max;GA_max;SAGA_max;ACO_max]');
% tick = {'JSFP','GRM','SAP','PSO','GA','SAGA','ACOSAGA'};
boxplot(Total_max)

set(gca,'XTick',4:l:l*(Nml-1)+4,'XTickLabel',MinNm:10:MaxNm);
set(gca,'YGrid','on');
h = findobj(gca,'Tag','Box');
for k=1:Nml
    for j=1:l
        
        box_color((k-1)*l+j)=patch(get(h((k-1)*l+j),'XData'),get(h((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end

hold on
plot(1:l:l*(Nml-1)+1,mean(HCGSAP_max,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    plot(j:l:l*(Nml-1)+j,eval(['mean(',algo_cell{j},'_max,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end
algo_cell1 = {'HCGSAP';'HCGDMCA';'JSFP';'GRMFMI';'PGSAP'};
xlabel('导弹数量');
ylabel('全局效用收敛值（归一化后）');
legend(box_color(l:-1:1),algo_cell1)
%% 时间
figure(2)
boxplot(Total_time)

set(gca,'XTick',3:l:l*(Nml-1)+3,'XTickLabel',MinNm:10:MaxNm);
set(gca,'YGrid','on');
h = findobj(gca,'Tag','Box');
for k=1:Nml
    for j=1:l
        
        box_color((k-1)*l+j)=patch(get(h((k-1)*l+j),'XData'),get(h((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end
legend(box_color(l:-1:1),algo_cell2)
% hold on
% plot(1:5,median(JSFP_max,2)','.-','color',C(1,:),'MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(1:l:l*(Nml-1)+1,mean(HCGSAP_time,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    plot(j:l:l*(Nml-1)+j,eval(['mean(',algo_cell{j},'_time,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end
algo_cell2 = {'HCGSAP';'HCGDMCA';'JSFP';'GRMFMI';'PGSAP'};
xlabel('导弹数量');
ylabel('算法运行时间/s');
legend(box_color(l:-1:1),algo_cell2)