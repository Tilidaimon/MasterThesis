clear
clc
algo_cell = {'JSFP';'GRMFMI';'SAP';'ACOSAGA';'PSO';'GA';'SAGA'};
[l,~]=size(algo_cell);
filename = 'TimeTest3/';


C =linspecer(l,'qualitative'); 

% 相同目标数量，不同导弹数量
MinNt = 10;
Nm = 55;
Ntl = 5;
for i=1:l
    eval([algo_cell{i},'_max=zeros(Ntl,100);'])
    eval([algo_cell{i},'_time=zeros(Ntl,100);'])
    k=1;
    for Nt = MinNt:10:Nm
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_max']);
        eval([algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max(find(',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max==0))=mean('...
            ,algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max);']);
        eval([algo_cell{i},'_max(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max;'])
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_time']);
        eval([algo_cell{i},'_time(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_time;'])
        k = k+1;
    end
end
GRMFMI_max(GRMFMI_max==min(min(GRMFMI_max)))=mean(GRMFMI_max(4,:));
for k=1:5
    GRMFMI_time(k,:) = GRMFMI_time(k,:)/(k*10);
    JSFP_time(k,:) = JSFP_time(k,:)/(k*10);
end
axes('NextPlot','replacechildren','ColorOrder',C); 

%% 结果
figure(1)


axes('NextPlot','replacechildren','ColorOrder',C); 
Total_max = zeros(100,Ntl*l);
Total_time = zeros(100,Ntl*l);
for j=1:Ntl
    for i=1:l
        Total_max(:,(j-1)*l+i) = eval([algo_cell{i},'_max(',num2str(j),',:)'';']);
        Total_time(:,(j-1)*l+i) = eval([algo_cell{i},'_time(',num2str(j),',:)'';']);
    end
end
% boxplot([JSFP_max;GRM_max;SAP_max;PSO_max;GA_max;SAGA_max;ACO_max]');
% tick = {'JSFP','GRM','SAP','PSO','GA','SAGA','ACOSAGA'};
boxplot(Total_max)

set(gca,'XTick',4:l:l*(Ntl-1)+4,'XTickLabel',MinNt:10:Nm);
set(gca,'YGrid','on');
h1 = findobj(gca,'Tag','Box');
for k=1:Ntl
    for j=1:l
        box_color1((k-1)*l+j)=patch(get(h1((k-1)*l+j),'XData'),get(h1((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end
algo_cell1 = {'JSFP';'GRMFMI';'SAP';'ACOSAGA';'PSO';'GA';'SAGA'};
% legend(box_color(l:-1:1),algo_cell1)
hold on

p1(1)=plot(1:l:l*(Ntl-1)+1,mean(JSFP_max,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    p1(j)=plot(j:l:l*(Ntl-1)+j,eval(['mean(',algo_cell{j},'_max,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end

xlabel('目标数量');
ylabel('全局效用收敛值（归一化后）');

algo_cell2 = {'JSFP平均值';'GRMFMI平均值';'SAP平均值';'ACOSAGA平均值';'PSO平均值';'GA平均值';'SAGA平均值'};
legend(box_color1(l:-1:1),algo_cell1)
ah1=axes('position',get(gca,'position'),'visible','off');
legend(ah1,p1,algo_cell2)



%% 时间
figure(2)

boxplot(Total_time)

set(gca,'XTick',4:l:l*(Ntl-1)+4,'XTickLabel',MinNt:10:Nm);
set(gca,'YGrid','on');
h2 = findobj(gca,'Tag','Box');
for k=1:Ntl
    for j=1:l   
        box_color2((k-1)*l+j)=patch(get(h2((k-1)*l+j),'XData'),get(h2((k-1)*l+j),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end
algo_cell1 = {'JSFP';'GRMFMI';'SAP';'ACOSAGA';'PSO';'GA';'SAGA'};


% hold on
% plot(1:5,median(JSFP_max,2)','.-','color',C(1,:),'MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
p2(1)=plot(1:l:l*(Ntl-1)+1,mean(JSFP_time,2)','.-','color',C(1,:),'MarkerEdgeColor',C(1,:));
for j=2:l
    p2(j)=plot(j:l:l*(Ntl-1)+j,eval(['mean(',algo_cell{j},'_time,2)''']),'.-','color',C(j,:),'MarkerEdgeColor',C(j,:));
end
algo_cell2 = {'JSFP平均值';'GRMFMI平均值';'SAP平均值';'ACOSAGA平均值';'PSO平均值';'GA平均值';'SAGA平均值'};

xlabel('目标数量');
ylabel('算法运行时间/s');
ylim([0,7])
legend(box_color2(l:-1:1),algo_cell1)
ah2=axes('position',get(gca,'position'),'visible','off');
legend(ah2,p2,algo_cell2)

