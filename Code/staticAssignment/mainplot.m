algo_cell = {'JSFP';'GRM';'SAP';'ACO';'PSO';'GA';'SAGA';'HCGSAP';'HCGDMEA'};
[l,~]=size(algo_cell);
filename = 'TimeTest2/';


C =linspecer(l-2,'qualitative'); 

% 相同目标数量，不同导弹数量
MaxNm = 50;
Nt = 10;
Nml = (MaxNm-Nt)/10+1;
for i=1:l
    eval([algo_cell{i},'_max=zeros(Nml,100);'])
    eval([algo_cell{i},'_time=zeros(Nml,100);'])
    k=1;
    for Nm = Nt:10:MaxNm
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_max']);
        eval([algo_cell{i},'_max(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max;'])
        load([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_time']);
        eval([algo_cell{i},'_time(k,:)=',algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_time;'])
        k = k+1;
    end
end


figure(1)
%boxplot(SAP_max');

p(1) = plot(mean(JSFP_max,2)','.-','color',C(1,:),'MarkerSize',15);
hold on

for j=2:3
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['plot(mean(',algo_cell{j},'_max,2)'',''o-'',''color'',C(',num2str(j),',:));'])
end
for j=4:5
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['p(',numstr(j),')=plot(mean(',algo_cell{j},'_max,2)'',''^-'',''color'',C(',num2str(j),',:));'])
end
for j=6:l
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['p(',numstr(j),')=plot(mean(',algo_cell{j},'_max,2)'',''*-'',''color'',C(',num2str(j),',:));'])
end
legend(p(1:l-2),algo_cell{1:l-2})
set(gca,'XTick',1:Nml,'XTickLabel',Nt:10:MaxNm)
xlabel('导弹数量')
ylabel('全局效用平均值（归一化后）')

figure(2)
p(1) = plot(mean(JSFP_time,2)','.-','color',C(1,:),'MarkerSize',15);
hold on

for j=2:3
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['p(',numstr(j),')=plot(mean(',algo_cell{j},'_time,2)'',''o-'',''color'',C(',num2str(j),',:));'])
end
for j=4:5
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['p(',numstr(j),')=plot(mean(',algo_cell{j},'_time,2)'',''^-'',''color'',C(',num2str(j),',:));'])
end
for j=6:l-2
    %eval(['boxplot(',algo_cell{j},'_max'');']);
    eval(['p(',numstr(j),')=plot(mean(',algo_cell{j},'_time,2)'',''*-'',''color'',C(',num2str(j),',:));'])
end
legend(p(1:l-2),algo_cell{1:l-2})
set(gca,'XTick',1:MaxNm/10,'XTickLabel',Nt:10:MaxNm)
xlabel('导弹数量')
ylabel('平均运行时间/s')

figure(3)


axes('NextPlot','replacechildren','ColorOrder',C); 
boxplot([JSFP_max;GRM_max;SAP_max;HCGSAP_max;HCGDMEA_max;PSO_max;GA_max;SAGA_max;ACO_max]');
tick = {'JSFP','GRM','SAP','HCGSAP','HCGDMEA','PSO','GA','SAGA','ACOSAGA'};
set(gca,'XTick',3:5:5*(l-1)+3,'XTickLabel',tick);
set(gca,'YGrid','on');
h = findobj(gca,'Tag','Box');
for j=1:l
    for k=1:Nml
        box_color((j-1)*Nml+k)=patch(get(h((j-1)*Nml+k),'XData'),get(h((j-1)*Nml+k),'YData'),C(l+1-j,:),'FaceAlpha',0.5);
    end
end
legend(box_color([5*(l-1)+3:-5:3]),algo_cell)
% hold on
% plot(1:5,median(JSFP_max,2)','.-','color',C(1,:),'MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');

xlabel('算法名称');
ylabel('全局效用平均值（归一化后）');