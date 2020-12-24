algo_cell = {'GRM';'SAP';'HCGSAP';'HCGDMEA';'ACO';'PSO';'GA';'SAGA'};
[l,~]=size(algo_cell);

num_missiles=50;
num_targets=40;
filename = [num2str(num_missiles),'vs',num2str(num_targets),'/MonteCarlo/'];
% for i=1:l
%     load([filename,algo_cell{i},'_max'],[algo_cell{i},'_max']);
% end

figure(1)
boxplot([GRM_max;SAP_max;HCGSAP_max;HCGDMEA_max;PSO_max;GA_max;SAGA_max;ACO_max]');
tick = {'GRM','SAP','HCGSAP','HCGDMEA','PSO','GA','SAGA','ACOSAGA'};
set(gca,'XTickLabel',tick);


% boxplot([HCGSAP_max;HCGDMEA_max]');

xlabel('算法名称');
ylabel('目标函数值/参照值');