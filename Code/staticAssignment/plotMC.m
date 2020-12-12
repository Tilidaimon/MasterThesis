algo_cell = {'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
[l,~]=size(algo_cell);
filename = [num2str(model.num_missiles),'vs',num2str(model.num_targets),'/MonteCarlo/'];
for i=1:l
    load([filename,algo_cell{i},'_max'],[algo_cell{i},'_max']);
end

figure(1)
boxplot([SAP_max;PSO_max;ACO_max;SAGA_max;GA_max]');
tick = {'SAPPG','PSO','ACOVFSAGA','SAGA','GA'};
set(gca,'XTickLabel',tick);

