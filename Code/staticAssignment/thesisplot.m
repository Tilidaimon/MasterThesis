num_missiles=100;
num_targets=100;
% algo_cell = {'JSFP';'GRM';'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
% [l,~]=size(algo_cell);
% filename = [num2str(num_missiles),'vs',num2str(num_targets),'/Convergence/'];

algo_cell = {'JSFP';'GRM';'PGSAP';'HCGSAP';'HCGDMEA';'ACO';'PSO';'GA';'SAGA'};
[l,~]=size(algo_cell);
filename = ['第四章/',num2str(model.num_missiles),'vs',num2str(model.num_targets)];

% for i=1:l
%     load([filename,algo_cell{i},'_global'],[algo_cell{i},'_global_save']);
%     load([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end

[lsap,~] = size(PGSAP_global_save);
[~,lpso] = size(PSO_global_save);
[ljsfp,~] = size(JSFP_global_save);
[lgrm,~] = size(GRM_global_save);
MaxE = min(max([lsap,lpso,ljsfp,lgrm]),500);
minE = 150;
e = minE:MaxE;

% algo_cell = {'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
% [l,~]=size(algo_cell);
% filename = [num2str(model.num_missiles),'vs',num2str(model.num_targets),'/'];
% for i=1:l
%     load([filename,algo_cell{i},'_global'],[algo_cell{i},'_global_save']);
%     load([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end



figure(1);



if MaxE>ljsfp
    JSFP_expand = [JSFP_global_save(1:end),JSFP_global_save(end)*ones(1,MaxE-lsap)];
else
    JSFP_expand = JSFP_global_save(1:MaxE);
end
plot(e,JSFP_expand(minE:MaxE));
hold on

if MaxE>lgrm
    GRM_expand = [GRM_global_save(1:end),GRM_global_save(end)*ones(1,MaxE-lsap)];
else
    GRM_expand = GRM_global_save(1:MaxE);
end
plot(e,GRM_expand(minE:MaxE));

if MaxE>lsap
    SAP_expand = [PGSAP_global_save(1:end),PGSAP_global_save(end)*ones(1,MaxE-lsap)];
else
    SAP_expand = PGSAP_global_save(1:MaxE);
end
plot(e,SAP_expand(minE:MaxE));

plot(e,HCGSAP_global_save(minE:MaxE));

plot(e,HCGDMEA_global_save(minE:MaxE));


PSO_expand = [PSO_global_save(1:end),PSO_global_save(end)*ones(1,MaxE-lpso)];
plot(e,PSO_expand(minE:MaxE));


plot(e,GA_global_save(minE:MaxE));

plot(e,SAGA_global_save(minE:MaxE));

plot(e,ACO_global_save(minE:MaxE))

% axis([1,MaxE,min([min(PGSAP_global_save),min(GRM_global_save),min(JSFP_global_save)]),1.05]);
% plot(e,IGABestSave(1:MaxE));

xlabel('迭代步数');
ylabel('目标函数值/参照值');


legend('JSFP','GRMFMI','PGSAP','HCGSAP','HCGDMEA','PSO','GA','SAGA','ACOSAGA');