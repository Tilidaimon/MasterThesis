[~,lsap] = size(SAP_global_save);
[~,lpso] = size(PSO_global_save);
[~,ljsfp] = size(JSFP_global_save);
[~,lgrm] = size(GRM_global_save);
MaxE = min(max([lsap,lpso,ljsfp,lgrm]),200);
e = 1:MaxE;

% algo_cell = {'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
% [l,~]=size(algo_cell);
% filename = [num2str(model.num_missiles),'vs',num2str(model.num_targets),'/'];
% for i=1:l
%     load([filename,algo_cell{i},'_global'],[algo_cell{i},'_global_save']);
%     load([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end

figure(1);
plot(e,ACO_global_save(1:MaxE))
hold on
plot(e,SAGA_global_save(1:MaxE));
% plot(e,IGABestSave(1:MaxE));
plot(e,GA_global_save(1:MaxE));

if MaxE>lsap
    SAP_expand = [SAP_global_save,SAP_global_save(lsap)*ones(1,MaxE-lsap)];
else
    SAP_expand = SAP_global_save(1:MaxE);
end
plot(e,SAP_expand);

if MaxE>ljsfp
    JSFP_expand = [JSFP_global_save,JSFP_global_save(ljsfp)*ones(1,MaxE-lsap)];
else
    JSFP_expand = JSFP_global_save(1:MaxE);
end
plot(e,JSFP_expand(1:MaxE));

if MaxE>lgrm
    GRM_expand = [GRM_global_save,GRM_global_save(lgrm)*ones(1,MaxE-lsap)];
else
    GRM_expand = GRM_global_save(1:MaxE);
end
plot(e,GRM_expand(1:MaxE));




PSO_expand = [PSO_global_save,PSO_global_save(lpso)*ones(1,MaxE-lpso)];
plot(e,PSO_expand);
axis([1,MaxE,min([min(SAP_global_save),min(GRM_global_save),min(JSFP_global_save)]),1.05]);
legend('ACOVFSAGA','SAGA','GA','SAP','JSFP','GRMFMI','PSO')