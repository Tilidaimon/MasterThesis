clear
clc

%% 导弹与目标位置初始化
Nm = 20;
Nt = 10;
adiN = Nm-Nt;
model = DynamicMissileAndTarget(Nm,Nt);
model = RandomSetSituation(model);
bmax = 8;
model.target_require_num_list = ones(1,model.num_targets);
    assign_flag = zeros(model.num_targets,1);
    leader = randi(model.num_targets);
    model.target_require_num_list(leader) = model.target_require_num_list(leader) + bmax;
    assign_flag(leader) = 1;
    residual = adiN - bmax;
    for i = 1:model.num_targets
        if assign_flag(i) == 0 && residual>0
            adiMissiles = randi(residual);
            while adiMissiles>bmax
                adiMissiles = randi(residual);
            end
            model.target_require_num_list(i) = model.target_require_num_list(i) + adiMissiles;
            residual = residual - adiMissiles;
        end
    end
    % 虚拟目标增广序列
    k=1;
    model.order_targets = zeros(model.num_missiles,1);
    for i=1:model.num_targets
        for j=1:model.target_require_num_list(i)
            model.order_targets(k) = i;
            k = k+1;
        end
    end
    
    % 目标价值
    model.Targets.value = 1000*ones(model.num_targets,1);
    model.Targets.value = model.target_require_num_list'.*model.Targets.value;
    model.extend_value = model.Targets.value(model.order_targets);

%% 穷举法求精确解(Exhaustive Attack method)
% tic;
% N = max(model.num_missiles, sum(model.target_require_num_list));
% allPerms = perms(1:N);
% numOfAll = size(allPerms,1);
% fEAM = zeros(numOfAll,1);
% for i=1:numOfAll
%     plan = allPerms(i,:);
%     fEAM(i) =  GlobalUtility(model, plan);
% end
% % plot(fEAM,'o');
% [maxF,maxIndex] = max(fEAM);
% bestEAMPlan = allPerms(maxIndex,:);
% bestEAM_assign = DecodeFromExtendPlan(model,bestEAMPlan);
% clear allPerms;
% toc

%% 分配
Epis = 500;

% JSFP算法
[JSFP_assign, JSFP_global_save, JSFP_time] = JSFPfunction(model,Epis);


% GRMFMI算法

[GRM_assign, GRM_global_save, GRM_time] = GRMfunction(model,Epis);

% SAP算法
[PGSAP_assign, PGSAP_global_save, PGSAP_time] = SAP_lagrange_function2(model,Epis);

e=3;
% HCGSAP算法
[HCGSAP_assign, HCGSAP_global_save, HCGSAP_time] = HCGFunctionSAP(e,model,Epis);

% HCG算法
[HCGDMEA_assign, HCGDMEA_global_save, HCGDMEA_time] = HCGFunction(model,Epis);

% [HCG2_assign, HCG2_global_save, HCG2_time] = HCGFunction2(model,Epis);

% PSO算法
[PSO_extend_assign, PSO_global_save, PSO_time] = PSOfunction(model);
PSO_assign = DecodeFromExtendPlan(model,PSO_extend_assign);

numLimit = max(model.num_missiles, sum(model.target_require_num_list));
s = zeros(20,numLimit); % 预分配内存
for i=1:20
    s(i,:) = randperm(numLimit);
end

% SAGA算法
[SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time] = SAGAfunction(s,model,20,500,0.4,0.4, "PMX", "EM",100,0.95);
SAGA_assign = DecodeFromExtendPlan(model,SAGA_extend_plan);

% ACOVFSAGA
popsize=20;steps=500;Pcross=0.4;Pmutate=0.4; 
Top = 50; Dop = 0.95;
[ACOVFSAGA_plan, ACO_global_save,ACO_mean_save, ACO_time] = ACOVFSAGA(model,2*popsize,steps,Pcross,Pmutate,Top,Dop);
ACO_plan = DecodeFromExtendPlan(model,ACOVFSAGA_plan);


% GA
[GA_extend_plan,GA_global_save,GA_mean_save, GA_time] = GA1(model,popsize,steps,Pcross,Pmutate, "PMX", "EM");
GA_assign = DecodeFromExtendPlan(model,GA_extend_plan);


% [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "EM");
% IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);



% 归一化
maxU = max([JSFP_global_save(end),GRM_global_save(end),PGSAP_global_save(end),HCGSAP_global_save(end),...
    HCGDMEA_global_save(end),PSO_global_save(end),SAGA_global_save(end),ACO_global_save(end),GA_global_save(end)]);

JSFP_global_save = JSFP_global_save/maxU;
SAGA_global_save = SAGA_global_save/maxU;
PGSAP_global_save = PGSAP_global_save/maxU;
GRM_global_save = GRM_global_save/maxU;
PSO_global_save = PSO_global_save/maxU;
ACO_global_save = ACO_global_save/maxU;
GA_global_save = GA_global_save/maxU;
HCGSAP_global_save = HCGSAP_global_save/maxU;
HCGDMEA_global_save = HCGDMEA_global_save/maxU;


algo_cell = {'JSFP';'GRM';'PGSAP';'HCGSAP';'HCGDMEA';'ACO';'PSO';'GA';'SAGA'};
[l,~]=size(algo_cell);
% filename = ['第四章/',num2str(model.num_missiles),'vs',num2str(model.num_targets)];
% for i=1:l
%     save([filename,algo_cell{i},'_global'],[algo_cell{i},'_global_save']);
%     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end



% tic
% [HCG_assign, HCG_global, HCG_global_series] = HCGFunction(model,Epis);
% HCG_time = toc;
% tic
% [HCG2_assign, HCG2_global, HCG2_global_series] = HCGFunction2(model,Epis);
% HCG_time = toc;
% tic
% [HCGSAP_assign, HCGSAP_global, HCGSAP_global_series] = HCGFunctionSAP(model,Epis);
% HCG2_time = toc;
% Ug_SAP_lagrange2 = GlobalUtility(SAP_lagrange2_assign,Target.value,T_go,J_opt);
% Ug_HCG = GlobalUtility(HCG_assign,Target.value,T_go,J_opt);
% Ug_HCGSAP = GlobalUtility(HCGSAP_assign,Target.value,T_go,J_opt);

% episode = 1:Epis;
% plot(episode,SAP_lagrange2_global_series,episode,HCG_global_series,episode,HCG2_global_series,episode,HCGSAP_global_series)%,episode,HCGZ_global_series);
% legend('SAP lagrange','HCG','HCG2','HCGSAP')%,'HCGZ')

% % Ug_SAP(r) = GlobalUtility(SAP_assign,Target.value,T_go,J_opt);
% Ug_SAP_lagrange(r) = GlobalUtility(SAP_lagrange_assign,Target.value,T_go,J_opt);
% Ug_SAP_lagrange2(r) = GlobalUtility(SAP_lagrange2_assign,Target.value,T_go,J_opt);
% Ug_HCG(r) = GlobalUtility(HCG_assign,Target.value,T_go,J_opt);
% Ug_HCG2(r) = GlobalUtility(HCG2_assign,Target.value,T_go,J_opt);
% 
% % [Time_SAP(r),J_SAP(r)] = Time_and_Energy(SAP_assign,T_go,J_opt);
% [Time_SAP_lagrange(r),J_SAP_lagrange(r)] = Time_and_Energy(SAP_lagrange_assign,T_go,J_opt);
% [Time_SAP_lagrange2(r),J_SAP_lagrange2(r)] = Time_and_Energy(SAP_lagrange2_assign,T_go,J_opt);
% [Time_HCG(r),J_HCG(r)] = Time_and_Energy(HCG_assign,T_go,J_opt);
% [Time_HCG2(r),J_HCG2(r)] = Time_and_Energy(HCG2_assign,T_go,J_opt);
% % 
% end
% 
% figure(1)
% boxplot([Ug_SAP_lagrange2;Ug_HCG;Ug_HCGSAP]');
% tick = {'SAPlagrange2','HCG','HCGSAP'};
% set(gca,'XTickLabel',tick);
% 
% figure(2)
% boxplot([Time_SAP_lagrange;Time_SAP_lagrange2;Time_HCG;Time_HCG2]');
% tick = {'SAPlagrange','SAPlagrange2','HCG','HCG2'};
% set(gca,'XTickLabel',tick);
% 
% figure(3)
% boxplot([J_SAP_lagrange;J_SAP_lagrange2;J_HCG;J_HCG2]');
% tick = {'SAPlagrange','SAPlagrange2','HCG','HCG2'};
% set(gca,'XTickLabel',tick);

