clear
clc

%% 导弹与目标位置初始化
% Nm = 100;
% Nt = 100;
Mcur = 10;
Ncur = 10;
JSFP_mean_time = zeros(Mcur,Ncur);
GRM_mean_time = zeros(Mcur,Ncur);
SAP_mean_time = zeros(Mcur,Ncur);
PSO_mean_time = zeros(Mcur,Ncur);
ACO_mean_time = zeros(Mcur,Ncur);
SAGA_mean_time = zeros(Mcur,Ncur);
GA_mean_time = zeros(Mcur,Ncur);
HCGSAP_mean_time = zeros(Mcur,Ncur);
HCGDMEA_mean_time = zeros(Mcur,Ncur);

JSFP_25tile_time = zeros(Mcur,Ncur);
GRM_25tile_time = zeros(Mcur,Ncur);
SAP_25tile_time = zeros(Mcur,Ncur);
PSO_25tile_time = zeros(Mcur,Ncur);
GA_25tile_time = zeros(Mcur,Ncur);
SAGA_25tile_time = zeros(Mcur,Ncur);
ACO_25tile_time = zeros(Mcur,Ncur);
HCGSAP_25tile_time = zeros(Mcur,Ncur);
HCGDMEA_25tile_time = zeros(Mcur,Ncur);

JSFP_75tile_time = zeros(Mcur,Ncur);
GRM_75tile_time = zeros(Mcur,Ncur);
SAP_75tile_time = zeros(Mcur,Ncur);
PSO_75tile_time = zeros(Mcur,Ncur);
GA_75tile_time = zeros(Mcur,Ncur);
SAGA_75tile_time = zeros(Mcur,Ncur);
ACO_75tile_time = zeros(Mcur,Ncur);
HCGSAP_75tile_time = zeros(Mcur,Ncur);
HCGDMEA_75tile_time = zeros(Mcur,Ncur);

JSFP_max_time = zeros(Mcur,Ncur);
GRM_max_time = zeros(Mcur,Ncur);
SAP_max_time = zeros(Mcur,Ncur);
PSO_max_time = zeros(Mcur,Ncur);
ACO_max_time = zeros(Mcur,Ncur);
SAGA_max_time = zeros(Mcur,Ncur);
GA_max_time = zeros(Mcur,Ncur);
HCGSAP_max_time = zeros(Mcur,Ncur);
HCGDMEA_max_time = zeros(Mcur,Ncur);


JSFP_min_time = zeros(Mcur,Ncur);
GRM_min_time = zeros(Mcur,Ncur);
SAP_min_time = zeros(Mcur,Ncur);
PSO_min_time = zeros(Mcur,Ncur);
ACO_min_time = zeros(Mcur,Ncur);
SAGA_min_time = zeros(Mcur,Ncur);
GA_min_time = zeros(Mcur,Ncur);
HCGSAP_min_time = zeros(Mcur,Ncur);
HCGDMEA_min_time = zeros(Mcur,Ncur);


JSFP_median_time = zeros(Mcur,Ncur);
GRM_median_time = zeros(Mcur,Ncur);
SAP_median_time = zeros(Mcur,Ncur);
PSO_median_time = zeros(Mcur,Ncur);
GA_median_time = zeros(Mcur,Ncur);
SAGA_median_time = zeros(Mcur,Ncur);
ACO_median_time = zeros(Mcur,Ncur);
HCGSAP_median_time = zeros(Mcur,Ncur);
HCGDMEA_median_time = zeros(Mcur,Ncur);

JSFP_25tile = zeros(Mcur,Ncur);
GRM_25tile = zeros(Mcur,Ncur);
SAP_25tile = zeros(Mcur,Ncur);
PSO_25tile = zeros(Mcur,Ncur);
GA_25tile = zeros(Mcur,Ncur);
SAGA_25tile = zeros(Mcur,Ncur);
ACO_25tile = zeros(Mcur,Ncur);
HCGSAP_25tile = zeros(Mcur,Ncur);
HCGDMEA_25tile = zeros(Mcur,Ncur);

JSFP_75tile = zeros(Mcur,Ncur);
GRM_75tile = zeros(Mcur,Ncur);
SAP_75tile = zeros(Mcur,Ncur);
PSO_75tile = zeros(Mcur,Ncur);
GA_75tile = zeros(Mcur,Ncur);
SAGA_75tile = zeros(Mcur,Ncur);
ACO_75tile = zeros(Mcur,Ncur);
HCGSAP_75tile = zeros(Mcur,Ncur);
HCGDMEA_75tile = zeros(Mcur,Ncur);

JSFP_mean = zeros(Mcur,Ncur);
GRM_mean = zeros(Mcur,Ncur);
SAP_mean = zeros(Mcur,Ncur);
PSO_mean = zeros(Mcur,Ncur);
GA_mean = zeros(Mcur,Ncur);
SAGA_mean = zeros(Mcur,Ncur);
ACO_mean = zeros(Mcur,Ncur);
HCGSAP_mean = zeros(Mcur,Ncur);
HCGDMEA_mean = zeros(Mcur,Ncur);

JSFP_median = zeros(Mcur,Ncur);
GRM_median = zeros(Mcur,Ncur);
SAP_median = zeros(Mcur,Ncur);
PSO_median = zeros(Mcur,Ncur);
GA_median = zeros(Mcur,Ncur);
SAGA_median = zeros(Mcur,Ncur);
ACO_median = zeros(Mcur,Ncur);
HCGSAP_median = zeros(Mcur,Ncur);
HCGDMEA_median = zeros(Mcur,Ncur);

JSFP_maxu = zeros(Mcur,Ncur);
GRM_maxu = zeros(Mcur,Ncur);
SAP_maxu = zeros(Mcur,Ncur);
PSO_maxu = zeros(Mcur,Ncur);
GA_maxu = zeros(Mcur,Ncur);
SAGA_maxu = zeros(Mcur,Ncur);
ACO_maxu = zeros(Mcur,Ncur);
HCGSAP_maxu = zeros(Mcur,Ncur);
HCGDMEA_maxu = zeros(Mcur,Ncur);

JSFP_minu = zeros(Mcur,Ncur);
GRM_minu = zeros(Mcur,Ncur);
SAP_minu = zeros(Mcur,Ncur);
PSO_minu = zeros(Mcur,Ncur);
GA_minu = zeros(Mcur,Ncur);
SAGA_minu = zeros(Mcur,Ncur);
ACO_minu = zeros(Mcur,Ncur);
HCGSAP_minu = zeros(Mcur,Ncur);
HCGDMEA_minu = zeros(Mcur,Ncur);

for Nt = 10:10:100
    for Nm=Nt:10:100
        model = DynamicMissileAndTarget(Nm,Nt);
        model = RandomSetSituation(model);
        
        
        %% 分配
        Epis = 200;
        
        Times = 100;
        JSFP_time = zeros(1,Times);
        GRM_time = zeros(1,Times);
        SAP_time = zeros(1,Times);
        HCGSAP_time = zeros(1,Times);
        HCGDMEA_time = zeros(1,Times);
        PSO_time = zeros(1,Times);
        ACO_time = zeros(1,Times);
        SAGA_time = zeros(1,Times);
        GA_time = zeros(1,Times);
        
        JSFP_max = zeros(1,Times);
        GRM_max = zeros(1,Times);
        SAP_max = zeros(1,Times);
        HCGSAP_max = zeros(1,Times);
        HCGDMEA_max = zeros(1,Times);
        PSO_max = zeros(1,Times);
        ACO_max = zeros(1,Times);
        SAGA_max = zeros(1,Times);
        GA_max = zeros(1,Times);
        
        
%         IGA_time = zeros(1,Times);
        for r = 1:Times
            %JSFP算法
            [JSFP_assign, JSFP_global_save, JSFP_time(r)] = JSFPfunction(model,Epis);
            JSFP_max(r) = max(JSFP_global_save);
            
            % GRMFMI算法
            [GRM_assign, GRM_global_save, GRM_time(r)] = GRMfunction(model,Epis);
            GRM_max(r) = max(GRM_global_save);
            
            % SAP算法
            tic
            [SAP_lagrange2_assign, SAP_global_save,SAP_time(r)] = SAP_lagrange_function2(model,Epis);
            SAP_max(r) = max(SAP_global_save);
            
            % PSO算法
            tic
            [PSO_extend_assign, PSO_global_save, PSO_time(r)] = PSOfunction(model);
            %PSO_assign = DecodeFromExtendPlan(model,PSO_extend_assign);
            PSO_max(r) = max(PSO_global_save);
            
            % HCGSAP算法
            [HCGSAP_assign, HCGSAP_global_save, HCGSAP_time(r)] = HCGFunctionSAP(model,Epis);
            HCGSAP_max(r) = max(HCGSAP_global_save);
            
            % HCG算法
            [HCGDMEA_assign, HCGDMEA_global_save, HCGDMEA_time(r)] = HCGFunction(model,Epis);
            HCGDMEA_max(r) = max(HCGDMEA_global_save);

            
            
            numLimit = max(model.num_missiles, sum(model.target_require_num_list));
            s = zeros(20,numLimit); % 预分配内存
            for i=1:20
                s(i,:) = randperm(numLimit);
            end
            
            % SAGA算法
            [SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time(r)] = SAGAfunction(s,model,20,200,0.4,0.4, "PMX", "EM",100,0.95);
            % SAGA_assign = DecodeFromExtendPlan(model,SAGA_extend_plan);
            SAGA_max(r) = max(SAGA_global_save);
            
            % ACOVFSAGA
            popsize=20;steps=200;Pcross=0.4;Pmutate=0.4;
            Top = 50; Dop = 0.95;
            [ACOVFSAGA_plan, ACO_global_save,ACO_mean_save, ACO_time(r)] = ACOVFSAGA(model,2*popsize,steps,Pcross,Pmutate,Top,Dop);
            % ACO_plan = DecodeFromExtendPlan(model,ACOVFSAGA_plan);
            ACO_max(r) = max(ACO_global_save);
            
            % GA
            [GA_extend_plan,GA_global_save,GA_mean_save, GA_time(r)] = GA1(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
            % GA_assign = DecodeFromExtendPlan(model,GA_extend_plan);
            GA_max(r) = max(GA_global_save);
            
%             [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
%             % IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);
%             IGA_max(r) = max(IGA_global_save);
            
            %归一化
            
%             IGA_max(r) = IGA_max(r)/maxU;
            
            
        end
        
        maxU = max([JSFP_max,GRM_max,SAP_max,HCGSAP_max,HCGDMEA_max,PSO_max,ACO_max,SAGA_max,GA_max]);
        
        JSFP_max = JSFP_max/maxU;
        GRM_max = GRM_max/maxU;
        SAP_max = SAP_max/maxU;
        
        PSO_max = PSO_max/maxU;
        ACO_max = ACO_max/maxU;
        SAGA_max = SAGA_max/maxU;
        GA_max = GA_max/maxU;
        
        indexM = Nm/10;
        indexN = Nt/10;
        
        % 最优值第一四分位数
        JSFP_25tile(indexM,indexN) = prctile(JSFP_max,25);
        GRM_25tile(indexM,indexN) = prctile(GRM_max,25);
        SAP_25tile(indexM,indexN) = prctile(SAP_max,25);
        PSO_25tile(indexM,indexN) = prctile(PSO_max,25);
        ACO_25tile(indexM,indexN) = prctile(ACO_max,25);
        SAGA_25tile(indexM,indexN) = prctile(SAGA_max,25);
        GA_25tile(indexM,indexN) = prctile(GA_max,25);
        HCGSAP_25tile(indexM,indexN) = prctile(HCGSAP_max,25);
        HCGDMEA_25tile(indexM,indexN) = prctile(HCGDMEA_max,25);
        
        % 最优值第三四分位数
        JSFP_75tile(indexM,indexN) = prctile(JSFP_max,75);
        GRM_75tile(indexM,indexN) = prctile(GRM_max,75);
        SAP_75tile(indexM,indexN) = prctile(SAP_max,75);
        PSO_75tile(indexM,indexN) = prctile(PSO_max,75);
        ACO_75tile(indexM,indexN) = prctile(ACO_max,75);
        SAGA_75tile(indexM,indexN) = prctile(SAGA_max,75);
        GA_75tile(indexM,indexN) = prctile(GA_max,75);
        HCGSAP_75tile(indexM,indexN) = prctile(HCGSAP_max,75);
        HCGDMEA_75tile(indexM,indexN) = prctile(HCGDMEA_max,75);
        
        % 最优值中位数
        JSFP_median(indexM,indexN) = median(JSFP_max);
        GRM_median(indexM,indexN) = median(GRM_max);
        SAP_median(indexM,indexN) = median(SAP_max);
        PSO_median(indexM,indexN) = median(PSO_max);
        ACO_median(indexM,indexN) = median(ACO_max);
        SAGA_median(indexM,indexN) = median(SAGA_max);
        GA_median(indexM,indexN) = median(GA_max);
        HCGSAP_median(indexM,indexN) = median(HCGSAP_max);
        HCGDMEA_median(indexM,indexN) = median(HCGDMEA_max);
        
        % 最优值平均值
        JSFP_mean(indexM,indexN) = mean(JSFP_max);
        GRM_mean(indexM,indexN) = mean(GRM_max);
        SAP_mean(indexM,indexN) = mean(SAP_max);
        PSO_mean(indexM,indexN) = mean(PSO_max);
        ACO_mean(indexM,indexN) = mean(ACO_max);
        SAGA_mean(indexM,indexN) = mean(SAGA_max);
        GA_mean(indexM,indexN) = mean(GA_max);
        HCGSAP_mean(indexM,indexN) = mean(HCGSAP_max);
        HCGDMEA_mean(indexM,indexN) = mean(HCGDMEA_max);
        
        % 最大值
        JSFP_maxu(indexM,indexN) = max(JSFP_max);
        GRM_maxu(indexM,indexN) = max(GRM_max);
        SAP_maxu(indexM,indexN) = max(SAP_max);
        PSO_maxu(indexM,indexN) = max(PSO_max);
        ACO_maxu(indexM,indexN) = max(ACO_max);
        SAGA_maxu(indexM,indexN) = max(SAGA_max);
        GA_maxu(indexM,indexN) = max(GA_max);
        HCGSAP_maxu(indexM,indexN) = max(HCGSAP_max);
        HCGDMEA_maxu(indexM,indexN) = max(HCGDMEA_max);
        
        % 最小值
        JSFP_minu(indexM,indexN) = min(JSFP_max);
        GRM_minu(indexM,indexN) = min(GRM_max);
        SAP_minu(indexM,indexN) = min(SAP_max);
        PSO_minu(indexM,indexN) = min(PSO_max);
        ACO_minu(indexM,indexN) = min(ACO_max);
        SAGA_minu(indexM,indexN) = min(SAGA_max);
        GA_minu(indexM,indexN) = min(GA_max);
        HCGSAP_minu(indexM,indexN) = min(HCGSAP_max);
        HCGDMEA_minu(indexM,indexN) = min(HCGDMEA_max);
        
        % 时间第一四分位数
        JSFP_25tile_time(indexM,indexN) = prctile(JSFP_time,25);
        GRM_25tile_time(indexM,indexN) = prctile(GRM_time,25);
        SAP_25tile_time(indexM,indexN) = prctile(SAP_time,25);
        PSO_25tile_time(indexM,indexN) = prctile(PSO_time,25);
        ACO_25tile_time(indexM,indexN) = prctile(ACO_time,25);
        SAGA_25tile_time(indexM,indexN) = prctile(SAGA_time,25);
        GA_25tile_time(indexM,indexN) = prctile(GA_time,25);
        HCGSAP_25tile_time(indexM,indexN) = prctile(HCGSAP_time,25);
        HCGDMEA_25tile_time(indexM,indexN) = prctile(HCGDMEA_time,25);
        
        % 时间第三四分位数
        JSFP_75tile_time(indexM,indexN) = prctile(JSFP_time,75);
        GRM_75tile_time(indexM,indexN) = prctile(GRM_time,75);
        SAP_75tile_time(indexM,indexN) = prctile(SAP_time,75);
        PSO_75tile_time(indexM,indexN) = prctile(PSO_time,75);
        ACO_75tile_time(indexM,indexN) = prctile(ACO_time,75);
        SAGA_75tile_time(indexM,indexN) = prctile(SAGA_time,75);
        GA_75tile_time(indexM,indexN) = prctile(GA_time,75);
        HCGSAP_75tile_time(indexM,indexN) = prctile(HCGSAP_time,75);
        HCGDMEA_75tile_time(indexM,indexN) = prctile(HCGDMEA_time,75);
        
        % 时间中位数
        JSFP_median_time(indexM,indexN) = median(JSFP_time);
        GRM_median_time(indexM,indexN) = median(GRM_time);
        SAP_median_time(indexM,indexN) = median(SAP_time);
        PSO_median_time(indexM,indexN) = median(PSO_time);
        ACO_median_time(indexM,indexN) = median(ACO_time);
        SAGA_median_time(indexM,indexN) = median(SAGA_time);
        GA_median_time(indexM,indexN) = median(GA_time);
        HCGSAP_median_time(indexM,indexN) = median(HCGSAP_time);
        HCGDMEA_median_time(indexM,indexN) = median(HCGDMEA_time);
        
        % 时间平均值
        JSFP_mean_time(indexM,indexN) = mean(JSFP_time);
        GRM_mean_time(indexM,indexN) = mean(GRM_time);
        SAP_mean_time(indexM,indexN) = mean(SAP_time);
        PSO_mean_time(indexM,indexN) = mean(PSO_time);
        ACO_mean_time(indexM,indexN) = mean(ACO_time);
        SAGA_mean_time(indexM,indexN) = mean(SAGA_time);
        GA_mean_time(indexM,indexN) = mean(GA_time);
        HCGSAP_mean_time(indexM,indexN) = mean(HCGSAP_time);
        HCGDMEA_mean_time(indexM,indexN) = mean(HCGDMEA_time);
%         IGA_mean_time(indexM,indexN) = mean(IGA_time);
        
        % 时间最大值
        JSFP_max_time(indexM,indexN) = max(JSFP_time);
        GRM_max_time(indexM,indexN) = max(GRM_time);
        SAP_max_time(indexM,indexN) = max(SAP_time);
        PSO_max_time(indexM,indexN) = max(PSO_time);
        ACO_max_time(indexM,indexN) = max(ACO_time);
        SAGA_max_time(indexM,indexN) = max(SAGA_time);
        GA_max_time(indexM,indexN) = max(GA_time);
        HCGSAP_max_time(indexM,indexN) = max(HCGSAP_time);
        HCGDMEA_max_time(indexM,indexN) = max(HCGDMEA_time);
        
        % 时间最小值
        JSFP_min_time(indexM,indexN) = min(JSFP_time);
        GRM_min_time(indexM,indexN) = min(GRM_time);
        SAP_min_time(indexM,indexN) = min(SAP_time);
        PSO_min_time(indexM,indexN) = min(PSO_time);
        ACO_min_time(indexM,indexN) = min(ACO_time);
        SAGA_min_time(indexM,indexN) = min(SAGA_time);
        GA_min_time(indexM,indexN) = min(GA_time);
        HCGSAP_min_time(indexM,indexN) = min(HCGSAP_time);
        HCGDMEA_min_time(indexM,indexN) = min(HCGDMEA_time);
    end
end

algo_cell = {'JSFP';'GRM';'SAP';'HCGSAP';'HCGDMEA';'ACO';'PSO';'GA';'IGA';'SAGA'};
[l,~]=size(algo_cell);
filename = 'TimeTest/';
for i=1:l
    save([filename,algo_cell{i},'_mean_time'],[algo_cell{i},'_mean_time']);
    save([filename,algo_cell{i},'_25tile_time'],[algo_cell{i},'_25tile_time']);
    save([filename,algo_cell{i},'_75tile_time'],[algo_cell{i},'_75tile_time']);
    save([filename,algo_cell{i},'_median_time'],[algo_cell{i},'_median_time']);
    save([filename,algo_cell{i},'_max_time'],[algo_cell{i},'_max_time']);
    save([filename,algo_cell{i},'_min_time'],[algo_cell{i},'_min_time']);
    
    save([filename,algo_cell{i},'_mean'],[algo_cell{i},'_mean']);
    save([filename,algo_cell{i},'_25tile'],[algo_cell{i},'_25tile']);
    save([filename,algo_cell{i},'_75tile'],[algo_cell{i},'_75tile']);
    save([filename,algo_cell{i},'_median'],[algo_cell{i},'_median']);
    save([filename,algo_cell{i},'_maxu'],[algo_cell{i},'_maxu']);
    save([filename,algo_cell{i},'_minu'],[algo_cell{i},'_minu']);
    %     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
end