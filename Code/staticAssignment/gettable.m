algo_cell = {'HCGSAP';'HCGDMEA';'JSFP';'GRMFMI';'SAP'};
[l,~]=size(algo_cell);

table = zeros(30,l);
for i=1:l
    data = eval([algo_cell{i},'_max']);
    table(1:3:30,i) = mean(data,2);
    table(2:3:30,i) = max(data')';
    table(3:3:30,i) = min(data');
end

timetable = zeros(30,l);
for i=1:l
    data = eval([algo_cell{i},'_time']);
    timetable(1:3:30,i) = mean(data,2);
    timetable(2:3:30,i) = max(data')';
    timetable(3:3:30,i) = min(data');
end