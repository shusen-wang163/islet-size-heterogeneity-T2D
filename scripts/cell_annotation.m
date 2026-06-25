%  Step 2.2.1 select islet cells and cell area large equal than 10, define  cell type by 20% marker
close all;
addpath(genpath('./'));

clear;clc;
load OriginalData/Data_MaskCell_CellProfiler.mat

% step 1: cell belong to islet and area large equal than 10
data   = data(data{:,"IsInIslet"}>0,:);
[s, r] = sort(data{:,"IsInIslet"});
data   = data(r,:);
I      = data{:,"area"}>=10;
data   = data(I,:);
I      = data{:,"Islet_Size"}>=10;
data   = data(I,:);

% step 2: islet cell, GLU INS SST PP area overlap large than 20%
Expressionthreshold = 0.2;

% step 2.2.2  islet cell exclude to be one type
% exclude overlap cell, each cell has to be GLU or INS or SST or PP positive
I_abundent_cells = find(max(data{:,["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]}')'>Expressionthreshold);
I_INS = find(data{I_abundent_cells,"EL_Overlap_141Pr_INS_Mask"} == max(data{I_abundent_cells,["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]}')');
data{I_abundent_cells(I_INS),["EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]} = -data{I_abundent_cells(I_INS),["EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]};
I_GLU = find(data{I_abundent_cells,"EL_Overlap_151Eu_GLU_Mask"} == max(data{I_abundent_cells,["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]}')');
data{I_abundent_cells(I_GLU),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]} = -data{I_abundent_cells(I_GLU),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]};
I_SST = find(data{I_abundent_cells,"EL_Overlap_159Tb_SST_Mask"} == max(data{I_abundent_cells,["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]}')');
data{I_abundent_cells(I_SST),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_153Eu_PP_Mask"]} = -data{I_abundent_cells(I_SST),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_153Eu_PP_Mask"]};
I_PP  = find(data{I_abundent_cells,"EL_Overlap_153Eu_PP_Mask"} == max(data{I_abundent_cells,["EL_Overlap_141Pr_INS_Mask","EL_Overlap_151Eu_GLU_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_153Eu_PP_Mask"]}')');
data{I_abundent_cells(I_PP),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_151Eu_GLU_Mask"]} = -data{I_abundent_cells(I_PP),["EL_Overlap_141Pr_INS_Mask","EL_Overlap_159Tb_SST_Mask","EL_Overlap_151Eu_GLU_Mask"]};

% data{data{:,"EL_Overlap_161Dy_CK19_Mask"}  > Expressionthreshold,'cell_type'} = 5;
% data{data{:,"EL_Overlap_146Nd_CD34_Mask"}  > Expressionthreshold,'cell_type'} = 6;
data{data{:,"EL_Overlap_151Eu_GLU_Mask"} > Expressionthreshold,'cell_type'}   = 1;
data{data{:,"EL_Overlap_141Pr_INS_Mask"} > Expressionthreshold,'cell_type'}   = 2;
data{data{:,"EL_Overlap_159Tb_SST_Mask"} > Expressionthreshold,'cell_type'}   = 3;
data{data{:,"EL_Overlap_153Eu_PP_Mask"}  > Expressionthreshold,'cell_type'}   = 4;


Alpha_cell_ids = data{data{:,"cell_type"}==1,'cell_id'};
Beta_cell_ids  = data{data{:,"cell_type"}==2,'cell_id'};
Delta_cell_ids = data{data{:,"cell_type"}==3,'cell_id'};
PP_cell_ids    = data{data{:,"cell_type"}==4,'cell_id'};


writetable(data,'Data_IsletCell_187090_cells_CellProfiler_type_20_size_10.csv');
save Data_IsletCell_187090_cells_CellProfiler_type_20_size_10.mat data

data_pancreas = grpstats(data,"PancreasID");
data_pancreas = data_pancreas(:,{'PancreasID', 'GroupCount'});
writetable(data_pancreas,'cell_number_per_pancreas.csv');

data_type = grpstats(data,"type");
data_type = data_type(:,{'type', 'GroupCount'});
writetable(data_type,'cell_number_per_type.csv');

data_cell_type = grpstats(data,"cell_type");
data_cell_type = data_cell_type(:,{'cell_type', 'GroupCount'});
writetable(data_type,'cell_number_per_cell_type.csv');



%% Islet : islet cell ratio based on threshold 20%

close all;
clear;
clc;
load Data_IsletCell_187090_cells_CellProfiler_type_20_size_10.mat
load OriginalData/Islet_5313_10pixel.mat

Expressionthreshold = 0.2;


for i=1:size(data_islet,1)
    i
islet_id = data_islet{i,'islet_id'};
I = data{:,"IsInIslet"} == islet_id;
% data_islet{i,"INS_cell_number"}     = sum(abs(data{I,"EL_Overlap_141Pr_INS_Mask"}) > Expressionthreshold);
% data_islet{i,"GLU_cell_number"}     = sum(abs(data{I,"EL_Overlap_151Eu_GLU_Mask"}) > Expressionthreshold);
% data_islet{i,"SST_cell_number"}     = sum(abs(data{I,"EL_Overlap_159Tb_SST_Mask"}) > Expressionthreshold);
% data_islet{i,"PP_cell_number"}      = sum(abs(data{I,"EL_Overlap_153Eu_PP_Mask"}) > Expressionthreshold);

data_islet{i,"GLU_cell_number"}     = sum(abs(data{I,"cell_type"}) == 1);
data_islet{i,"INS_cell_number"}     = sum(abs(data{I,"cell_type"}) == 2);
data_islet{i,"SST_cell_number"}     = sum(abs(data{I,"cell_type"}) == 3);
data_islet{i,"PP_cell_number"}      = sum(abs(data{I,"cell_type"}) == 4);

data_islet{i,"CD34_cell_number"}    = sum(abs(data{I,"EL_Overlap_146Nd_CD34_Mask"}) > Expressionthreshold);
data_islet{i,"CK19_cell_number"}    = sum(abs(data{I,"EL_Overlap_161Dy_CK19_Mask"}) > Expressionthreshold);
data_islet{i,"KI67cellnumber"}      = sum(abs(data{I,"EL_Overlap_168Er_KI67_Mask"}) > Expressionthreshold);
data_islet{i,"SMA19cellnumber"}     = sum(abs(data{I,"EL_Overlap_167Er_SMA_Mask"}) > Expressionthreshold);

end

    data_islet{:,"cellnumber"}       = data_islet{:,"INS_cell_number"} + data_islet{:,"GLU_cell_number"} + data_islet{:,"SST_cell_number"}+data_islet{:,"PP_cell_number"};

data_islet{:,'INS_ratio'} = data_islet{:,'INS_cell_number'}./data_islet{:,'cellnumber'};
data_islet{:,'GLU_ratio'} = data_islet{:,'GLU_cell_number'}./data_islet{:,'cellnumber'};
data_islet{:,'SST_ratio'} = data_islet{:,'SST_cell_number'}./data_islet{:,'cellnumber'};
data_islet{:,'PP_ratio'}  = data_islet{:,'PP_cell_number'}./data_islet{:,'cellnumber'};


save Islet_5313_10pixel_cell_ratio_threshold_20.mat data_islet
writetable(data_islet,'Islet_5313_10pixel_cell_ratio_threshold_20.csv');

data_pancreas = grpstats(data_islet,"PancreasID");
data_pancreas = data_pancreas(:,{'PancreasID', 'GroupCount'});
writetable(data_pancreas,'Islet_number_per_pancreas_threshold_20.csv');

data_type = grpstats(data_islet,"type");
data_type = data_type(:,{'type', 'GroupCount'});
writetable(data_type,'Islet_number_per_type_threshold_20.csv');

