% clear; clc;
% load('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\GDSC\Raw_GDSC_X.mat')
% Gene_Expression_Cell_lines=Raw_GDSC_X(2:end,1);
% Gene_Expression_Data = cell2mat(Raw_GDSC_X(2:end,2:end));
[num,txt,raw]=xlsread('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\GDSC\GDSC_Cell_AUC.xlsx');
Cell_Cancer_type=txt(7:end,3);
Cancer_cell_type={'blood','nervous_system','soft_tissue','bone','skin','urogenital_system',...
    'lung','pancreas','aero_digestive_tract','breast','kidney',...
    'digestive_system','thyroid'};
for Drug=1:140
    Drug1_AUC=num(8:end,Drug+4);
    Drug_AUC=Drug1_AUC(~any(isnan(Drug1_AUC),2),:);
    Cell_Cancer_type_Drug=Cell_Cancer_type(~any(isnan(Drug1_AUC),2),:);
    cancer_cell_number=[];
    for j=1:length(Cancer_cell_type)
        cancer_cell_number{j,1}=strmatch(Cancer_cell_type(j), Cell_Cancer_type_Drug);
        AUC{j}=1-Drug_AUC(cancer_cell_number{j,1},1);
        LLL(j,Drug)=length(AUC{j});
        MMM(j,Drug)=mean(AUC{j});
        VVV(j,Drug)=var(AUC{j});
    end
end
% MM'