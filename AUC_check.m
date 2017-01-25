clear; clc;
load('C:\Users\razrahma\Dropbox\Probabilistic RF\MATLAB Files\Gene_Expression_Cell_lines_Full.mat')
load('C:\Users\razrahma\Dropbox\Probabilistic RF\MATLAB Files\Gene_Expression_Full.mat')
for Drug=1:24
    [num,txt,raw]=xlsread('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\CCLE\CCLE_Drug_Sensitivity_Raziur.xlsx',Drug);
    Cell_type=[];
    for i=1:size(num,1)
        data=char(raw(i+1,1));
        C = strsplit(data,'_');
        Cell_type{i,1}= strjoin(C(2:end),'_');
    end
    Cancer_cell={'CENTRAL_NERVOUS_SYSTEM','PROSTATE','URINARY_TRACT',   'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',...
        'KIDNEY','THYROID','SOFT_TISSUE','SKIN','SALIVARY_GLAND','OVARY',   'LUNG','BONE','ENDOMETRIUM',...
        'PANCREAS','BREAST','STOMACH','LARGE_INTESTINE','LIVER','UPPER_AERODIGESTIVE_TRACT',...
        'AUTONOMIC_GANGLIA','BILIARY_TRACT','PLEURA','OESOPHAGUS'};
    cancer_cell_number=[];
    for j=1:length(Cancer_cell)
        cancer_cell_number{j,1}=strmatch(Cancer_cell(j), Cell_type);
        AUC{j}=num(cancer_cell_number{j,1},12)/8;
        LLL(j,Drug)=length(AUC{j});
        MMM(j,Drug)=mean(AUC{j});
        VVV(j,Drug)=var(AUC{j});
    end
end
% MM'