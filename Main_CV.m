clear; clc;
load('Gene_Expression_Cell_lines_Full.mat')
load('Gene_Expression_Full.mat')
[num,txt,raw]=xlsread('CCLE_Drug_Sensitivity_Raziur.xlsx',9);
for i=1:size(num,1)
    data=char(raw(i+1,1));
    C = strsplit(data,'_');
    Cell_type{i,1}= strjoin(C(2:end),'_');
end
Cancer_cell={'CENTRAL_NERVOUS_SYSTEM','PROSTATE','URINARY_TRACT',   'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',...
    'KIDNEY','THYROID','SOFT_TISSUE','SKIN','SALIVARY_GLAND','OVARY',   'LUNG','BONE','ENDOMETRIUM',...
    'PANCREAS','BREAST','STOMACH','LARGE_INTESTINE','LIVER','UPPER_AERODIGESTIVE_TRACT',...
    'AUTONOMIC_GANGLIA','BILIARY_TRACT','PLEURA','OESOPHAGUS'};
for j=1:length(Cancer_cell)
    cancer_cell_number{j}=strmatch(Cancer_cell(j), Cell_type);
end
HLT_AUC=num(cancer_cell_number{1},12)/8;
HLT_cell=txt(cancer_cell_number{1}+1,1);
LUNG_AUC=num(cancer_cell_number{8},12)/8;
LUNG_cell=txt(cancer_cell_number{8}+1,1);
%%
Ind=[];
for ii=1:length(HLT_cell)
    if strmatch(HLT_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind=[Ind; strmatch(HLT_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind=[Ind; 0];
    end
    
end
L = Ind(Ind~=0);
HLT_gene=Gene_Expression_Data(:,L)';
LL=find(Ind==0);
HLT_AUC(LL)=[];
HLT_cell(LL)=[];
Perm_HLT = randperm(length(HLT_AUC));
HLT_AUC=HLT_AUC(Perm_HLT,:);
HLT_gene=HLT_gene(Perm_HLT,:);
HLT_AUC=[HLT_AUC 2*ones(length(HLT_AUC),1)];

Ind2=[];
for ii=1:length(LUNG_cell)
    if strmatch(LUNG_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind2=[Ind2; strmatch(LUNG_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind2=[Ind2; 0];
    end
end
L = Ind2(Ind2~=0);
LUNG_gene=Gene_Expression_Data(:,L)';
LL=find(Ind2==0);
LUNG_AUC(LL)=[];
LUNG_cell(LL)=[];
Perm_LUNG = randperm(length(LUNG_AUC));
LUNG_AUC=LUNG_AUC(Perm_LUNG,:);
LUNG_gene=LUNG_gene(Perm_LUNG,:);
LUNG_AUC=[LUNG_AUC ones(length(LUNG_AUC),1)];
%%
n_tree=100;
mtree=10;
min_leaf=2;
F=3;
finalX=[LUNG_gene; HLT_gene];
finalY=[LUNG_AUC; HLT_AUC];
[Xtrain,Xtest,Ytrain,Ytest,FoldedIndex]=CreateFoldedDataMRF(finalX,finalY,F);

for FF=1:F
    X1=Xtrain{FF};
    Y1=Ytrain{FF};
    Xt{FF}=Xtest{FF};%50-100
    Yt{FF}=Ytest{FF};
    
    model{FF} = build_forest(X1, Y1, n_tree, mtree, 1,3, min_leaf, 1);
    for C=1:2
        C_ind=find(C==Ytest{FF}(:,2));
        Xtt=Xt{FF}(C_ind,:);
        Y_hatVV{FF} = forest_predict(Xtt, model{FF}, C);
        Y_hatV{FF}(C_ind,:)=Y_hatVV{FF};
    end
end
Yactual=[];
YpredV=[];
for i=1:F
    Yactual=[Yactual;Yt{i}];
    YpredV=[YpredV;Y_hatV{i}];
end
[Fold,Ind23]=sort(cell2mat(FoldedIndex));
finalYY=Yactual(Ind23,:);
finalYpred=YpredV(Ind23,:);
[mean((finalY(:,1)-finalYpred(:,1)).^2) mean((finalY(:,1)-finalYpred(:,2)).^2)]