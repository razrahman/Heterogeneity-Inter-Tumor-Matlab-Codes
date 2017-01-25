clear; clc;
load('C:\Users\razrahma\Dropbox\Probabilistic RF\MATLAB Files\Gene_Expression_Cell_lines_Full.mat')
load('C:\Users\razrahma\Dropbox\Probabilistic RF\MATLAB Files\Gene_Expression_Full.mat')
Drug=3;
[num,txt,raw]=xlsread('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\CCLE\CCLE_Drug_Sensitivity_Raziur.xlsx',Drug);
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
Cancer2_AUC=num(cancer_cell_number{8},12)/8;
Cancer2_cell=txt(cancer_cell_number{8}+1,1);
Cancer1_AUC=num(cancer_cell_number{10},12)/8;
Cancer1_cell=txt(cancer_cell_number{10}+1,1);
%%
Ind=[];
for ii=1:length(Cancer2_cell)
     if strmatch(Cancer2_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind=[Ind; strmatch(Cancer2_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind=[Ind; 0];
    end
    
end
L = Ind(Ind~=0);
Cancer2_gene=Gene_Expression_Data(:,L)';
LL=find(Ind==0);
Cancer2_AUC(LL)=[];
Cancer2_cell(LL)=[];
Perm_Cancer2 = randperm(length(Cancer2_AUC));
Cancer2_AUC=Cancer2_AUC(Perm_Cancer2,:);
Cancer2_gene=Cancer2_gene(Perm_Cancer2,:);
Cancer2_AUC=[Cancer2_AUC 2*ones(length(Cancer2_AUC),1)];

Ind2=[];
for ii=1:length(Cancer1_cell)
    if strmatch(Cancer1_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind2=[Ind2; strmatch(Cancer1_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind2=[Ind2; 0];
    end
end
L = Ind2(Ind2~=0);
Cancer1_gene=Gene_Expression_Data(:,L)';
LL=find(Ind2==0);
Cancer1_AUC(LL)=[];
Cancer1_cell(LL)=[];
Perm_Cancer1 = randperm(length(Cancer1_AUC));
Cancer1_AUC=Cancer1_AUC(Perm_Cancer1,:);
Cancer1_gene=Cancer1_gene(Perm_Cancer1,:);
Cancer1_AUC=[Cancer1_AUC ones(length(Cancer1_AUC),1)];
%% Modeling and Prediction
n_tree=500;
mtree=10;
min_leaf=2;
Sample_size_Cancer1=size(Cancer1_AUC,1);
Sample_size_Cancer2=size(Cancer2_AUC,1);
Test_size=10;
combined_gene_train=[Cancer1_gene(1:Sample_size_Cancer1-Test_size,:); Cancer2_gene(1:Sample_size_Cancer2-Test_size,:)];
combined_AUC_train=[Cancer1_AUC(1:Sample_size_Cancer1-Test_size,:); Cancer2_AUC(1:Sample_size_Cancer2-Test_size,:)];
combined_gene_test=[Cancer1_gene(Sample_size_Cancer1-Test_size+1:end,:); Cancer2_gene(Sample_size_Cancer2-Test_size+1:end,:)];
combined_AUC_test=[Cancer1_AUC(Sample_size_Cancer1-Test_size+1:end,:); Cancer2_AUC(Sample_size_Cancer2-Test_size+1:end,:)];

% [ranked_combined2, ~]=relieff(combined_gene_train2,combined_AUC_train2, 10);
model_combined = build_forest(combined_gene_train,combined_AUC_train, n_tree, mtree, 1,3, min_leaf, 1);

Y_hat_combined_modified = [forest_predict(combined_gene_test(1:Test_size,:),model_combined,ones(Test_size,1));...
                           forest_predict(combined_gene_test(Test_size+1:end,:),model_combined,2*ones(Test_size,1))];
                       
Result=[mean((combined_AUC_test(:,1)-Y_hat_combined_modified(:,1)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_modified(:,1)))...
     mean((combined_AUC_test(:,1)-Y_hat_combined_modified(:,2)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_modified(:,2)))]    
% %% Fit discriminant analysis classifier
% MdlLinear = fitcdiscr(combined_gene_train,combined_AUC_train(:,2));
% LDA_classify=predict(MdlLinear,combined_gene_test);
% %% Fit binary classification decision tree for multiclass classification
% Decision_Tree = fitctree(combined_gene_train,combined_AUC_train(:,2),'MinLeafSize',min_leaf);
% Decision_Tree_classify=predict(Decision_Tree,combined_gene_test);
% 
% Classify=[Y_hat_combined_modified(:,3) LDA_classify Decision_Tree_classify]