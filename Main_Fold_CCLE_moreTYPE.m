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
Cancer1_AUC=[num(cancer_cell_number{1},12)/8; num(cancer_cell_number{10},12)/8]; 
Cancer1_cell=[txt(cancer_cell_number{1}+1,1); txt(cancer_cell_number{10}+1,1)];
Cancer2_AUC=[num(cancer_cell_number{8},12)/8];
Cancer2_cell=[txt(cancer_cell_number{8}+1,1)];
%%
Ind1=[];
for ii=1:length(Cancer1_cell)
    if strmatch(Cancer1_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind1=[Ind1; strmatch(Cancer1_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind1=[Ind1; 0];
    end
end
L = Ind1(Ind1~=0);
Cancer1_gene=Gene_Expression_Data(:,L)';
LL=find(Ind1==0);
Cancer1_AUC(LL)=[]; Cancer1_cell(LL)=[];
Cancer1_AUC=[Cancer1_AUC ones(length(Cancer1_AUC),1)];

Ind2=[];
for ii=1:length(Cancer2_cell)
    if strmatch(Cancer2_cell(ii), Gene_Expression_Cell_lines)~=0
        Ind2=[Ind2; strmatch(Cancer2_cell(ii), Gene_Expression_Cell_lines)];
    else
        Ind2=[Ind2; 0];
    end
    
end
L = Ind2(Ind2~=0);
Cancer2_gene=Gene_Expression_Data(:,L)';
LL=find(Ind2==0);
Cancer2_AUC(LL)=[]; Cancer2_cell(LL)=[];
Cancer2_AUC=[Cancer2_AUC 2*ones(length(Cancer2_AUC),1)];
%% Modeling and Prediction
n_tree=100; mtree=10; F=4; min_leaf=2;
% Sample_size_Cancer1=size(Cancer1_AUC,1);
% Sample_size_Cancer2=size(Cancer2_AUC,1);
finalX=[Cancer1_gene; Cancer2_gene];
finalY=[Cancer1_AUC; Cancer2_AUC];
[Xtrain,Xtest,Ytrain,Ytest,FoldedIndex]=CreateFoldedDataMRF(finalX,finalY,F);
Y_hat_combined=[]; Y_hat_combined_LDA=[]; Y_hat_combined_DT=[];
combined_AUC_test=[];
for FF=1:F
    combined_gene_train=Xtrain{FF};
    combined_AUC_train=Ytrain{FF};
    combined_gene_test=Xtest{FF};
    combined_AUC_test=[combined_AUC_test; Ytest{FF}];
    
%     [ranked_combined, ~]=relieff(combined_gene_train,combined_AUC_train(:,1), 10); % (:,ranked_combined(1:500))
    model_combined = build_forest(combined_gene_train,combined_AUC_train, n_tree, mtree, 1,1, min_leaf, 1);
    
    Y_hat_combined = [Y_hat_combined; forest_predict(combined_gene_test,model_combined,[])];
        %% Fit discriminant analysis classifier
    MdlLinear = fitcdiscr(combined_gene_train,combined_AUC_train(:,2));
    LDA_classify=predict(MdlLinear,combined_gene_test);
    Y_hat_combined_LDA = [Y_hat_combined_LDA; forest_predict(combined_gene_test,model_combined,LDA_classify)];
    %% Fit binary classification decision tree for multiclass classification
    Decision_Tree = fitctree(combined_gene_train,combined_AUC_train(:,2),'MinLeafSize',min_leaf);
    Decision_Tree_classify= predict(Decision_Tree,combined_gene_test);
    Y_hat_combined_DT = [Y_hat_combined_DT; forest_predict(combined_gene_test,model_combined,Decision_Tree_classify)];
end
Result=[mean((combined_AUC_test(:,1)-Y_hat_combined(:,1)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined(:,1)))...
        mean((combined_AUC_test(:,1)-Y_hat_combined(:,2)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined(:,2)))]
Result_LDA=[mean((combined_AUC_test(:,1)-Y_hat_combined_LDA(:,1)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_LDA(:,1)))...
        mean((combined_AUC_test(:,1)-Y_hat_combined_LDA(:,2)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_LDA(:,2)))]
Result_DT=[mean((combined_AUC_test(:,1)-Y_hat_combined_DT(:,1)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_DT(:,1)))...
        mean((combined_AUC_test(:,1)-Y_hat_combined_DT(:,2)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined_DT(:,2)))]
    
Classify=[combined_AUC_test(:,2) Y_hat_combined(:,3) Y_hat_combined_LDA(:,3) Y_hat_combined_DT(:,3)];
for i=1:3
    Misclassification(i,:)=numel(find(Classify(:,1)~=Classify(:,i+1)));
end