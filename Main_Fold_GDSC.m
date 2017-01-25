% clear; clc;
% load('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\GDSC\Raw_GDSC_X.mat')
% Gene_Expression_Cell_lines=Raw_GDSC_X(2:end,1);
% Gene_Expression_Cell_lines=strcat(Gene_Expression_Cell_lines,'_GDSC');
% Gene_Expression_Data = cell2mat(Raw_GDSC_X(2:end,2:end));
% [Gene_Expression_Cell_lines,Ind_temp]=unique(Gene_Expression_Cell_lines);
% Gene_Expression_Data=Gene_Expression_Data(Ind_temp,:);
%%
[num,txt,raw]=xlsread('C:\Users\razrahma\Dropbox\Data CCLE and GDSC\GDSC\GDSC_Cell_AUC.xlsx');
Cell_Cancer_type=txt(7:end,3); 
Cell_name=strcat(txt(7:end,1),'_GDSC');
Cancer_cell_type={'blood','nervous_system','soft_tissue','bone','skin','urogenital_system',...
    'lung','pancreas','aero_digestive_tract','breast','kidney',...
    'digestive_system','thyroid'};
Drug=64;
Drug1_AUC=num(8:end,Drug+4);
Drug_AUC=Drug1_AUC(~any(isnan(Drug1_AUC),2),:);
Cell_Cancer_type_Drug=Cell_Cancer_type(~any(isnan(Drug1_AUC),2),:);
Cell_name_Cancer_type_Drug=Cell_name(~any(isnan(Drug1_AUC),2),:);
cancer_cell_number=[]; cancer_cell_name=[];
for j=1:length(Cancer_cell_type)
    cancer_cell_number{j,1}=strmatch(Cancer_cell_type(j), Cell_Cancer_type_Drug);
    AUC_cancer_type{j}=1-Drug_AUC(cancer_cell_number{j,1},1);
    cancer_cell_name{j}=Cell_name_Cancer_type_Drug(cancer_cell_number{j,1},1);
end
Cancer1_AUC =[AUC_cancer_type{1}];%  
Cancer1_cell=[cancer_cell_name{1}];%
Cancer2_AUC=AUC_cancer_type{5};%[AUC_cancer_type{7}; AUC_cancer_type{10}];% 
Cancer2_cell=cancer_cell_name{5};%[cancer_cell_name{7}; cancer_cell_name{10}];%  
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
Cancer1_gene=Gene_Expression_Data(L,:);
LL=find(Ind1==0);
Cancer1_AUC(LL)=[];
Cancer1_cell(LL)=[];
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
Cancer2_gene=Gene_Expression_Data(L,:);
LL=find(Ind2==0);
Cancer2_AUC(LL)=[];
Cancer2_cell(LL)=[];
Cancer2_AUC=[Cancer2_AUC 2*ones(length(Cancer2_AUC),1)];
%% Modeling and Prediction 
n_tree=100; mtree=10; F=3; min_leaf=2;
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
    
    [ranked_combined, ~]=relieff(combined_gene_train,combined_AUC_train(:,1), 10); % (:,ranked_combined(1:500))
    model_combined = build_forest(combined_gene_train(:,ranked_combined(1:1000)),combined_AUC_train, n_tree, mtree, 1,1, min_leaf, 1);
    
    Y_hat_combined = [Y_hat_combined; forest_predict(combined_gene_test(:,ranked_combined(1:1000)),model_combined,[])];
        %% Fit discriminant analysis classifier
    MdlLinear = fitcdiscr(combined_gene_train(:,ranked_combined(1:1000)),combined_AUC_train(:,2));
    LDA_classify=predict(MdlLinear,combined_gene_test(:,ranked_combined(1:1000)));
    Y_hat_combined_LDA = [Y_hat_combined_LDA; forest_predict(combined_gene_test(:,ranked_combined(1:1000)),model_combined,LDA_classify)];
    %% Fit binary classification decision tree for multiclass classification
    Decision_Tree = fitctree(combined_gene_train(:,ranked_combined(1:1000)),combined_AUC_train(:,2),'MinLeafSize',min_leaf);
    Decision_Tree_classify= predict(Decision_Tree,combined_gene_test(:,ranked_combined(1:1000)));
    Y_hat_combined_DT = [Y_hat_combined_DT; forest_predict(combined_gene_test(:,ranked_combined(1:1000)),model_combined,Decision_Tree_classify)];
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