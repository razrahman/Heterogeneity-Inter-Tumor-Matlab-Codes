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
Cancer2_AUC=[Cancer2_AUC -1*ones(length(Cancer2_AUC),1)];

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
T=500;%[5 10 20 30 50 75 100 200 300 500];
for TT=1:length(T)
    n_tree=T(TT);%10;
    mtree=10;
    F=4;
    min_leaf=2;
    Sample_size_Cancer1=size(Cancer1_AUC,1);
    Sample_size_Cancer2=size(Cancer2_AUC,1);
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
        model_combined = build_forest(combined_gene_train,combined_AUC_train, n_tree, mtree, 1,3, min_leaf, 1);
        Tree_class=repmat(-Ytest{FF}(:,2),1,n_tree);
        Temp=forest_predict(combined_gene_test,model_combined,[]);
        for ii=1:size(Ytest{FF},1)
            Tree_class(ii,Temp{ii+4})=Ytest{FF}(ii,2);
        end
        Y_hat_combined = [Y_hat_combined; [Temp{1:4}]];
    end
    Result=[mean((combined_AUC_test(:,1)-Y_hat_combined(:,1)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined(:,1)))...
        mean((combined_AUC_test(:,1)-Y_hat_combined(:,2)).^2) mean(abs(combined_AUC_test(:,1)-Y_hat_combined(:,2)))]
    Y_hat_combined(:,5)=combined_AUC_test(:,2);
    Y_hat_combined(:,6)=Y_hat_combined(:,4).*(Y_hat_combined(:,5)==Y_hat_combined(:,3))...
        +(10-Y_hat_combined(:,4)).*(Y_hat_combined(:,5)~=Y_hat_combined(:,3));
    Classify=[combined_AUC_test(:,2) Y_hat_combined(:,3)];% Y_hat_combined_LDA(:,3) Y_hat_combined_DT(:,3)];
    for i=1%:3
        Misclassification(i,:)=numel(find(Classify(:,1)~=Classify(:,i+1)));
    end
    %% Bayes Error
    b=mean(Y_hat_combined(:,6))/n_tree
    s1=0;
    for K=0:n_tree/2-1
        %     s1=s1+nchoosek(T,k)*(1/a^k*(a/(1+a))^T-a^k*(1+a)^(-T));
        s1=s1+nchoosek(n_tree,K)*(  ( (1-b)^K * b^n_tree /b^K) - (b^K * (1-b)^n_tree / (1-b)^K )  );
    end
    Epsilon(TT,1)=1/2*(1-s1);
    Epsilon(TT,2)=Misclassification/size(combined_AUC_test,1);
end