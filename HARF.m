function Final_Sensitivity=HARF(Cancer1_gene,Cancer1_Sensitivity,Cancer2_gene,Cancer2_Sensitivity,Test_gene,n_tree,mtree,min_leaf)
% HARF performs Heterogeneity Aware Random Forest
%     Input:
%         Cancer1_gene
%         Cancer1_Sensitivity
%         Cancer2_gene
%         Cancer2_Sensitivity
%         Test_gene
%         n_tree
%         mtree
%         min_leaf
% %   Output:
%         Final_Sensitivity: Predicted sensitivities using Heterogenity Aware Random Forest (HARF), where 1st column values corresponds to classification
%         using HARF and 2nd column corresponds to classification using Linear Discriminant Analysis
%
% %   Example:
%        Cancer1_gene=randn(100,10000);
%        Cancer1_Sensitivity =.25.*rand(100,1);
%        Cancer2_gene=randn(100,10000);
%        Cancer2_Sensitivity =.5.*rand(100,1);
%        Test_gene=randn(50,10000);
%        Final_Sensitivity=HARF(Cancer1_gene,Cancer1_Sensitivity,Cancer2_gene,Cancer2_Sensitivity,Test_gene,[],[],[]);


Cancer1_Sensitivity=[Cancer1_Sensitivity ones(length(Cancer1_Sensitivity),1)];
Cancer2_Sensitivity=[Cancer2_Sensitivity 2*ones(length(Cancer2_Sensitivity),1)];

finalX=[Cancer1_gene; Cancer2_gene];
finalY=[Cancer1_Sensitivity; Cancer2_Sensitivity];

if isempty(n_tree)
    n_tree=size(finalX,1)/2;
end
if isempty(mtree)
    mtree=10;
end
if isempty(min_leaf)
    min_leaf=3;
end
%% Modeling and Prediction
[ranked_combined, ~]=relieff(finalX,finalY(:,1), 10); % (:,ranked_combined(1:500))
combined_gene_train1=finalX(:,ranked_combined(1:500));
combined_gene_test1=Test_gene(:,ranked_combined(1:500));
model_combined = build_forest(finalX,finalY, n_tree, mtree, 1,1, min_leaf, 1);

HARF_sensitivity = forest_predict(Test_gene,model_combined,[]);

% % Fit discriminant analysis classifier
MdlLinear = fitcdiscr(combined_gene_train1,finalY(:,2));
LDA_classify=predict(MdlLinear,combined_gene_test1);
LDA_sensitivity = forest_predict(Test_gene,model_combined,LDA_classify);

Final_Sensitivity=[HARF_sensitivity(:,1) LDA_sensitivity(:,1)];
