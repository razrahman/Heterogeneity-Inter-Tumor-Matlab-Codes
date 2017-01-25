function [TrainingData,TestingData,OriginalTrain,OriginalTest,FoldedIndex]=CreateFoldedDataMRF(ParcentEC50,Original,F)

DrugNumber=size(Original,1);
Index=[1:DrugNumber];


for i=1:F-1
   FoldedIndex{i}=randsample(Index,floor(DrugNumber/F));
   Index=setdiff(Index,FoldedIndex{i});
   if i==F-1
       FoldedIndex{i+1}=Index;
   end
end


for Fold=1:F
    TestingIndex{Fold}=FoldedIndex{Fold};
    TrainingIndex{Fold}=setdiff([1:DrugNumber],TestingIndex{Fold});
    TrainingData{Fold}=ParcentEC50(TrainingIndex{Fold},:);
    TestingData{Fold}=ParcentEC50(TestingIndex{Fold},:);
    OriginalTrain{Fold}=Original(TrainingIndex{Fold},:);
    OriginalTest{Fold}=Original(TestingIndex{Fold},:);
end

