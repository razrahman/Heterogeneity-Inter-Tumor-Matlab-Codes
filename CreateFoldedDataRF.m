function [TrainingData,TestingData,OriginalTrain,OriginalTest]=CreateFoldedDataRF(ParcentEC50,Original,F,FoldedIndex)

DrugNumber=size(Original,1);

if F>1
    for Fold=1:F
        TestingIndex{Fold}=FoldedIndex{Fold};
        TrainingIndex{Fold}=setdiff([1:DrugNumber],TestingIndex{Fold});
        TrainingData{Fold}=ParcentEC50(TrainingIndex{Fold},:);
        TestingData{Fold}=ParcentEC50(TestingIndex{Fold},:);
        OriginalTrain{Fold}=Original(TrainingIndex{Fold},:);
        OriginalTest{Fold}=Original(TestingIndex{Fold},:);
    end
else
        Fold=1;
        TrainingIndex{Fold}=FoldedIndex{Fold};
        TestingIndex{Fold}=setdiff([1:DrugNumber],TrainingIndex{Fold});
        TrainingData{Fold}=ParcentEC50(TrainingIndex{Fold},:);
        TestingData{Fold}=ParcentEC50(TestingIndex{Fold},:);
        OriginalTrain{Fold}=Original(TrainingIndex{Fold},:);
        OriginalTest{Fold}=Original(TestingIndex{Fold},:);
end
    