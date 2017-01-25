function [TrainingData,TestingData,OriginalTrain,OriginalTest,FoldedIndex]...
    =CreateFoldedDataResub(ParcentEC50,Original)

FoldedIndex{1}=1:size(Original,1);

    TrainingData{1}=ParcentEC50;
    TestingData{1}=ParcentEC50;
    OriginalTrain{1}=Original;
    OriginalTest{1}=Original;
