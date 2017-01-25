clear; clc; close all;
Tree_sizes=[5     10    20    50   100  200  500];
RF=        [20    14    9     6    5    4    4]/68*100;
Bayes=     [16.78 13.19 11.16 6.93 4.00 1.54 2.99e-4];
plot(Tree_sizes,[RF' Bayes'])
xlabel('Number of trees')
ylabel('Percentage of Misclassification')
legend('HARF Misclassification Error','Bayes Error')
% title('')