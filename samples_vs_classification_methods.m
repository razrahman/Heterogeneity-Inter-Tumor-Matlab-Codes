clear; clc; close all;

HARF=[22 19 14 12 9 8 8];
DT=[30 27 23 19 19 13 12];
LDA=[25 20 17 10 8 6 4];
KNN=[26 19 16 11 9 6 5];
figure
plot(20:20:140,HARF,'-g*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','c','MarkerFaceColor',[0.5,0.5,0.5]); hold on
plot(20:20:140,DT,'--bs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','m','MarkerFaceColor',[0.5,0.5,0.5])
plot(20:20:140,LDA,':rd','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
% plot(20:20:140,KNN)
xticks([20 40 60 80 100 120 140])
xticklabels({'(10,10)','(20,20)','(30,30)','(40,40)','(50,50)','(60,60)','(69,74)'})
%yticks([-1 -0.5 0 0.5 1])
xlabel('Number of samples of cancer types HLT & Lung')
ylabel('Percentage of misclassification')
axis([10 150 0 40])
legend('HARF','Decision Tree','Linear Discriminant Analysis')