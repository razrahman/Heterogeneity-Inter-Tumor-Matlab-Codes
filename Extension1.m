CC1=forest_predict(combined_gene_test(1:Test_size,:),model_combined,ones(Test_size,1));
CC2=forest_predict(combined_gene_test(Test_size+1:end,:),model_combined,2*ones(Test_size,1));
Tree_class1=2*ones(Test_size,n_tree);
Tree_class2=1*ones(Test_size,n_tree);
for ii=1:Test_size
    Tree_class1(ii,CC1{ii+4})=1;
    Tree_class2(ii,CC2{ii+4})=2;
end
for jj=1:n_tree
    C1_sum(jj)=length(find(1==Tree_class1(:,jj)));
    C2_sum(jj)=length(find(2==Tree_class2(:,jj)));
end
Trees=[5 10 20 30 50 75 100 200 300 500];
for tt=1:length(Trees)
    LL=nchoosek(1:Trees(tt),2);
    Dependency1=[];
    for kk=1:length(LL)
        Merge=[Tree_class1(:,LL(kk,:)); Tree_class2(:,LL(kk,:))];
        Dependency1(kk,:)=corr(Merge(:,1),Merge(:,2));
    end
    Dependency=[LL Dependency1];
    Res(tt,:)=[mean(Dependency1) var(Dependency1)];
end