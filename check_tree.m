function Y_hat = check_tree(Model,Train_data,Train_out,Class, Test_data)

forest = Model.forest;
n_tree = size(forest,2);
variable_num = Model.variable_num;

TT=zeros(size(Train_data,1),n_tree);
Ind_one=find(1==Train_out(:,2));
Ind_two=find(2==Train_out(:,2));
for i = 1:n_tree
    t = forest{1,i};
    
    for ii = 1:size(Train_data,1)
        x = Train_data(ii,:);
    
        leaf_info = predict(x,t);
        TT(ii,i) = mode(leaf_info(:,end));
    end
    
    if (numel(find(TT(Ind_one,i)~=1))/length(Ind_one)*100 >= numel(find(TT(Ind_two,i)~=2))/length(Ind_two)*100)
        TTT(i,:)=2;
    else
        TTT(i,:)=1;
    end
end
%%
Alpha=ones(n_tree,1)/n_tree;
Y_hat = zeros(size(Test_data,1),variable_num);
Y_hat_X = zeros(size(Test_data,1),variable_num);
for i = 1:n_tree
    t = forest{1,i};
    Y_predict = single_tree_predict(Test_data, t, variable_num);
    
        parfor j = 1:variable_num
            Y_hat_X(:,j) = Y_hat_X(:,j) + Alpha(i)*Y_predict(:,j);
        end
    Y_predict_class(:,i)=Y_predict(:,end);
    Y_predict_hat{i}=Y_predict(:,1:end-1);
end

for ii=1:size(Test_data,1)
    Class_test(ii,:)=Class(ii,:);
    Ind_class{ii}=find(Class(ii,:)==TTT);
    
    Class_test_number(ii,:)=length(Ind_class{ii});
    Temp=cell2mat(Y_predict_hat(Ind_class{ii}));
    for kk=1:variable_num
        Y_hat(ii,kk)=mean(Temp(ii,kk:variable_num:size(Temp,2)));
    end
end

Y_hat = [Y_hat Y_hat_X Class_test Class_test_number];% Ind_class];