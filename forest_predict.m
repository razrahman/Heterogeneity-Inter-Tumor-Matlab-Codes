function Y_hat = forest_predict(XX_hat, Model, Class)

forest = Model.forest;
n_tree = size(forest,2);
variable_num = Model.variable_num;

Alpha=ones(n_tree,1)/n_tree;
Y_hat = zeros(size(XX_hat,1),variable_num);
Y_hat_X = zeros(size(XX_hat,1),variable_num);
for i = 1:n_tree
    t = forest{1,i};
    Y_predict = single_tree_predict(XX_hat, t, variable_num);
    
        parfor j = 1:variable_num
            Y_hat_X(:,j) = Y_hat_X(:,j) + Alpha(i)*Y_predict(:,j);
        end
    Y_predict_class(:,i)=Y_predict(:,end);
    Y_predict_hat{i}=Y_predict(:,1:end-1);
end

for ii=1:size(XX_hat,1)
    if isempty(Class)
        Class_test(ii,:)=mode(Y_predict_class(ii,:));
        Ind_class{ii}=find(Class_test(ii,:)==Y_predict_class(ii,:));
    elseif ~isempty(Class)
        Class_test(ii,:)=Class(ii,:);
        Ind_class{ii}=find(Class(ii,:)==Y_predict_class(ii,:));
    end
    Class_test_number(ii,:)=length(Ind_class{ii});
    Temp=cell2mat(Y_predict_hat(Ind_class{ii}));
    for kk=1:variable_num
        Y_hat(ii,kk)=mean(Temp(ii,kk:variable_num:size(Temp,2)));
    end
end

Y_hat = [Y_hat Y_hat_X Class_test Class_test_number];% Ind_class];



