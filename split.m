function [index_left, index_right, which_feature, threshold_feature, DD, IN, indexX] = split(X, Y, index, V_inv, mtree,NN,copula,Command, Alpha, min_leaf)


x = X(index,:);
y = Y(index,:);

f = size(x,2); % number of features
% f_1 = randi([1, f]);
f_1 = sort(randperm(f, mtree));
N = size(x,1); % number of samples
% score = zeros(N-1,f);
min_score = [];
id = cell(N-1,length(f_1));
% threshold = zeros(N-1,f);
for j = 1:length(f_1)
  
    xj = x(:,f_1(j));
    [s, idx] = sort(xj); % sort the features
    
    %for k1 = min_leaf:(N-min_leaf) % calculate D(N)
    for k1 = 1:(N-1) % calculate D(N)
       
        %k=k1-min_leaf+1; %%
        k=k1;
        yk_left = y(idx(1:k1),:);
        yk_right = y(idx(k1+1:end),:);
        D_left= Multi_D_mod(yk_left(:,1:end-1),V_inv,Command);

        D_right = Multi_D_mod(yk_right(:,1:end-1), V_inv,Command);
        
        
        if Command==3
            D1L(j, k)=D_left(1);
            D2L(j, k)=D_left(2);
            D1R(j, k)=D_right(1);
            D2R(j, k)=D_right(2);

            D=D1L(j,k)+D1R(j,k)+Alpha*(D2L(j,k)+D2R(j,k));
            %D=(D1L(j,k))^Alpha+(D1R(j,k))^Alpha+((D2L(j,k))^Alpha+(D2R(j,k))^Alpha);
            
        elseif Command==1 | Command==2
            D=D_left+D_right;
            D1L(j,k)=D_left;
            D2L(j,k)=0;
            D1R(j,k)=D_right;
            D2R(j,k)=0;
        end
            
        if Command==1 | Command==2 | Command==3
            
            
            if j==1&&k==1
                min_score = D;
                which_feature = f_1(1);
                indexX=[j k];
                threshold_feature = (s(k) + s(k+1))/2;
                ij_i = 1;
                ij_j = 1;
            end

            if D < min_score
                min_score = D;
                which_feature = f_1(j); 
                indexX=[j k];
  % which feature determines the split
                threshold_feature = (s(k) + s(k+1))/2; % threshold calculation
                ij_i = k;
                ij_j = j;
            end
        
        end
        
%         score(k,j) = D;
        id{k,j} = idx(1:k1);
      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if Command==3
%         XX=[D1L(:)+D1R(:) D2L(:)+D2R(:)];
%         front=paretoGroup(XX);
%         in1=find(front==1);
%         d1 = abs(XX(in1,1)-XX(in1,2))/sqrt(2);  %%distance from y=x line
%         impPoint=in1(min(find(d1==min(d1))));   %% Taking the point closest to y=x line
%         %impPoint=randsample(in1,1);   %%taking random pareto frontier point
% 
%         Z=zeros(size(front));
%         Z(impPoint)=1;
%         Z=vec2mat(Z, size(D1L,1))';
% 
%         [j k]=find(Z==1);
%         which_feature=f_1(j);
%         threshold_feature = (s(k) + s(k+1))/2;
%         ij_i = k;
%         ij_j = j;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DD={D1L D2L;D1R D2R};

index_left = id{ij_i,ij_j};
index_right = 1:N;
index_right(index_left) =[];

index_left = index(sort(index_left));
index_right = index(sort(index_right));
IN={index_left index_right};

