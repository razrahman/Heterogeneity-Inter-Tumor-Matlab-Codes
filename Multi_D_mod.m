function D = Multi_D_mod(y,V_inv,Command)


       ybar1 = sum(y,1)/size(y,1);

       if Command==1 % RF        
            ybar = sum(y)/size(y,1);
           D = sum((y - ybar).^2);
            
        elseif Command==2 %%Using MRF
            yhat = y - ybar1(ones(size(y,1), 1), :); % yhat = y - repmat(ybar,size(y,1),1);
            D = sum(diag((yhat*V_inv*yhat')));
       end
    


