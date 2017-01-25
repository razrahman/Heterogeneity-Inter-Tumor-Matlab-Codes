clear; clc; close all;
L=20;
W=100;
s2=zeros(W,L);
alpha=zeros(W,L);
Y=zeros(1,L);
Z=zeros(1,W);

for U=1:W
    for j=1:L
        
        a=1+(j-1)/5;
        
        s1=0;
        T=2*U+1;
        
        for i=0:(T+1)/2-1
            s1=s1+nchoosek(T,i)*(1/a^i*(a/(1+a))^T-a^i*(1+a)^(-T));
        end
        
        
        
        s2(U,j)=1/2*(1-s1);
        %alpha(U,j)=a;
        Y(j)=a;
        Z(U)=T;
        
    end
end

mesh(Y,Z,s2)
%plot(alpha,s2)