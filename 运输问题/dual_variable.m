function [ udual,vdual ] = dual_variable(A,price)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
s=size(A);
m=s(1);
n=s(2);
udual=zeros(1,m);
vdual=zeros(1,n);
udual(:,:)=inf;
vdual(:,:)=inf;
udual(1)=1;
while 1
    flag=1;
    for i=1:m
        if udual(i)==inf
            flag=0;
        end
    end
    for j=1:n
        if vdual(j)==inf
            flag=0;
        end
    end
    if flag
        break
    end
    for i= 1:m
        for j=1:n
            if A(i,j)==inf|| (udual(i)<inf && vdual(j)<inf)
                continue
            end
            if udual(i)<inf
                vdual(j)=price(i,j)-udual(i);
            end
            if vdual(j)<inf
                udual(i)=price(i,j)-vdual(j);
            end
        end
    end
end
            



end

