function [ totalcost ] = calculate( A,price)
%UNTITLED6 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
totalcost=0;
s=size(A);
m=s(1);
n=s(2);
for i=1:m
    for j=1:n
        if A(i,j)<inf && price(i,j)>0
            totalcost=totalcost+A(i,j)*price(i,j);
        end
    end
end


end

