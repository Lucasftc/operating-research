function [judge] = check_sigma(udual,vdual,price)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
s=size(price);
m=s(1);
n=s(2);
judge=zeros(m,n);
for i=1:m
    for j=1:n
        judge(i,j)=price(i,j)-udual(i)-vdual(j);
    end
end

end

