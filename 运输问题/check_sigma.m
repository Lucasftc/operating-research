function [judge] = check_sigma(udual,vdual,price)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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

