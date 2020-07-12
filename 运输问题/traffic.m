function [ A,totalcost ] = traffic( price,prod,sell )
%traffic 表上作业法解决方案
%
[m,n]=size(price);
A=zeros(length(prod),length(sell));
A=greedy(A,price,prod,sell);
[udual,vdual]=dual_variable(A,price);
judge=check_sigma(udual,vdual,price);
[min1,row1]=min(judge);
[min2,col1]=min(min1);
while(min2)<0
    i=row1(col1);
    j=col1;
    [result,point1,point2]=find_close_path(A,i,j);
    A=adjust_distribution(A,point1,point2);
    [udual,vdual]=dual_variable(A,price);
    judge=check_sigma(udual,vdual,price);
    [min1,row1]=min(judge);
    [min2,col1]=min(min1);
end
totalcost=calculate(A,price);
return 

end

