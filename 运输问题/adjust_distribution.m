function [ A ] = adjust_distribution( A,point1,point2 )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
[p_min,index]=min(point1(:,3));
    point1(:,3)=point1(:,3)-p_min;
    point2(:,3)=point2(:,3)+p_min;

    for k=1:length(point1(:,1))
        A(point1(k,1),point1(k,2))=point1(k,3);
    end

    for k=1:length(point2(:,1))
        A(point2(k,1),point2(k,2))=point2(k,3);
    end
    A(point1(index,1),point1(index,2))=inf;
end

