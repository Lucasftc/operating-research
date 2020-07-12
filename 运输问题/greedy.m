function [ A ] = greedy( A,price,prod_copy,sell_copy )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
A(:,:)=inf;
prod=prod_copy;
sell=sell_copy;
while max(prod)>0 && max(sell)>0 
    [minvalue,row]=min(price);
    [minvalue,col]=min(minvalue);
    row=row(col);
    if prod(row)<sell(col)
        A(row,col)=prod(row);
        sell(col)=sell(col)-prod(row);
        prod(row)=0;
        price(row,:)=inf;
    else if prod(row)>sell(col)
            A(row,col)=sell(col);
            prod(row)=prod(row)-sell(col);
            sell(col)=0;
            price(:,col)=inf;
        else
            A(row,col)=prod(row);
            price(:,col)=inf;
            price(row,:)=inf;
            if(prod(row)==prod_copy(row) && sell(col)==sell_copy(col))
                for j=1 : length(sell)
                    if A(row,j)==-1
                        A(row,j)=0
                        break
                    end
                end
            end
            prod(row)=0;
            sell(col)=0;
            
        end
    end
end
return 

                
            
    




end

