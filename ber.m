function [rate,count,total] = ber(data1,data2)

total = size(data1,1);
count = sum(double(data1 ~= data2));
rate = count/total;

end