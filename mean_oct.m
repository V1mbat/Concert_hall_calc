function reduction_oct = mean_oct(f_c, reduction)

temp = zeros(1,length(f_c));
for i = 1:length(f_c)
    f_up = f_c(i)*sqrt(2);
    f_low = f_c(i)/sqrt(2);
    
    temp(i) = mean(reduction(round(f_low):round(f_up)));
    
end
reduction_oct = temp;
end