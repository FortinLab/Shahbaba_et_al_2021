function [first_step, second_step] = StepsWindow_SM(time_length,time_bin,rate)

steps = [1:rate:time_length];
add_step = repmat(1,1,size(steps,2)); 
add_step = add_step .* time_bin;
ind = find(steps==(max(steps)-time_bin));  
second_step = steps(1,1:(ind)) + add_step(1,1:(ind)); 
first_step = steps(1,1:(ind)); 

end 


 

