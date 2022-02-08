function probability_matrix = BayesianAnalysis_SM(mean_firing,number_spikes,first_step,time_bin) 

%the terms below explain my variables. 
%%number_spikes = ni (number of spikes in the time window)
%%mean_firing_rate = fi(time,odors)
%%normalization_factor = C
%%posteriors = product of the first term and second term 

probability_matrix = zeros([size(first_step,2) size(first_step,2)]);

time_window = time_bin/1000; %converting the time window in seconds 

%%computing the first term of the equation:
product_first_term = []; final_first_term = []; 
fact = factorial(number_spikes);
for j = 1:size(first_step,2)   
    for i = 1:size(mean_firing,1)   
        first_term(i,:) = ((time_window*mean_firing(i,:)).^number_spikes(i,j))./fact(i,j);
        if size(first_term,1)>1
            product_first_term = prod(first_term,1);  
        else
            product_first_term = first_term;   
        end
    end    
    final_first_term(j,:) = product_first_term; 
end

%%computing the second term of the equation 
sum_firing_rate_each_neuron = sum(mean_firing,1);
second_term = exp((-1)*time_window*sum_firing_rate_each_neuron);

for i = 1: size(final_first_term,1)
    first_second = final_first_term(i,:).*second_term; 
    posteriors(i,:) = first_second; 
end    

%the line below is for P(time,odor) term. Since the rat spend equal time in
%each time window, I am multiplying the posterior values to 1/time_window
%size. 
posteriors = posteriors.*(1/(size(first_step,2)/4)); 

%probability_matrix = posteriors; 
%%computing the normalization factor. this makes the sum for the
%%probability distributions during each time window 1. 

normalization_factor = 1./sum(posteriors,2);
normalization_factor = repmat(normalization_factor,1,size(first_step,2)); 
probability_matrix = posteriors.*normalization_factor;

end 

