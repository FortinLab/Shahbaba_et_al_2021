function mean_firing_all_neurons = fi_SM(spiking_values_wo_gaussian,first_step,second_step,time_bin) 

sum_spikes = zeros(size(first_step,2),1);
sum_spikes_all_n =zeros(size(spiking_values_wo_gaussian,1),size(first_step,2)) ; 
for n = 1:size(spiking_values_wo_gaussian,1)
    for i = 1:size(first_step,2)
        first_index = first_step(:,i);
        second_index = second_step(:,i);
        sum_spikes(i,:) = sum(spiking_values_wo_gaussian(n,first_index:second_index));
    end
    sum_spikes_all_neurons(n,:)  = (sum_spikes);
end


time_bin_seconds = [time_bin/1000];
time_bin_seconds = repmat(time_bin_seconds,size(sum_spikes_all_neurons,1),size(sum_spikes_all_neurons,2)); 
mean_firing_all_neurons = sum_spikes_all_neurons./time_bin_seconds;  

end