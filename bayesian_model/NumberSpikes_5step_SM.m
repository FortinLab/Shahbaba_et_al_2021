function number_spikes_perTrial = NumberSpikes_5step_SM(spiking_values_wo_gaussian,non_zero,first_step,second_step) %, superchris_sorted_terada

number_spikes = zeros(size(first_step,2),1);
number_spikes_all_neurons =zeros(size(spiking_values_wo_gaussian,1), size(first_step,2)) ; 
for n = 1:size(spiking_values_wo_gaussian,1)
    for i = 1:size(first_step,2)
        first_index = first_step(:,i);
        second_index = second_step(:,i);
        number_spikes(i,:) = sum(spiking_values_wo_gaussian(n,first_index:second_index));
    end
    number_spikes_all_neurons(n,:)  = (number_spikes);
end

number_spikes_perTrial = number_spikes_all_neurons(non_zero,:); 

end