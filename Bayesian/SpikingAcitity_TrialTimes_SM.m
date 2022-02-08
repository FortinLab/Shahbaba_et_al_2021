function spiking_acitivity = SpikingAcitity_TrialTimes_SM(ensembleMatrix_unit,ind_times) 

%this includes the spiking acitivity of all neurons during the trial times.

for i=1:length(ensembleMatrix_unit(1,:))
    neuron_activity = ensembleMatrix_unit(:,i);  
    for j = 1:size(ind_times,1)
        neuronActivity_timeWindow(:,j) = neuron_activity(ind_times(j,:));        
    end
    spiking(:,i) = (neuronActivity_timeWindow(:,j)); 
end
spiking_acitivity = transpose(spiking);

end