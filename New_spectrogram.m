function [final_avg, t, f] = New_spectrogram(thisEpoch1, thisEpoch2,num_channel,win,noverlap,freq)

k=0;

for i=1:size(thisEpoch1.samples,3)
    
    [~,f,t,p_baseline] = spectrogram(thisEpoch1.samples(num_channel,:,i), 16, 'power');
    [~,f,t,p_MI] = spectrogram(thisEpoch2.samples(num_channel,:,i),16, 'power');
    
    power_avg_baseline = mean(p_baseline');
    
    power_ratio{i} = p_MI./(power_avg_baseline');
    
    if rem(i,33)==0
    k=k+1;
    avg{k} = sum(cat(3,power_ratio{(k-1)*33+1:end}),3)./33;
    
    end
    
end
    
final_avg=sum(cat(3,avg{:}),3)./(size(thisEpoch1.samples,3));
    

end