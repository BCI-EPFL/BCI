function [pow, freq]=pwelch_decoder(num_channel,epoched_data,f,fs)

for i=1:size(epoched_data.data,3)
    [a,b]=pwelch(squeeze(epoched_data.data(:,num_channel,i)),0.5*fs,0.5*0.5*fs,f,fs);
    pow{i}=a;
    freq=b;
end


end