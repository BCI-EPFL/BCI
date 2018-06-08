function [pow, freq]=pwelch_for_each_channel(num_channel,epoched_data,f,fs)
    pow=0;
    for i=1:size(epoched_data.data,3)
        [a,b]=pwelch(squeeze(epoched_data.data(:,num_channel,i)),fs,0.5*fs,f,512);
        pow=pow+a;
    end
    pow=pow/i;
    freq=b;
end