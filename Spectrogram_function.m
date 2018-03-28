function [~,f,t,p] = Spectrogram_function(thisEpoch,win,noverlap,freq,fs,'power')

for i=1:size(thisEpoch.data,3)

    [~,f,t,p] = spectrogram(thisEpoch,win,noverlap,freq,fs,'power');
    [~,f,t,p] = spectrogram(thisEpoch,win,noverlap,freq,fs,'power');

end