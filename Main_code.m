clear;
close;

addpath(genpath('biosig\t200_FileAccess'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));

load('channel_location_16_10-20_mi');

filename = 's_run1_offlineMIterm_20180703154501.gdf';
[s, h]= sload(filename);



session.fs=h.SampleRate;
session.data=(s)';
session.channels={chanlocs16.labels};
session.Event_type=h.EVENT.TYP;
session.Event_pos=h.EVENT.POS;


epoch_baseline=epoch_struct(session,200,0,3);

%% pwelch for each channel

figure;
for i=1:16
pwelch_baseline_1channel=pwelch(epoch_baseline.data(:,i,:), 0.5*epoch_baseline.fs, 0.5*0.5*epoch_baseline.fs);

pwelch_baseline_avg{i}=mean(pwelch_baseline_1channel, 2);

plot(10*log10(pwelch_baseline_avg{1,i}));
hold on;
end

%% time filtering (non va pero')
[b,a] = butter(2, [5 40]/512/2);
filt_data=filter(b,a,session.data');
session_filt=session;
session_filt.data=filt_data';
filt_epoch_baseline=epoch_struct(session_filt,200,0,3);

%%
figure
for i=1:16
filt_pwelch_baseline_1channel=pwelch(filt_epoch_baseline.data(:,i,:), 0.5*filt_epoch_baseline.fs, 0.5*0.5*filt_epoch_baseline.fs);

filt_pwelch_baseline_avg{i}=mean(filt_pwelch_baseline_1channel, 2);

plot(10*log10(filt_pwelch_baseline_avg{1,i}));
hold on;
end

%% temporal filtering on the raw data
data_filter=zeros(16,330752);

for i=1:size(chanlocs16,2)
data=session.data(i,:);
[b,a]=butter(4,[5 40] /512/2); %bandpass
data_filter(i,:)=filter(b,a,data);
figure;
hold on
plot(data);
hold on
plot(data_filter(i,:));
end 


%% sdp after temporal filtering

session.data=data_filter;
baseline=epoch_struct(session, 200, 0,3);
MI=epoch_struct(session, 400, 0,3);


for i=1:size(chanlocs16,2)
[pwr_b, fre_b]=pwelch(squeeze(epoch_baseline.data(:,i,:)),0.5*session.fs,0.5*0.5*session.fs,500,session.fs);

plot(fre_b,10*log(pwr_b));%10log
hold on
legend('baseline channel 1','baseline channel 2','baseline channel 3','baseline channel 4','baseline channel 5','baseline channel 6','baseline channel 7','baseline channel 8','baseline channel 9','baseline channel 10','baseline channel 11','baseline channel 12','baseline channel 13','baseline channel 14','baseline channel 15','baseline channel 16');
xlabel('frequencies');
ylabel('spectral density');
title('spectral density');
end 
%%
for i=1:size(chanlocs16,2)
[pwr_MI,fre_MI]=pwelch(squeeze(MI.data(:,i,:)),0.5*session.fs,0.5*0.5*session.fs,500,session.fs);
plot(fre_MI,pwr_MI);
hold on 
legend('motor imagery channel 1','motor imagery channel 2','motor imagery channel 3','motor imagery channel 4','motor imagery channel 5','motor imagery channel 6','motor imagery channel 7','motor imagery channel 8','motor imagery channel 9','motor imagery channel 10','motor imagery channel 11','motor imagery channel 12','motor imagery channel 13','motor imagery channel 14','motor imagery channel 15','motor imagery channel 16');
xlabel('frequencies');
ylabel('spectral density');
title('spectral density');
end 


%% spatial filtering on the raw data.CAR
medium=mean(mean(s));
signal_car=zeros(size(s,1),size(s,2));
for i=1:size(s,2)-1
    signal_car(:,i)=s(:,1)-medium;
end

plot(signal_car(:,9))
hold on
plot(s(:,9))
title ('car filter and raw signal');
labels('CAR signal','raw signal');

%%laplacian filter
laplacian=s(:,1:16)*lap;
subplot(3,1,1)
plot(s(:,9))
subplot(3,1,2)
plot(laplacian(:,9))
subplot(3,1,3)
plot(signal_car(:,9))
  

%topoplot






















