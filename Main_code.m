clear;
close;

addpath(genpath('biosig'));
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

% pwelch_baseline_avg=zeros(16,130);

for i=1:16
pwelch_baseline_1channel=pwelch(epoch_baseline.data(:,i,:), 0.5*epoch_baseline.fs, 0.5*0.5*epoch_baseline.fs);

pwelch_baseline_avg{i}=mean(pwelch_baseline_1channel, 2);

plot(10*log10(pwelch_baseline_avg{1,i}));
hold on;
end


