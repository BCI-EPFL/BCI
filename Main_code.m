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


a=epoch_struct(session,h,200,2,3);



