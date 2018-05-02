function [ struct_window ] = epoch_window( session,Align_Event,timebEvent,timeaEvent,windowlength,overlap)
  
struct_window.Event_pos=session.event.position(session.event.name==Align_Event);
struct_window.sampleRate=512;
 
samplesBeforeAlignEvents=abs(timebEvent)/overlap; % they should be both multiply by the sample frequences
samplesAfterAlignEvents=abs((timeaEvent-windowlength))/overlap;
struct_window.labels=session.event.name(session.event.name==Align_Event);
 
struct_window.samples=[];
 
 for i=1:size(struct_window.labels,1)
 
 A=session.data(1:16,:,(struct_window.Event_pos(i)-samplesBeforeAlignEvents):struct_window.Event_pos(i));
 %B=session.data(1:16,:,struct_window.Event_pos(i));
 C=session.data(1:16,:,(struct_window.Event_pos(i)+1:struct_window.Event_pos(i)+samplesAfterAlignEvents));
 struct_window.samples=cat(3,struct_window.samples,A,C);
 end
 
 struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*(samplesBeforeAlignEvents+samplesAfterAlignEvents+1),1)*Align_Event;
  end