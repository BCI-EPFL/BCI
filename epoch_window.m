function [ struct_window ] = epoch_window( session,Align_Event,timebEvent,timeaEvent,windowlength,overlap)
  


 struct_window.labels=session.event.name(session.event.name==Align_Event);
 struct_window.Event_pos=session.event.position(session.event.name==Align_Event);
 
 samplesBeforeAlignEvents=abs(timebEvent)/overlap;
 samplesAfterAlignEvents=abs((timeaEvent-windowlength))/overlap;
 
 struct_window.samples=[];
 
 for i=1:size(struct_window.labels,1)
 
 A=session.data(1:16,:,(struct_window.Event_pos(i)-samplesBeforeAlignEvents):(struct_window.Event_pos(i)-1));
 %B=session.data(1:16,:,struct_window.Event_pos(i));
 C=session.data(1:16,:,struct_window.Event_pos(i):(struct_window.Event_pos(i)+samplesAfterAlignEvents));
       struct_window.samples=cat(3,struct_window.samples,A,C);
 end
 

  end
