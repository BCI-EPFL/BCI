function [ struct_window ] = epoch_window( session,Align_Event,timebEvent,timeaEvent,windowlength,overlap)
  
struct_window.Event_pos=session.event.position(session.event.name==Align_Event);
struct_window.sampleRate=512;

if (timeaEvent ~=0 )
    struct_window.Event_pos=struct_window.Event_pos+windowlength/overlap; %we have to shift because "backward" has the 0 (i.e. the event) on the right


%     if timebEvent==0
%     samplesBeforeAlignEvents=abs(timebEvent)/overlap; % they should be both multiply by the sample frequences
%     samplesAfterAlignEvents=abs((timeaEvent-windowlength))/overlap;
%     %samplesAfterAlignEvents=abs(timeaEvent)/overlap;
%     struct_window.labels=session.event.name(session.event.name==Align_Event);
%  
%     struct_window.samples=[];
%  
%     for i=1:size(struct_window.labels,1)
%  
%         A=session.data(1:16,:,(struct_window.Event_pos(i)-samplesBeforeAlignEvents):struct_window.Event_pos(i));
%         C=session.data(1:16,:,((struct_window.Event_pos(i)+1):(struct_window.Event_pos(i)+samplesAfterAlignEvents)));
%         %C=session.data(1:16,:,struct_window.Event_pos(i):(struct_window.Event_pos(i)+samplesAfterAlignEvents));
%         struct_window.samples=cat(3,struct_window.samples,A,C);
%         %struct_window.samples=cat(3,struct_window.samples,C);
%     end
%     struct_window.duration=samplesAfterAlignEvents+samplesBeforeAlignEvents+1;
%     struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*(samplesBeforeAlignEvents+samplesAfterAlignEvents+1),1)*Align_Event;
% end
%     if timebEvent>0
    %samplesAfterAlignEvents=abs((timeaEvent-timebEvent))/overlap;
   
    struct_window.labels=session.event.name(session.event.name==Align_Event);
    struct_window.samples=[];
    
    for i=1:size(struct_window.labels,1)
 
       % A=session.data(1:16,:,(struct_window.Event_pos(i)-samplesBeforeAlignEvents):struct_window.Event_pos(i));
        C=session.data(1:16,:,((struct_window.Event_pos(i)+timebEvent/overlap):(struct_window.Event_pos(i)+(timeaEvent-windowlength)/overlap)));
        %C=session.data(1:16,:,struct_window.Event_pos(i):(struct_window.Event_pos(i)+samplesAfterAlignEvents));
        %struct_window.samples=cat(3,struct_window.samples,A,C);
        struct_window.samples=cat(3,struct_window.samples,C);
    end
    struct_window.duration=(timeaEvent-timebEvent-windowlength)/overlap+1;
    struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*struct_window.duration,1)*Align_Event;
    
 
    else  
     samplesBeforeAlignEvents=abs(timebEvent+windowlength)/overlap; 
     samplesAfterAlignEvents=abs((timeaEvent))/overlap;
     struct_window.labels=session.event.name(session.event.name==Align_Event);
 
    struct_window.samples=[];
 
    for i=1:size(struct_window.labels,1)
 
        A=session.data(1:16,:,(struct_window.Event_pos(i)-samplesBeforeAlignEvents):struct_window.Event_pos(i));
        %C=session.data(1:16,:,((struct_window.Event_pos(i)+1):(struct_window.Event_pos(i)+samplesAfterAlignEvents)));
        %C=session.data(1:16,:,struct_window.Event_pos(i):(struct_window.Event_pos(i)+samplesAfterAlignEvents));
        %struct_window.samples=cat(3,struct_window.samples,A,C);
        struct_window.samples=cat(3,struct_window.samples,A);
    end
    struct_window.duration=samplesAfterAlignEvents+samplesBeforeAlignEvents+1;
    struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*(samplesBeforeAlignEvents+samplesAfterAlignEvents+1),1)*Align_Event;

     
     
     % 
 %struct_window.duration=samplesAfterAlignEvents+samplesBeforeAlignEvents+1;
 %struct_window.duration=samplesAfterAlignEvents+1;
 %struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*(samplesBeforeAlignEvents+samplesAfterAlignEvents+1),1)*Align_Event;
 %struct_window.labels=ones(size(session.event.name(session.event.name==Align_Event),1)*(samplesAfterAlignEvents+1),1)*Align_Event;
  end
