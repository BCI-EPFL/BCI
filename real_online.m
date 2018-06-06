%{
Inputs: whole signal 4th run (s{1,4}=Session), mu, sigma, classifier, FeaturesIndeces


%}

alpha=0.96;
f=[4:2:40];

%%Divide the run into the (30) trials
TrialStart=Session.event.position(Session.event.name==400);%400=MI initiation
TrialEnd=Session.event.position(Session.event.name==700);%700=Relax phase

Trials=[];
for i=1:(length(TrialStart))
    Trials{i}=Session.data(TrialStart(i):TrialEnd(i),:)';
end



%%Divide each trial into 1-second windows shifted by 0.0625s
%(1 second=512 points; 0.0625 s = 32 points)
Aux=Trials;
Trials=[];
for i=1:length(Aux)
j = 1;
    for t = 512:32:(size(Aux{i},2)-mod(size(Aux{i},2),32)) %we remove the last window (not of 1 second)
        for idxChannels = 1:16
            Trials{i}.Windows(idxChannels,:,j)=Aux{i}(idxChannels,t-512+1:t);
            %bufferedData = testData(idxChannels,t-512+1:t);
        end 
        j= j+1;
    end
end

%% Real pseudoonline
Evidence=[];
for i=1:length(Trials)
   Evidence{i}=zeros(1,size(Trials{i}.Windows,3)+1);
   Evidence{1,i}(1)=0.5; %for each trial we start by 0.5
   Trials{i}.Features=zeros(size(Trials{i}.Windows,3),304);
   
   for j=1:size(Trials{i}.Windows,3)
       %%CAR filtering
       mean_channels=[];
       mean_channels = mean(Trials{i}.Windows(:,:,j));
       Trials{i}.Windows(:,:,j) = Trials{i}.Windows(:,:,j)- mean_channels;
       
       %%pwelch (already reshaped)
       for idxChannels=1:16
            [Trials{i}.Features(j,19*(idxChannels-1)+1:19*idxChannels),~] = pwelch(Trials{i}.Windows(idxChannels,:,j),0.5*512,0.4375*512,f,512);
       end
       for idxChannels=1:16
            Trials{i}.Features(j,19*(idxChannels-1)+1:19*idxChannels) = log(Trials{i}.Features(j,19*(idxChannels-1)+1:19*idxChannels));
       end 
       
       %%normalization
       Trials{i}.Features(j,:)=(Trials{i}.Features(j,:)-mu)./sigma;
       
       %%predict (posterior prob)
       [Trials{i}.Predicted(j),Trials{i}.PosteriorProb(j,:),~]=predict(classifier,Trials{i}.Features(j,FeaturesIndeces));

       
       %%smoothing
       Evidence{1,i}(j+1)=Evidence{1,i}(j)*alpha+(1-alpha)*Trials{i}.PosteriorProb(j,2);
       fprintf('Trial %d Window %d \n',i,j);

   end
   
end

%% Concatenate all the smoothed probabilities
SmoothedTotal=[];
for i=1:length(Trials)
   SmoothedTotal=cat(2,SmoothedTotal,Evidence{1,i}); 
end

figure
sz = 4;
c = [0,0,1];
scatter(1:(30/length(SmoothedTotal)):(31-30/length(SmoothedTotal)),SmoothedTotal,sz,c,'filled')
hold on
cont=1;
for i=1:length(Trials)
    vline(cont,'r','');
    cont=cont+length(Evidence{1,i})*30/length(SmoothedTotal);
end
xlabel('Trials')
ylabel('Smoothed probability')
