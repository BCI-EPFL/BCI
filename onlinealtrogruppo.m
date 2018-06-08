%% ONLINE %%

% parameters
timeFrame = 1; % one second buffer length
nonoverlap = 0.0625; % 62.5ms frame shift

samplingRate=512;
psdWindow = 0.5*samplingRate;
psdNOverlap = 0.25*samplingRate;
f = [4:2:40];

alphaEvidenceAccumulation = 0.96; % smoothing parameter

% extract the true stopping points
stopId = 555;
%stopTimesInSR = sessions{1,4}.EVENT.POS(sessions{1,4}.EVENT.TYP == 555);
stopTimesInSR = s{1,4}.event.position(s{1,4}.event.name == 555);

OptFeat=10; 
indexes = indexPower(1:OptFeat);%look for the first f best features
newFeatures = TrainingData(:,indexes);%only the features selected
classifier_online = fitcdiscr(newFeatures,TLabels(1:8730),'DiscrimType',type);

testData = spatFilteredData{1,4}; % data for the pseudo online classification
i = 1;
for t = 512:32:size(testData,2)
    for idxChannels = 1:NumChannels
        bufferedData = testData(idxChannels,t-512+1:t);
        [pxx(:,idxChannels,i),~] = pwelch(bufferedData,psdWindow,psdNOverlap,f,512);
        %pxx(:,idxChannels,i) = log(pxx(:,idxChannels,i)); % take log for correct scaling
    end 
    i= i+1;
end

% reshape pxx to have same dimensionality as feature matrix
pxx = reshape(pxx,[19*16,size(pxx,3)]);
pxx = pxx';
% extract best features
pxx_bestFeat = pxx(:,indexes);
% normalize wrt newfeatures matrix
pxx_bestFeat = zscore(pxx_bestFeat);

% predict (with score output)
[~,rawScore_bufferedData] = predict(classifier_online,pxx_bestFeat);
% smooth score 
smoothScore_bufferedData = zeros(size(pxx_bestFeat,1),2);
smoothScore_bufferedData(1,:) = 0.5; % initialize first element at 0.5 for equal probs
for idx = 2:size(pxx_bestFeat,1)
    smoothScore_bufferedData(idx,:) = alphaEvidenceAccumulation*smoothScore_bufferedData(idx-1,:) ...
        + (1-alphaEvidenceAccumulation)*rawScore_bufferedData(idx,:);
end

x = 0.5*ones(size(pxx_bestFeat,1),1);
x_axisVector = 1/16:1/16:(size(smoothScore_bufferedData,1))/16;
stopTimesInSeconds = stopTimesInSR./16;

% general figure over whole session
figure(j)
hold on
plot(x_axisVector,smoothScore_bufferedData(:,1))
plot(x_axisVector,x,'k--');
for idxLine = 1:30
    plot([stopTimesInSeconds(idxLine) stopTimesInSeconds(idxLine)],[0 1],'r');
end
hold off
title(['smoothed pseude online classification for subject: ',testPerson,', using alpha = ',num2str(alphaEvidenceAccumulation)],'FontSize',24)
xlabel('time [s]','FontSize',22)
ylabel('decoder probability of motor termination','FontSize',22)
legend('motor termination','threshold 50%','stop-times')
j = j-1;
k = 50;
%% subsampling around stop times + figure averaging

windowAroundStop = -2:1/16:3;
threshShort = 0.5*ones(size(windowAroundStop,2),1);

for idxTrial = 1:30
    startIdx = stopTimesInSR(idxTrial)+windowAroundStop(1)*16;
    stopIdx = stopTimesInSR(idxTrial)++windowAroundStop(end)*16;
    smoothScoreAroundStop(:,idxTrial) = smoothScore_bufferedData([startIdx:stopIdx],1);
end

averageSmoothScoreAroundStop = mean(smoothScoreAroundStop,2);

figure(k)
for idxSP = 1:30
    subplot(5,6,idxSP)
    hold on;
    plot(windowAroundStop,smoothScoreAroundStop(:,idxSP))
    plot(windowAroundStop,threshShort,'k--');
    hold off;
    title (['Trial: ',num2str(idxSP)]);
    ylim([0 1])
    xlabel('time [s]')
    ylabel('probability MT')
    legend('motor termination')%,'threshold 50%','stop-times')
end
k = k-1;

idxGoodTrials = [1,3,8,25];
idxBadTrials = [2,15];

figure(1)
for idx = 1:4
    subplot(2,2,idx)
    hold on;
    plot(windowAroundStop,smoothScoreAroundStop(:,idxGoodTrials(idx)))
    plot(windowAroundStop,threshShort,'k--');
    hold off;
    title (['Trial: ',num2str(idxGoodTrials(idx))],'FontSize',24);
    ylim([0 1])
    xlabel('time [s]','FontSize',22)
    ylabel('probability MT','FontSize',22)
    legend('motor termination')%,'threshold 50%','stop-times')
end

figure(2)
for idx = 1:2
    subplot(2,1,idx)
    hold on;
    plot(windowAroundStop,smoothScoreAroundStop(:,idxBadTrials(idx)))
    plot(windowAroundStop,threshShort,'k--');
    hold off;
    title (['Trial: ',num2str(idxBadTrials(idx))],'FontSize',24);
    ylim([0 1])
    xlabel('time [s]','FontSize',22)
    ylabel('probability MT','FontSize',22)
    legend('motor termination')%,'threshold 50%','stop-times')
end