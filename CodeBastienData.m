% no division on train and test
% took all the data and to the 
addpath(genpath('DataBastien'));
load('ak4/training4.mat');
% load('ak6/training5');
% load('ak6/training4');




%do the rank feat to all the data
for j=1:size(training.features.spectral,3)
training_data(:,j)=reshape(training.features.spectral(:,:,j), [34*120,1]);

end

training_labels=reshape(training.label, [34*120,1]);

[ind, power_feat] = rankfeat(training_data, training_labels,  'fisher');


for i=1:size(power_feat,2)
averageFisher(1,i)=power_feat(ind==i);
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisherReshape=(reshape(averageFisher,[19 16]))';

%plot of Fisher's scores.
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisherReshape)
axis tight
title('Fisher scores - MI initiation vs MI termination');
xlabel('Frequency [Hz]');
ylabel('channels');
set(gca,'yTick',1:16,'YTickLabel', {chanlocs16.labels});
h=colorbar;


%% ROC CURVE


training_data1=training_data(1:(0.75*size(training_data,1)),:);
training_labels1=training_labels(1:0.75*size(training_labels,1));

test_data1=training_data((0.75*size(training_data,1)):end,:);
test_labels1=training_labels((0.75*size(training_labels,1)):end);

classifierBastien=fitcdiscr(training_data1,training_labels1,'discrimtype', 'linear');
[yhat,PosteriorProb,~]=predict(classifierBastien,test_data1);

[X,Y] = perfcurve(test_labels1,PosteriorProb(:,2),1);
plot(X,Y)
