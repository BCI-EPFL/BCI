function spectrogramData = preprocess_spectrogram(runProcessedData,params_spectrogram)

spectrogramData = runProcessedData;

nRun = numel(runProcessedData);
for iRun = 1:nRun
    thisRun = runProcessedData{iRun};
    % Compute fast spectrogram
    [psd, freqgrid] = proc_spectrogram(thisRun.data', params_spectrogram.wlength, params_spectrogram.wshift, params_spectrogram.pshift, thisRun.rate, params_spectrogram.mlength);
    psd = log(psd);
    
    % Selecting desired frequencies
    [freqs, idfreqs] = intersect(freqgrid,params_spectrogram.freq);
    psd = psd(:, idfreqs, :);
    
    % Event alignement and subsampling
    cevents     = thisRun.event;
    events.name = cevents.name;
    events.position = proc_pos2win(cevents.position, params_spectrogram.wshift*thisRun.rate, 'backward', params_spectrogram.mlength*thisRun.rate);
    
    % Save structure for each run
    spectrogramData{iRun}.event = events;
    spectrogramData{iRun}.freq = freqs;
    spectrogramData{iRun}.data = permute(psd,[3 2 1]);
    spectrogramData{iRun}.rate = 1/params_spectrogram.wshift;
end
end

