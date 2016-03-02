% clear all variables to begin with
clear; clc; clf;
opengl hardware; % use hardware openGL for speed
drawnow;

% load the audio files
fprintf('R E A D Y...\n');
[ready.dat, Fready] = wavread('ready.wav');
[relax.dat, Frelax] = wavread('close_your_eyes.wav');
[conc.dat, Fconc] = wavread('concentrate.wav');

Fs = 44100;
% need a high enough value so that alpha power below baseline can be played
Fc = 500; 
Fi = 500;

fprintf('R E A D Y...\n');


% play sound to indicate the start of the demo
readyobj=audioplayer(ready.dat, Fready);
playblocking(readyobj);

% number of passes to average for baseline
blPass = 15;

% number of passes for stimulus period
stPass = 75;

% total number of passes for the demo
totPass = 90;

% FieldTrip buffer source
filename = 'buffer://localhost:1972';

% read the header for the first time to determine number of channels and sampling rate
hdr = ft_read_header(filename, 'cache', true);

count      = 0;
prevSample = 0;
blocksize  = hdr.Fs;
chanindx   = 1:hdr.nChans; % all channels

params.tapers = [1 1]; % tapers
params.pad = -1; % no padding
params.Fs = hdr.Fs; % sampling frequency
params.trialave = 1; % average over trials
params.fpass = [0 50];

alphaUpperLimit = repmat(12, 1, totPass+1);
alphaLowerLimit = repmat(8, 1, totPass+1);
timeRange = linspace(0,totPass,totPass+1);
freqRange = linspace(0,50,51);

% create the plots
colorLimsRawTF = [-3 3];
colorLimsChangeTF = [-15 15];

hRawSpectrum = subplot('Position',[0.075 0.55 0.525 0.4]);
hold on;
plot(hRawSpectrum, timeRange, alphaUpperLimit, 'k--');
plot(hRawSpectrum, timeRange, alphaLowerLimit, 'k--');
xlabel(hRawSpectrum, 'Time (s)'); ylabel(hRawSpectrum, 'Frequency');
xlim(hRawSpectrum, [1 totPass]);
ylim(hRawSpectrum, [0 50]);
title(hRawSpectrum, 'Log (Raw power spectrum)');
caxis(hRawSpectrum, colorLimsRawTF);
colorbar('peer',hRawSpectrum);

hChangeSpectrum = subplot('Position',[0.075 0.05 0.525 0.4]);
hold on;
plot(hChangeSpectrum, timeRange, alphaUpperLimit, 'k--');
plot(hChangeSpectrum, timeRange, alphaLowerLimit, 'k--');
xlabel(hChangeSpectrum, 'Time (s)'); ylabel(hChangeSpectrum, 'Frequency');
xlim(hChangeSpectrum, [1 totPass]);
ylim(hChangeSpectrum, [0 50]);
title(hChangeSpectrum, 'Change in power spectrum (dB)');
caxis(hChangeSpectrum, colorLimsChangeTF);
colorbar('peer',hChangeSpectrum);

halpha = subplot('Position',[0.675 0.08 0.250 0.6]);
hold on;
xlabel(halpha, 'time(s)'); ylabel(halpha, 'increment factor');
title(halpha, '% of modulation');
xlim(halpha, [blPass stPass]);
ylim(halpha,[-2 2]);


blCount = blPass;

while (count < totPass)

  % determine number of samples available in buffer
  hdr = ft_read_header(filename, 'cache', true);

  % see whether new samples are available
  newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

  if newsamples >= blocksize

    % determine the samples to process
    begsample  = prevSample+1;
    endsample  = prevSample+blocksize ;

    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    blCount     = blCount - 1;
    %fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

    % read data segment from buffer
    nextdat = ft_read_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx);

    [power(count,:),freq] = mtspectrumc(nextdat',params);

    % plot the raw TF spectrum continously
    if (count > 1)
        subplot(hRawSpectrum);
        pcolor(1:size(power,1), freq, double(power'));
        shading interp;
        plot(hRawSpectrum, timeRange, alphaUpperLimit, 'k--');
        plot(hRawSpectrum, timeRange, alphaLowerLimit, 'k--');
        xlabel(hRawSpectrum, 'Time (s)'); ylabel(hRawSpectrum, 'Frequency');
        xlim(hRawSpectrum, [1 totPass]);
        ylim(hRawSpectrum, [0 50]);
        title(hRawSpectrum, 'Log (Raw power spectrum)');
        caxis(hRawSpectrum, colorLimsRawTF);
        colorbar('peer',hRawSpectrum);
        

        drawnow;
    end
    
    % calculate power in various bands (baseline/test)
    if (blCount >= 0)
        
      if (blCount == 0)
            mLogBL = mean(log10(power(5:15,:)));
            
            for kk=1:count
                dPower(kk,:) = log10(power(kk,:)) - mLogBL;
            end



            % cue user's attention
            fprintf('R E L A X...\n');
             % Relax now and close your eyes
            relaxobj = audioplayer(relax.dat,Frelax);
            playblocking(relaxobj);
            
            % advance last remembered position by 1 sec to account for
            % message playback above
            prevSample = prevSample + 2*blocksize; 
        end
        

        continue;
    end

    if (count == stPass)
        fprintf('open your eyes and relax...\n');
        concobj = audioplayer(conc.dat,Fconc);
        playblocking(concobj);
                    
        % advance last remembered position by 2 sec to account for
        % message playback above
        prevSample = prevSample + 2*blocksize; 
    end
    
    nextPower = log10(power(count,:)) - mLogBL;
    dPower = [dPower; nextPower];

    % get the frequency
    smoothKernel = repmat(1/10,1,10);
    epochsToAvg = length(smoothKernel);
    incrFact(count) = mean(dPower(end-epochsToAvg+1:end,8:15)'*smoothKernel');
    incrFactRaw(count) = mean(dPower(end,8:15));
    fprintf(['count ' num2str(count) ': incrRaw = ' num2str(incrFactRaw(count)) ', incrFact = ' num2str(incrFact(count)) '\n']);
    stFreq = round(Fc + incrFact(count) * Fi);
    
    % play the sound after baseline epoch is over
    soundTone = sine_tone(Fs,1,stFreq);
    sound(soundTone,Fs);
    
    subplot(hChangeSpectrum);
    pcolor(1:size(dPower,1), freq, 10*double(dPower(1:size(dPower,1),:)'));
    shading interp;
    plot(hChangeSpectrum, timeRange, alphaUpperLimit, 'k--');
    plot(hChangeSpectrum, timeRange, alphaLowerLimit, 'k--');
    xlabel(hChangeSpectrum, 'Time (s)'); ylabel(hChangeSpectrum, 'Frequency');
    xlim(hChangeSpectrum, [1 totPass]);
    ylim(hChangeSpectrum, [0 50]);
    caxis(hChangeSpectrum, colorLimsChangeTF);
    
    title(hChangeSpectrum, 'Change in power spectrum');
   
     subplot(halpha);
     plot(incrFactRaw,'b');
     hold on;
     plot(incrFact,'r');
     xlabel(halpha, 'time(s)'); ylabel(halpha, 'increment factor');
     title(halpha, '% of modulation');

     xlim(halpha, [0 totPass]);

    drawnow;

  end % if 
end % while true


% reset figure
fprintf('End of the demo');


analysisRange = 26:75;

blPowerArray = mean(power(5:15,8:15),2);
stPowerArray = mean(power(analysisRange,8:15),2);

changeArray = stPowerArray/mean(blPowerArray) - 1;

quot = 100*mean(changeArray);
fluct = std(stPowerArray)/mean(stPowerArray);
% show a message box
msgbox(['Your relaxation quotient is ' num2str(quot) ' % and your sustenance quotient is ' num2str(fluct)], 'EEG Demo', 'help');
save 'alphacontrol.mat';


figure
plot(1:count,incrFact,'r');
hold on;
clear incrFact;
load('alphaAwareness.mat','incrFact');
plot(1:count,incrFact,'b');
legend('alpha suppression','alpha enhancement');
