function [wavdata] = readWave(fhand, iChan)
maxTimeTicks = CEDS64ChanMaxTime( fhand, 1 )+1;
sr = CEDS64TimeBase(fhand);
margin = int64(0.1 / sr);
readUntil = maxTimeTicks - margin;

clearvars margin;

[ lastTick, fVals, ~ ] = CEDS64ReadWaveF( fhand, iChan, maxTimeTicks, 0, maxTimeTicks);

% initialize a vector to store the data
wavdata = zeros(maxTimeTicks, 1, 'single');
wavdata(1:lastTick) = fVals;

% read until the last time ticks
while lastTick < readUntil
    % if not at the end of channel yet, read agian
    
    % starting tick is lastick
    [ fRead, fVals, firstTick ] = CEDS64ReadWaveF( fhand, iChan, maxTimeTicks, lastTick+1, maxTimeTicks);
    
    % if some values are missing, add zeros in between
    wavdata(firstTick+1 : firstTick + fRead) = fVals; 
    
    % get next tick to read from
    lastTick = firstTick + fRead;
    
    % if data chunk is less than 1 sec and 90% of data are read, terminate
    if fRead * sr < 1 && maxTimeTicks * 0.9 < lastTick 
        wavdata = wavdata(1:lastTick);
        break;
    end
    
end

end