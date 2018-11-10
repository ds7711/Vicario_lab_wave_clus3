function fillret = create_wavemark(fhand, temp_spiketiming, temp_spikewave, iChan, channel_name, param)
% create a wavemark channel based on timing and spikes

spike_ndpts = param.ndpts;
dRate = param.dRate_wavemark;

sr_25k_ratio = CEDS64ChanDiv( fhand, iChan );
shift_in_tick = param.shift * sr_25k_ratio;

wmarkerbuffer(length(temp_spiketiming), 1) = CEDWaveMark();
for nnn = 1 : length(temp_spiketiming)
    tmp_tick = CEDS64SecsToTicks(fhand, temp_spiketiming(nnn));
    wmarkerbuffer(nnn).SetTime(tmp_tick + shift_in_tick);
    if param.normalize_wav_mag
        tmp_data = transpose(temp_spikewave(nnn, :)) / max(abs(temp_spikewave(nnn, :)));
        tmp_data = tmp_data * param.int16_factor;
    else
        tmp_data = transpose(temp_spikewave(nnn, :)) * param.int16_factor / 3.0;

    end
    wmarkerbuffer(nnn).SetData(tmp_data);
end
% create the wavemark channel
wmarkchan = CEDS64GetFreeChan( fhand );

createret = CEDS64SetExtMarkChan(fhand, wmarkchan, dRate, 6, spike_ndpts, 1, sr_25k_ratio);
if createret ~= 0, warning('wavemarker channel not created correctly'); 
    disp(wav_name);
end;

CEDS64ChanTitle( fhand, wmarkchan, channel_name);
if param.normalize_wav_mag
    tmp_comment = ['normalized waveform * ', int2str(param.int16_factor), ' to integer'];
else
    tmp_comment = ['unnormalized waveform * ', int2str(param.int16_factor / 3.0), ' to integer'];
end
CEDS64ChanComment( fhand, wmarkchan, tmp_comment);
% CEDS64ChanUnits( fhand, wmarkchan, 'KA' );
fillret = CEDS64WriteExtMarks(fhand, wmarkchan, wmarkerbuffer);
if fillret < 0, warning('wave-marker channel not filled correctly'); end;
end