% Show the spectrogram, using "bach_fugue_short.mat"
load("bach_fugue_short.mat");
% short ver doesn't have the third voice
fs = 11025;
bpm = 120;
beats_per_second = bpm/60;
seconds_per_beat = 1/beats_per_second;
seconds_per_pulse = seconds_per_beat/4;
X = 1; % change X with E(t) 
% 3 voices in the song use the loop to put the data into xx
fugue = cell(1,2);
for i=1:2
    % data type problem zeros(1,N) need N to be an positive integral!
    xx = zeros(1,ceil(sum(theVoices(i).durations*seconds_per_pulse)+1)*fs);
    n1=theVoices(i).startPulses(1);
    % use loop to put the E(t) for each note!
    for j = 1:length(theVoices(i).noteNumbers)
        keynum = theVoices(i).noteNumbers(j);
        dur = theVoices(i).durations(j)*seconds_per_pulse;
        note = key2note(X,keynum,dur);
        % ADSR envelope array length setting
        Envelope_length = length(note);
        A = round(0.15*Envelope_length);
        D = round(0.1*Envelope_length);
        R = round(0.2*Envelope_length);
        S = Envelope_length - A - R - D;
        % linspace to plot the figure 5 
        Attack = linspace(0,1,A);
        Decay = linspace(1,0.8,D);      
        Sustain   = linspace(0.8,0.7,S); 
        Release = linspace(0.7,0,R);
        % Add the arrays into the envelope
        ADSR_Envelope = [Attack, Decay, Sustain, Release];
        ADSR_Envelope = ADSR_Envelope(1:length(note));
        % Adjust the X=1 into E(t)
        note = ADSR_Envelope.*note;
        % Add the notes like playscale.m
        n2 = n1 + length(note) - 1;    
        xx(n1:n2) = xx(n1:n2) + note;
        n1 = theVoices(i).startPulses(j)+n2;
    end
    % the {} in array is used to save arrays in an array!
    fugue{i} = xx;
end
% Add the 3 fugue voices into a single song by sum()
bach_fugue_short = sum([fugue{:}],1);
spectrogram(bach_fugue_short);

% Make the song, using data, "bach_fugue.mat". 
load("bach_fugue.mat");
fs = 11025;
bpm = 120;
beats_per_second = bpm/60;
seconds_per_beat = 1/beats_per_second;
seconds_per_pulse = seconds_per_beat/4;
X = 1; % change X with E(t) 
% 3 voices in the song use the loop to put the data into xx
fugue = cell(1,3);
for i=1:3
    % data type problem zeros(1,N) need N to be an positive integral!
    xx = zeros(1,ceil(sum(theVoices(i).durations*seconds_per_pulse)+1)*fs);
    n1=theVoices(i).startPulses(1);
    % use loop to put the E(t) for each note!
    for j = 1:length(theVoices(i).noteNumbers)
        keynum = theVoices(i).noteNumbers(j);
        dur = theVoices(i).durations(j)*seconds_per_pulse;
        note = key2note(X,keynum,dur);
        % ADSR envelope array length setting
        Envelope_length = length(note);
        A = round(0.15*Envelope_length);
        D = round(0.1*Envelope_length);
        R = round(0.2*Envelope_length);
        S = Envelope_length - A - R - D;
        % linspace to plot the figure 5 
        Attack = linspace(0,1,A);
        Decay = linspace(1,0.8,D);      
        Sustain   = linspace(0.8,0.7,S); 
        Release = linspace(0.7,0,R);
        % Add the arrays into the envelope
        ADSR_Envelope = [Attack, Decay, Sustain, Release];
        ADSR_Envelope = ADSR_Envelope(1:length(note));
        % Adjust the X=1 into E(t)
        note = ADSR_Envelope.*note;
        % Add the notes like playscale.m
        n2 = n1 + length(note) - 1;    
        xx(n1:n2) = xx(n1:n2) + note;
        n1 = theVoices(i).startPulses(j)+n2;
    end
    % the {} in array is used to save arrays in an array!
    fugue{i} = xx;
end
% Add the 3 fugue voices into a single song by sum(A,1)
bach_fugue = sum([fugue{:}],1);
soundsc(bach_fugue)