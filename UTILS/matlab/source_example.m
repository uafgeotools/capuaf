% write sac file
% The outputted sac file can be used as an input to CAP
% 
% Copy this file inside the directory containing RTZ sac waveforms
%
% Example command: 
% cap.pl -H0.02 -P1/20/60 -p1 -S3/10/0 -T15/120 -D1/1/0.5 -C0.2/0.5/0.025/0.06 -W1 -Mscak_41 -m4.5 -I100000 -Zweight111_subset.dat -Y1 -R0/0 -Lsin_source.sac 20090407201255351
%
% See example: cap/EXAMPLES/run_example_alaska.sh (example 10)

nsamples = 100;
sstart = 0;
send = pi; 
DELTA = 0.02; % sac header
sacfilename = './sin_source.sac';
srcdelay = 0;
amp = 0.02; % this is not clear (magnitude dependent?)

% this is not the discritization in time
% time vector is sampled in mksac.m using nsamples and DELTA
samp = linspace(sstart,send,nsamples);

% create source data vectore
src = amp*sin(samp);

% optional (for ploting)
tvec = 0:DELTA:((nsamples-1)*DELTA);
plot(tvec,src)
xlabel('time (s)');

% save sac
mksac(sacfilename,src,srcdelay,'DELTA',DELTA);