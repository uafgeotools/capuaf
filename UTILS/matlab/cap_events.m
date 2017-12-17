% Examples for generating event set for CAP inverison
% calls read_eq_AEICdb and read_mech_AEICfp
% Vipul Silwal
% 17-oct-2013

axmoos = [-154 -146 58 62.5];
axmfsz = [-151 -147.5 63.5 65.5];
axstn = [-158 -139 55 67];

% Events for CAP (call read_eq_AEICdb)
ax3 = [axmoos -10 200];
oran = [datenum(2007,8,15) datenum(2009,8,15)];
oran = datenum('2014/04/16 20:24:24.814'); % near denali earthquake
ax3 = [-149 -148 62.8863 63 -10 200]
Mwran = [3.5 10];
[otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
figure; scatter3(lon,lat,dep,6^2,'filled');
title('Events for MT inversion using CAP');
colorbar; axis tight;

% seismicity 
ax3 = [axstn -10 200];
oran = [datenum(1990,1,1) datenum(2015,1,1)];
Mwran = [3 10];
[otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
figure; scatter3(lon,lat,dep,6^2,'filled');
title('Events for MT inversion using CAP');
colorbar; axis tight;

% Events for CAP along with first motion focal mechanism (call
% read_mech_AEICfp)
oran = [datenum(2007,8,15) datenum(2009,8,15)];
ax3 = [axmoos -10 700];
Mwran = [3.5 10];
[otime,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AECfp(oran,ax3,Mwran);
%[~,isort] = sort(eid);
[~,isort] = sort(sdep,'ascend');
%[~,isort] = sort(Mw,'descend');
display_eq_list(isort,otime,slon,slat,sdep,Mw,eid);
%write_psmeca('scak_event.dat',otime,slat,slon,sdep,M,eid)

% subset example (MOOS)
oran = [datenum(2007,8,15) datenum(2009,8,15)];
%oran = [datenum(1900,8,15) datenum(2015,8,15)];
ax3 = [axmoos -10 700];
Mwran = [0 10];
[otime,slat,slon,sdep,M,M0,Mw,eid,depc] = read_mech_AEC(oran,ax3,Mwran);
[~,isort] = sort(otime,'ascend');
%[~,isort] = sort(Mw,'descend');
display_eq_list(isort,otime,slon,slat,sdep,Mw,eid,depc);

% Nenana events (AEC solutions)
oran = [datenum(1999,1,1) datenum(2100,1,1)];  % BEAAR to present
ax3 = [axmfsz -10 700];
Mwran = [0 10];
[otime,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AECfp(oran,ax3,Mwran);
%[~,isort] = sort(eid);
[~,isort] = sort(sdep,'ascend');
%[~,isort] = sort(Mw,'descend');
display_eq_list(isort,otime,slon,slat,sdep,Mw,eid);

% Display eid and focal mechanism 
[gamma,delta,M0,kappa,theta,sigma,K,N,S,mu,lam] = CMT2TT(M,0);
for ii=1:length(M)
    disp(sprintf('%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f',eid{ii},Mw(ii),kappa(ii),theta(ii),sigma(ii),delta(ii),gamma(ii),sdep(ii)));
end

%==========================================================================

% FOR COMPARING EVENTS IN DIFFERENT CATALOGS, SEE alaska_tomo_sources.m
% v 1302 will generate this list:
% list of events common to A and B (A values):
%   1   3 otime 2007-10-03 (2007276) 14:06:12 lon -151.29 lat  58.28 dep  40.00 km Mw 5.17 20071003140612 C200710031406A
%   2   5 otime 2007-11-28 (2007332) 23:57:03 lon -151.13 lat  61.91 dep  75.00 km Mw 4.82 20071128235703 C200711282357A
%   3   6 otime 2008-03-14 (2008074) 09:38:21 lon -152.64 lat  61.07 dep 165.00 km Mw 4.97 20080314093821 C200803140938A
%   4   7 otime 2008-03-27 (2008087) 23:07:45 lon -152.17 lat  59.01 dep  75.00 km Mw 5.26 20080327230745 C200803272307A
%   5  11 otime 2009-01-24 (2009024) 18:09:50 lon -152.89 lat  59.43 dep  95.00 km Mw 5.72 20090124180950 C200901241809A
%   6  17 otime 2009-04-30 (2009120) 04:54:57 lon -151.31 lat  58.99 dep  40.00 km Mw 4.88 20090430045457 C200904300454A
%   7  19 otime 2009-06-22 (2009173) 19:28:05 lon -150.70 lat  61.94 dep  80.00 km Mw 5.51 20090622192805 C200906221928A

% FOR SMALL EVENTS, RUN rs_moos.m
% display_eq_list.m: listing 7/7 events
%    1      1 otime 2009-06-10 08:14:52 lon -149.39 lat   61.38 dep  27.10 km M  2.58 20090610081452
%    2      2 otime 2009-07-03 19:23:32 lon -148.85 lat   61.13 dep   1.40 km M  2.66 20090703192332
%    3      3 otime 2009-08-14 18:31:23 lon -149.31 lat   61.41 dep  27.01 km M  2.84 20090814183123
%    4      4 otime 2008-05-21 08:43:31 lon -148.98 lat   60.48 dep  20.69 km M  2.65 20080521084331
%    5      5 otime 2008-05-08 01:37:31 lon -149.97 lat   60.77 dep  28.74 km M  2.63 20080508013731
%    6      6 otime 2008-06-10 21:20:58 lon -149.85 lat   60.11 dep  28.43 km M  2.59 20080610212058
%    7      7 otime 2008-06-13 02:33:28 lon -149.34 lat   61.99 dep   3.23 km M  2.52 20080613023328

% FOR MFSZ EVENTS, RUN nenana_seismicity.m
% Largest 10 events in EAST
% display_eq_list.m: listing 10/10 events
%    1      1 otime 2001-03-25 11:34:50 lon -149.25 lat   64.63 dep  22.20 km M  4.60 20010325113450
%    2      2 otime 2002-12-29 20:38:30 lon -148.61 lat   64.95 dep  17.63 km M  3.43 20021229203830
%    3      3 otime 2002-08-13 06:23:52 lon -149.12 lat   64.67 dep  17.12 km M  3.37 20020813062352
%    4      4 otime 2013-03-05 21:55:58 lon -148.73 lat   64.84 dep  17.27 km M  3.36 20130305215558
%    5      5 otime 2005-09-14 17:24:52 lon -149.42 lat   64.50 dep  20.96 km M  3.12 20050914172452
%    6      6 otime 2001-05-02 23:53:15 lon -149.24 lat   64.62 dep  22.05 km M  2.91 20010502235315
%    7      7 otime 2006-02-11 22:05:54 lon -148.87 lat   64.82 dep  17.53 km M  2.87 20060211220554
%    8      8 otime 2008-02-19 17:47:01 lon -148.77 lat   64.88 dep   7.80 km M  2.83 20080219174701
%    9      9 otime 2009-01-16 22:51:38 lon -148.64 lat   64.92 dep  18.19 km M  2.78 20090116225138
%   10     10 otime 2003-12-29 14:41:42 lon -148.86 lat   64.83 dep  16.28 km M  2.77 20031229144142
% Largest 10 events in WEST
% display_eq_list.m: listing 10/10 events
%    1      1 otime 2000-02-07 04:17:53 lon -149.21 lat   64.71 dep  10.00 km M  4.10 20000207041753
%    2      2 otime 2008-07-16 10:12:00 lon -149.53 lat   64.59 dep  30.52 km M  4.10 20080716101200
%    3      3 otime 2012-04-11 09:21:57 lon -148.95 lat   64.92 dep  19.34 km M  3.93 20120411092157
%    4      4 otime 2009-07-28 12:13:15 lon -149.49 lat   64.61 dep  22.68 km M  3.71 20090728121315
%    5      5 otime 1999-12-23 02:42:25 lon -149.05 lat   64.87 dep  27.09 km M  3.67 19991223024225
%    6      6 otime 2004-11-17 11:29:00 lon -149.10 lat   64.89 dep  18.81 km M  3.60 20041117112900
%    7      7 otime 1999-12-18 16:55:43 lon -149.04 lat   64.88 dep  31.88 km M  3.47 19991218165543
%    8      8 otime 2011-11-18 10:46:23 lon -148.94 lat   64.94 dep  11.12 km M  3.44 20111118104623
%    9      9 otime 2011-10-19 23:42:26 lon -149.34 lat   64.65 dep  19.63 km M  3.39 20111019234226
%   10     10 otime 2013-09-23 04:12:37 lon -149.00 lat   64.81 dep  15.28 km M  3.39 20130923041237
% Largest 10 events in NORTH
% display_eq_list.m: listing 10/10 events
%    1      1 otime 2013-07-12 07:59:17 lon -148.77 lat   65.09 dep  17.16 km M  3.40 20130712075917
%    2      2 otime 2008-07-20 08:09:37 lon -148.64 lat   65.12 dep  20.16 km M  3.14 20080720080937
%    3      3 otime 2008-11-09 20:09:39 lon -148.61 lat   65.10 dep  20.00 km M  3.04 20081109200939
%    4      4 otime 2007-08-25 02:28:19 lon -148.54 lat   65.11 dep  13.72 km M  2.89 20070825022819
%    5      5 otime 1999-10-25 05:19:03 lon -148.51 lat   65.03 dep   5.23 km M  2.56 19991025051903
%    6      6 otime 1999-04-06 14:27:50 lon -148.49 lat   65.22 dep  16.05 km M  2.50 19990406142750
%    7      7 otime 2010-05-21 08:38:10 lon -148.65 lat   65.09 dep  15.54 km M  2.29 20100521083810
%    8      8 otime 2009-03-23 19:25:10 lon -148.72 lat   65.09 dep   8.47 km M  2.19 20090323192510
%    9      9 otime 2009-10-03 16:53:39 lon -148.60 lat   65.05 dep  16.30 km M  2.17 20091003165339
%   10     10 otime 2009-10-03 16:51:19 lon -148.67 lat   65.08 dep  13.91 km M  2.14 20091003165119
% Largest 10 events in SOUTH
% display_eq_list.m: listing 10/10 events
%    1      1 otime 2000-11-29 10:35:47 lon -150.35 lat   63.90 dep  16.39 km M  5.80 20001129103547
%    2      2 otime 2000-12-06 18:40:26 lon -150.31 lat   63.89 dep  11.69 km M  5.00 20001206184026
%    3      3 otime 2001-06-30 09:41:42 lon -150.15 lat   64.04 dep  14.61 km M  4.40 20010630094142
%    4      4 otime 2002-04-06 15:07:09 lon -150.00 lat   64.15 dep  26.70 km M  3.80 20020406150709
%    5      5 otime 2000-12-02 16:25:05 lon -150.29 lat   63.98 dep  10.00 km M  3.60 20001202162505
%    6      6 otime 2000-12-13 08:20:56 lon -150.21 lat   63.93 dep  14.27 km M  3.53 20001213082056
%    7      7 otime 2000-12-06 19:33:41 lon -150.39 lat   63.87 dep   7.70 km M  3.50 20001206193341
%    8      8 otime 2004-04-05 09:20:35 lon -150.07 lat   64.16 dep  13.59 km M  3.43 20040405092035
%    9      9 otime 2000-12-02 04:25:07 lon -150.25 lat   63.93 dep  14.28 km M  3.40 20001202042507
%   10     10 otime 2000-11-29 15:19:19 lon -150.21 lat   63.92 dep  13.55 km M  3.30 20001129151919
% Largest 10 events in SOUTHWEST
% display_eq_list.m: listing 10/10 events
%    1      1 otime 2013-06-05 18:58:23 lon -149.68 lat   64.64 dep  13.89 km M  3.99 20130605185823
%    2      2 otime 2010-08-29 18:08:06 lon -149.61 lat   64.70 dep  16.69 km M  2.91 20100829180806
%    3      3 otime 2005-07-04 21:15:15 lon -149.63 lat   64.70 dep  19.09 km M  2.79 20050704211515
%    4      4 otime 2006-01-09 04:38:24 lon -149.85 lat   64.51 dep   6.62 km M  2.28 20060109043824
%    5      5 otime 2013-04-09 13:56:29 lon -149.93 lat   64.57 dep  19.70 km M  2.28 20130409135629
%    6      6 otime 2004-09-26 02:40:52 lon -149.96 lat   64.56 dep  17.43 km M  2.26 20040926024052
%    7      7 otime 2012-08-25 23:20:48 lon -149.58 lat   64.71 dep   8.97 km M  2.25 20120825232048
%    8      8 otime 2002-11-23 08:28:39 lon -149.68 lat   64.70 dep  18.28 km M  2.11 20021123082839
%    9      9 otime 2012-11-05 17:21:35 lon -149.61 lat   64.66 dep  16.44 km M  2.11 20121105172135
%   10     10 otime 2013-03-01 04:00:27 lon -149.67 lat   64.64 dep  21.90 km M  2.10 20130301040027
  
  
