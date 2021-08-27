/**************************************************************************************************
programmer    :Hiroya Hashimioto
file name     :Simulation for 1-sample RMST.sas
date created  :2020/06/24
     modified :2021/08/23
description   :Simulation experiments to evaluate interval estimation using variable transformation of RMST
               (generating survival time data according to Weibull distribution).
**************************************************************************************************/

title ;
footnote ;

proc datasets library=work kill nolist ;
quit ;

dm 'output; clear; log; clear;';
options nonumber nonotes nomprint ls=100 ps=9999 formdlim="-" center ;



%let execpath = " ";
%let path = " ";
%macro setexecpath;
  %let execpath = %sysfunc(getoption(sysin));
  %if %length(&execpath) = 0 %then
    %let execpath = %sysget(sas_execfilepath);
  data _null_;
    do i = length("&execpath") to 1 by -1;
      if substr("&execpath", i, 1) = "\" then do;
        call symput("path", substr("&execpath", 1, i));
        stop;
      end;
    end;
  run;
%mend setexecpath;
%setexecpath;

libname output "&path";


proc format;
  value result 1="above error" 2="coverage" 3="below error";


/***************************************************************
 Macro section
***************************************************************/



%macro rmst_sim(SIM_No, /* simulation number*/
                 INIT,   /* seed */
                 SIZE,   /* sample size*/
                 SHAPE,  /* shape parameter of Weible distribution */
                 SCALE,  /* scale parameter of Weible distribution */
                 TRUE,   /* ture value of RMST */
                );


data tmp1;
  call streaminit(&INIT);
  do SIM = 1 to &SIM_N;
    do I = 1 to &SIZE;
      T = rand("WEIBULL",&SHAPE,&SCALE);
      if T>1 then STATUS=0;
      else STATUS=1;
      output;
    end;
  end;
run;

ods listing close;
proc lifetest data=tmp1 timelim=1;
  by SIM;
  time T*STATUS(0);
  ods output Means=tmpa1 
             ProductLimitEstimates=tmpb1(where=(SURVIVAL ne . and T <= 1));
run;
ods output close;

proc sort data=tmpb1 out=tmpb2(keep=SIM FAILED);
  by SIM descending FAILED;
run;

proc sort data=tmpb2 out=tmpb3 nodupkey;
  by SIM;
run;

data tmp2;
  merge tmpa1 tmpb3;
  by SIM;
run;

%let c= 1-&ALPHA./2;

data tmp3;
  set tmp2;

  SE_KM=STDERR*sqrt((FAILED-1)/FAILED);

  NT_Gw_L= Mean - SE_KM*quantile('NORMAL',&c.);
  NT_Gw_U= Mean + SE_KM*quantile('NORMAL',&c.);
  CLL_Gw_L=exp(-exp( log(-log(Mean)) - quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*SE_KM ));
  CLL_Gw_U=exp(-exp( log(-log(Mean)) + quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*SE_KM ));
  ASR_Gw_L=(sin( arsin(sqrt(Mean)) - quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*SE_KM ))**2;
  ASR_Gw_U=(sin( arsin(sqrt(Mean)) + quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*SE_KM ))**2;
  Logit_Gw_L =-1/(1+exp( log(Mean/(1-Mean)) - quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*SE_KM )) +1;
  Logit_Gw_U =-1/(1+exp( log(Mean/(1-Mean)) + quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*SE_KM )) +1;

  NT_LM_L= Mean - StdErr*quantile('NORMAL',&c.);
  NT_KM_U= Mean + StdErr*quantile('NORMAL',&c.);
  CLL_KM_L=exp(-exp( log(-log(Mean)) - quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*StdErr ));
  CLL_KM_U=exp(-exp( log(-log(Mean)) + quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*StdErr ));
  ASR_KM_L=(sin( arsin(sqrt(Mean)) - quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*StdErr ))**2;
  ASR_KM_U=(sin( arsin(sqrt(Mean)) + quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*StdErr ))**2;
  Logit_KM_L =-1/(1+exp( log(Mean/(1-Mean)) - quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*StdErr )) +1;
  Logit_KM_U =-1/(1+exp( log(Mean/(1-Mean)) + quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*StdErr )) +1;

  if &TRUE < NT_Gw_L then NT_Gw_R=1;
    else if &TRUE < NT_Gw_U then NT_Gw_R=2;
    else NT_Gw_R=3;
  if &TRUE < CLL_Gw_L then CLL_Gw_R=1;
    else if &TRUE < CLL_Gw_U then CLL_Gw_R=2;
    else CLL_Gw_R=3;
  if &TRUE < ASR_Gw_L then ASR_Gw_R=1;
    else if &TRUE < ASR_Gw_U then ASR_Gw_R=2;
    else ASR_Gw_R=3;
  if &TRUE < Logit_Gw_L then Logit_Gw_R=1;
    else if &TRUE < Logit_Gw_U then Logit_Gw_R=2;
    else Logit_Gw_R=3;

  if &TRUE < NT_KM_L then NT_KM_R=1;
    else if &TRUE < NT_KM_U then NT_KM_R=2;
    else NT_KM_R=3;
  if &TRUE < CLL_KM_L then CLL_KM_R=1;
    else if &TRUE < CLL_KM_U then CLL_KM_R=2;
    else CLL_KM_R=3;
  if &TRUE < ASR_KM_L then ASR_KM_R=1;
    else if &TRUE < ASR_KM_U then ASR_KM_R=2;
    else ASR_KM_R=3;
  if &TRUE < Logit_KM_L then Logit_KM_R=1;
    else if &TRUE < Logit_KM_U then Logit_KM_R=2;
    else Logit_KM_R=3;

run;

ods listing close;
proc tabulate data=tmp3;
  class NT_Gw_R CLL_Gw_R ASR_Gw_R Logit_Gw_R NT_KM_R CLL_KM_R ASR_KM_R Logit_KM_R;
  var STRATUM;
  table STRATUM*n,(NT_Gw_R CLL_Gw_R ASR_Gw_R Logit_Gw_R NT_KM_R CLL_KM_R ASR_KM_R Logit_KM_R);
  ods output table=tmp4;
run;
ods output close;
ods listing;

data output.out&SIM_NO._&SIZE(keep=SIM_NO SIM N METHOD RESULT PROP);
  set tmp4;
  length METHOD $ 10.;
  SIM_NO=&SIM_NO;
  SIM=&SIM_N;
  N=&SIZE;
  if NT_Gw_R ne . then do;
    METHOD="NT_Gw";
    RESULT=NT_Gw_R; end;
  if CLL_Gw_R ne . then do;
    METHOD="CLL_Gw";
    RESULT=CLL_Gw_R; end;
  if ASR_Gw_R ne . then do;
    METHOD="ASR_Gw";
    RESULT=ASR_Gw_R; end;
  if Logit_Gw_R ne . then do;
    METHOD="Logit_Gw";
    RESULT=Logit_Gw_R; end;

  if NT_KM_R ne . then do;
    METHOD="NT_KM";
    RESULT=NT_KM_R; end;
  if CLL_KM_R ne . then do;
    METHOD="CLL_KM";
    RESULT=CLL_KM_R; end;
  if ASR_KM_R ne . then do;
    METHOD="ASR_KM";
    RESULT=ASR_KM_R; end;
  if Logit_KM_R ne . then do;
    METHOD="Logit_KM";
    RESULT=Logit_KM_R; end;
  PROP=STRATUM_N/&SIM_N;
  format RESULT result.;
run;

proc datasets lib=work;
  delete tmp:;
run;quit;

%mend;




/***************************************************************
 The part that executes the macro
 It is divided into 10000 times and executed because of the lack of memory.
***************************************************************/
%let SIM_N=100000; /* number of simulations */
%let ALPHA=0.05;   /* Type I error rate (2-sided) */

/*
  n=25
 */
%rmst_sim(11, 1501, 25, 0.7, 0.5067, 0.4262);
%rmst_sim(12, 1502, 25, 0.7, 0.5067, 0.4262);
%rmst_sim(13, 1503, 25, 0.7, 0.5067, 0.4262);
%rmst_sim(14, 1504, 25, 0.7, 0.5067, 0.4262);
%rmst_sim(15, 1505, 25, 0.7, 0.5067, 0.4262);

%rmst_sim(21, 2501, 25, 0.7, 1.6881, 0.6767);
%rmst_sim(22, 2502, 25, 0.7, 1.6881, 0.6767);
%rmst_sim(23, 2503, 25, 0.7, 1.6881, 0.6767);
%rmst_sim(24, 2504, 25, 0.7, 1.6881, 0.6767);
%rmst_sim(25, 2505, 25, 0.7, 1.6881, 0.6767);

%rmst_sim(31, 3501, 25, 0.7, 8.5230, 0.8785);
%rmst_sim(32, 3502, 25, 0.7, 8.5230, 0.8785);
%rmst_sim(33, 3503, 25, 0.7, 8.5230, 0.8785);
%rmst_sim(34, 3504, 25, 0.7, 8.5230, 0.8785);
%rmst_sim(35, 3505, 25, 0.7, 8.5230, 0.8785);

%rmst_sim(41, 4501, 25, 1.0, 0.6213, 0.4971);
%rmst_sim(42, 4502, 25, 1.0, 0.6213, 0.4971);
%rmst_sim(43, 4503, 25, 1.0, 0.6213, 0.4971);
%rmst_sim(44, 4504, 25, 1.0, 0.6213, 0.4971);
%rmst_sim(45, 4505, 25, 1.0, 0.6213, 0.4971);

%rmst_sim(51, 5501, 25, 1.0, 1.4427, 0.7213);
%rmst_sim(52, 5502, 25, 1.0, 1.4427, 0.7213);
%rmst_sim(53, 5503, 25, 1.0, 1.4427, 0.7213);
%rmst_sim(54, 5504, 25, 1.0, 1.4427, 0.7213);
%rmst_sim(55, 5505, 25, 1.0, 1.4427, 0.7213);

%rmst_sim(61, 6501, 25, 1.0, 4.4814, 0.8963);
%rmst_sim(62, 6502, 25, 1.0, 4.4814, 0.8963);
%rmst_sim(63, 6503, 25, 1.0, 4.4814, 0.8963);
%rmst_sim(64, 6504, 25, 1.0, 4.4814, 0.8963);
%rmst_sim(65, 6505, 25, 1.0, 4.4814, 0.8963);

%rmst_sim(71, 7501, 25, 1.5, 0.7281, 0.5850);
%rmst_sim(72, 7502, 25, 1.5, 0.7281, 0.5850);
%rmst_sim(73, 7503, 25, 1.5, 0.7281, 0.5850);
%rmst_sim(74, 7504, 25, 1.5, 0.7281, 0.5850);
%rmst_sim(75, 7505, 25, 1.5, 0.7281, 0.5850);

%rmst_sim(81, 8501, 25, 1.5, 1.2768, 0.7739);
%rmst_sim(82, 8502, 25, 1.5, 1.2768, 0.7739);
%rmst_sim(83, 8503, 25, 1.5, 1.2768, 0.7739);
%rmst_sim(84, 8504, 25, 1.5, 1.2768, 0.7739);
%rmst_sim(85, 8505, 25, 1.5, 1.2768, 0.7739);

%rmst_sim(91, 9001, 25, 1.5, 2.7182, 0.9166);
%rmst_sim(92, 9002, 25, 1.5, 2.7182, 0.9166);
%rmst_sim(93, 9003, 25, 1.5, 2.7182, 0.9166);
%rmst_sim(94, 9004, 25, 1.5, 2.7182, 0.9166);
%rmst_sim(95, 9005, 25, 1.5, 2.7182, 0.9166);


/*
  n=50
 */
%rmst_sim(11, 1001, 50, 0.7, 0.5067, 0.4262);
%rmst_sim(12, 1002, 50, 0.7, 0.5067, 0.4262);
%rmst_sim(13, 1003, 50, 0.7, 0.5067, 0.4262);
%rmst_sim(14, 1004, 50, 0.7, 0.5067, 0.4262);
%rmst_sim(15, 1005, 50, 0.7, 0.5067, 0.4262);

%rmst_sim(21, 2001, 50, 0.7, 1.6881, 0.6767);
%rmst_sim(22, 2002, 50, 0.7, 1.6881, 0.6767);
%rmst_sim(23, 2003, 50, 0.7, 1.6881, 0.6767);
%rmst_sim(24, 2004, 50, 0.7, 1.6881, 0.6767);
%rmst_sim(25, 2005, 50, 0.7, 1.6881, 0.6767);

%rmst_sim(31, 3001, 50, 0.7, 8.5230, 0.8785);
%rmst_sim(32, 3002, 50, 0.7, 8.5230, 0.8785);
%rmst_sim(33, 3003, 50, 0.7, 8.5230, 0.8785);
%rmst_sim(34, 3004, 50, 0.7, 8.5230, 0.8785);
%rmst_sim(35, 3005, 50, 0.7, 8.5230, 0.8785);

%rmst_sim(41, 4001, 50, 1.0, 0.6213, 0.4971);
%rmst_sim(42, 4002, 50, 1.0, 0.6213, 0.4971);
%rmst_sim(43, 4003, 50, 1.0, 0.6213, 0.4971);
%rmst_sim(44, 4004, 50, 1.0, 0.6213, 0.4971);
%rmst_sim(45, 4005, 50, 1.0, 0.6213, 0.4971);

%rmst_sim(51, 5001, 50, 1.0, 1.4427, 0.7213);
%rmst_sim(52, 5002, 50, 1.0, 1.4427, 0.7213);
%rmst_sim(53, 5003, 50, 1.0, 1.4427, 0.7213);
%rmst_sim(54, 5004, 50, 1.0, 1.4427, 0.7213);
%rmst_sim(55, 5005, 50, 1.0, 1.4427, 0.7213);

%rmst_sim(61, 6001, 50, 1.0, 4.4814, 0.8963);
%rmst_sim(62, 6002, 50, 1.0, 4.4814, 0.8963);
%rmst_sim(63, 6003, 50, 1.0, 4.4814, 0.8963);
%rmst_sim(64, 6004, 50, 1.0, 4.4814, 0.8963);
%rmst_sim(65, 6005, 50, 1.0, 4.4814, 0.8963);

%rmst_sim(71, 7001, 50, 1.5, 0.7281, 0.5850);
%rmst_sim(72, 7002, 50, 1.5, 0.7281, 0.5850);
%rmst_sim(73, 7003, 50, 1.5, 0.7281, 0.5850);
%rmst_sim(74, 7004, 50, 1.5, 0.7281, 0.5850);
%rmst_sim(75, 7005, 50, 1.5, 0.7281, 0.5850);

%rmst_sim(81, 8001, 50, 1.5, 1.2768, 0.7739);
%rmst_sim(82, 8002, 50, 1.5, 1.2768, 0.7739);
%rmst_sim(83, 8003, 50, 1.5, 1.2768, 0.7739);
%rmst_sim(84, 8004, 50, 1.5, 1.2768, 0.7739);
%rmst_sim(85, 8005, 50, 1.5, 1.2768, 0.7739);

%rmst_sim(91, 9001, 50, 1.5, 2.7182, 0.9166);
%rmst_sim(92, 9002, 50, 1.5, 2.7182, 0.9166);
%rmst_sim(93, 9003, 50, 1.5, 2.7182, 0.9166);
%rmst_sim(94, 9004, 50, 1.5, 2.7182, 0.9166);
%rmst_sim(95, 9005, 50, 1.5, 2.7182, 0.9166);


/*
  n=100
 */
%rmst_sim(11, 1101, 100, 0.7, 0.5067, 0.4262);
%rmst_sim(12, 1102, 100, 0.7, 0.5067, 0.4262);
%rmst_sim(13, 1103, 100, 0.7, 0.5067, 0.4262);
%rmst_sim(14, 1104, 100, 0.7, 0.5067, 0.4262);
%rmst_sim(15, 1105, 100, 0.7, 0.5067, 0.4262);

%rmst_sim(21, 2101, 100, 0.7, 1.6881, 0.6767);
%rmst_sim(22, 2102, 100, 0.7, 1.6881, 0.6767);
%rmst_sim(23, 2103, 100, 0.7, 1.6881, 0.6767);
%rmst_sim(24, 2104, 100, 0.7, 1.6881, 0.6767);
%rmst_sim(25, 2105, 100, 0.7, 1.6881, 0.6767);

%rmst_sim(31, 3101, 100, 0.7, 8.5230, 0.8785);
%rmst_sim(32, 3102, 100, 0.7, 8.5230, 0.8785);
%rmst_sim(33, 3103, 100, 0.7, 8.5230, 0.8785);
%rmst_sim(34, 3104, 100, 0.7, 8.5230, 0.8785);
%rmst_sim(35, 3105, 100, 0.7, 8.5230, 0.8785);

%rmst_sim(41, 4101, 100, 1.0, 0.6213, 0.4971);
%rmst_sim(42, 4102, 100, 1.0, 0.6213, 0.4971);
%rmst_sim(43, 4103, 100, 1.0, 0.6213, 0.4971);
%rmst_sim(44, 4104, 100, 1.0, 0.6213, 0.4971);
%rmst_sim(45, 4105, 100, 1.0, 0.6213, 0.4971);

%rmst_sim(51, 5101, 100, 1.0, 1.4427, 0.7213);
%rmst_sim(52, 5102, 100, 1.0, 1.4427, 0.7213);
%rmst_sim(53, 5103, 100, 1.0, 1.4427, 0.7213);
%rmst_sim(54, 5104, 100, 1.0, 1.4427, 0.7213);
%rmst_sim(55, 5105, 100, 1.0, 1.4427, 0.7213);

%rmst_sim(61, 6101, 100, 1.0, 4.4814, 0.8963);
%rmst_sim(62, 6102, 100, 1.0, 4.4814, 0.8963);
%rmst_sim(63, 6103, 100, 1.0, 4.4814, 0.8963);
%rmst_sim(64, 6104, 100, 1.0, 4.4814, 0.8963);
%rmst_sim(65, 6105, 100, 1.0, 4.4814, 0.8963);

%rmst_sim(71, 7101, 100, 1.5, 0.7281, 0.5850);
%rmst_sim(72, 7102, 100, 1.5, 0.7281, 0.5850);
%rmst_sim(73, 7103, 100, 1.5, 0.7281, 0.5850);
%rmst_sim(74, 7104, 100, 1.5, 0.7281, 0.5850);
%rmst_sim(75, 7105, 100, 1.5, 0.7281, 0.5850);

%rmst_sim(81, 8101, 100, 1.5, 1.2768, 0.7739);
%rmst_sim(82, 8102, 100, 1.5, 1.2768, 0.7739);
%rmst_sim(83, 8103, 100, 1.5, 1.2768, 0.7739);
%rmst_sim(84, 8104, 100, 1.5, 1.2768, 0.7739);
%rmst_sim(85, 8105, 100, 1.5, 1.2768, 0.7739);

%rmst_sim(91, 9101, 100, 1.5, 2.7182, 0.9166);
%rmst_sim(92, 9102, 100, 1.5, 2.7182, 0.9166);
%rmst_sim(93, 9103, 100, 1.5, 2.7182, 0.9166);
%rmst_sim(94, 9104, 100, 1.5, 2.7182, 0.9166);
%rmst_sim(95, 9105, 100, 1.5, 2.7182, 0.9166);


/*
  n=200
 */
%rmst_sim(11, 1201, 200, 0.7, 0.5067, 0.4262);
%rmst_sim(12, 1202, 200, 0.7, 0.5067, 0.4262);
%rmst_sim(13, 1203, 200, 0.7, 0.5067, 0.4262);
%rmst_sim(14, 1204, 200, 0.7, 0.5067, 0.4262);
%rmst_sim(15, 1205, 200, 0.7, 0.5067, 0.4262);

%rmst_sim(21, 2201, 200, 0.7, 1.6881, 0.6767);
%rmst_sim(22, 2202, 200, 0.7, 1.6881, 0.6767);
%rmst_sim(23, 2203, 200, 0.7, 1.6881, 0.6767);
%rmst_sim(24, 2204, 200, 0.7, 1.6881, 0.6767);
%rmst_sim(25, 2205, 200, 0.7, 1.6881, 0.6767);

%rmst_sim(31, 3201, 200, 0.7, 8.5230, 0.8785);
%rmst_sim(32, 3202, 200, 0.7, 8.5230, 0.8785);
%rmst_sim(33, 3203, 200, 0.7, 8.5230, 0.8785);
%rmst_sim(34, 3204, 200, 0.7, 8.5230, 0.8785);
%rmst_sim(35, 3205, 200, 0.7, 8.5230, 0.8785);

%rmst_sim(41, 4201, 200, 1.0, 0.6213, 0.4971);
%rmst_sim(42, 4202, 200, 1.0, 0.6213, 0.4971);
%rmst_sim(43, 4203, 200, 1.0, 0.6213, 0.4971);
%rmst_sim(44, 4204, 200, 1.0, 0.6213, 0.4971);
%rmst_sim(45, 4205, 200, 1.0, 0.6213, 0.4971);

%rmst_sim(51, 5201, 200, 1.0, 1.4427, 0.7213);
%rmst_sim(52, 5202, 200, 1.0, 1.4427, 0.7213);
%rmst_sim(53, 5203, 200, 1.0, 1.4427, 0.7213);
%rmst_sim(54, 5204, 200, 1.0, 1.4427, 0.7213);
%rmst_sim(55, 5205, 200, 1.0, 1.4427, 0.7213);

%rmst_sim(61, 6201, 200, 1.0, 4.4814, 0.8963);
%rmst_sim(62, 6202, 200, 1.0, 4.4814, 0.8963);
%rmst_sim(63, 6203, 200, 1.0, 4.4814, 0.8963);
%rmst_sim(64, 6204, 200, 1.0, 4.4814, 0.8963);
%rmst_sim(65, 6205, 200, 1.0, 4.4814, 0.8963);

%rmst_sim(71, 7201, 200, 1.5, 0.7281, 0.5850);
%rmst_sim(72, 7202, 200, 1.5, 0.7281, 0.5850);
%rmst_sim(73, 7203, 200, 1.5, 0.7281, 0.5850);
%rmst_sim(74, 7204, 200, 1.5, 0.7281, 0.5850);
%rmst_sim(75, 7205, 200, 1.5, 0.7281, 0.5850);

%rmst_sim(81, 8201, 200, 1.5, 1.2768, 0.7739);
%rmst_sim(82, 8202, 200, 1.5, 1.2768, 0.7739);
%rmst_sim(83, 8203, 200, 1.5, 1.2768, 0.7739);
%rmst_sim(84, 8204, 200, 1.5, 1.2768, 0.7739);
%rmst_sim(85, 8205, 200, 1.5, 1.2768, 0.7739);

%rmst_sim(91, 9201, 200, 1.5, 2.7182, 0.9166);
%rmst_sim(92, 9202, 200, 1.5, 2.7182, 0.9166);
%rmst_sim(93, 9203, 200, 1.5, 2.7182, 0.9166);
%rmst_sim(94, 9204, 200, 1.5, 2.7182, 0.9166);
%rmst_sim(95, 9205, 200, 1.5, 2.7182, 0.9166);



