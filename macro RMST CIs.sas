/**************************************************************************************************
programmer    :Hiroya Hashimioto
file name     :macro RMST CIs.sas
date created  :2020/08/27
     modified :
description   :SAS macro to estimate some confidence intervals considered for RMST.
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


/***************************************************************
 Macro section
***************************************************************/



%macro rmst_ci(data,       /* source data */
                time,       /* time to event */
                status,     /* variable indicating event or censored */
                cen_ind,    /* value indicating censoring */
                tau,        /* truncation time */
                alpha,      /* gives (1-alpha) confidence interval */
                output     /* the output dataset name */
                );

data mac1;
  set &data;
  t_tau=&time/&tau;
  if t_tau>1 then &status=&cen_ind;


ods listing close;
proc lifetest data=mac1 timelim=1;
  time t_tau * &status(&cen_ind);
  ods output Means=mac2
             ProductLimitEstimates=mac3(where=(survival ne . and t_tau <= 1));
run;
ods listing;

proc sort data=mac3 out=mac4(keep= failed);
  by descending failed;
run;

data mac5;
  set mac4(obs=1);
run;

data mac6;
  merge mac2 mac5;
run;

%let c= 1-&alpha./2;

data &output(drop=stratum mean stderr bias failed se_km);
  set mac6;

  Estimate = Mean*&tau;

  SE_KM=StdErr*sqrt((failed-1)/failed);

  NT_Gw_L    = (Mean - SE_KM*quantile('NORMAL',&c.))*&tau;
  NT_Gw_U    = (Mean + SE_KM*quantile('NORMAL',&c.))*&tau;
  CLL_Gw_L   = (exp(-exp( log(-log(Mean)) - quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*SE_KM )))*&tau;
  CLL_Gw_U   = (exp(-exp( log(-log(Mean)) + quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*SE_KM )))*&tau;
  ASR_Gw_L   = ((sin( arsin(sqrt(Mean)) - quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*SE_KM ))**2)*&tau;
  ASR_Gw_U   = ((sin( arsin(sqrt(Mean)) + quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*SE_KM ))**2)*&tau;
  Logit_Gw_L = (-1/(1+exp( log(Mean/(1-Mean)) - quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*SE_KM )) +1)*&tau;
  Logit_Gw_U = (-1/(1+exp( log(Mean/(1-Mean)) + quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*SE_KM )) +1)*&tau;

  NT_KM_L    = (Mean - StdErr*quantile('NORMAL',&c.))*&tau;
  NT_KM_U    = (Mean + StdErr*quantile('NORMAL',&c.))*&tau;
  CLL_KM_L   = (exp(-exp( log(-log(Mean)) - quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*StdErr )))*&tau;
  CLL_KM_U   = (exp(-exp( log(-log(Mean)) + quantile('NORMAL',&c.)*(1/log(Mean))*(1/Mean)*StdErr )))*&tau;
  ASR_KM_L   = ((sin( arsin(sqrt(Mean)) - quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*StdErr ))**2)*&tau;
  ASR_KM_U   = ((sin( arsin(sqrt(Mean)) + quantile('NORMAL',&c.)*(1/sqrt(1-Mean))*(1/(2*sqrt(Mean)))*StdErr ))**2)*&tau;
  Logit_KM_L = (-1/(1+exp( log(Mean/(1-Mean)) - quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*StdErr )) +1)*&tau;
  Logit_KM_U = (-1/(1+exp( log(Mean/(1-Mean)) + quantile('NORMAL',&c.)*(1/(Mean*(1-Mean)))*StdErr )) +1)*&tau;

run;

proc datasets lib=work nolist;
  delete mac:;
run;quit;

%mend;


data bmt_all;
  set sashelp.bmt;
  where group="ALL";
run;

%rmst_ci(bmt_all, T, Status, 0, 1000, 0.05, result);

