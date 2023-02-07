function [VAR] = mean_var(fname, nam, var_string)

%% ========================================================================
% read all relevant 3-D variables

var=ncread(fname,var_string);
s=size(var); n=s(1); m=s(2); l=s(3); nm=n*m;

VAR=mean(reshape(var,nm,l));                                           

% These here in case needed later
% varpr(1,1,:)=VAR;
% varpert=bsxfun(@minus,var,varpr);

