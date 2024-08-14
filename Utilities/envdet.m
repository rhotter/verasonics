% function [y1,y2]=envdet(a,x);
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%  * envdet: nonlinear filter, first order envelope follower.
%  * envdet.c, envdet.m mexfile
%  * a - forgetting factor
%  * x - input data row vector
%  * y1,y2 - outputs (order dependent, see below)
%  *
%  * of a 1xN matrix (inMatrix)
%  * outputs a 1xN matrix (outMatrix)
%  *
%  * The calling syntax is:
%  *
%  *>>		[yrev] = envdet(alpha , x);
%  *>>		[yfwd,yrev] = envdet(alpha , x);
%  *compute forward and reverse envelopes, where
%  * forward (causal) envelope is:
%  * for k=1:N,
%  *  y(k) <-- x(k)*alpha + y(k-1)*(1-alpha);
%  *  y(k) <-- max(y(k),x(k))
%  * Reverse (anticausal) defined similarly.
%  * Testing:
%  * >> x = [cos(linspace(-1,1,10)');.000001*ones(10,1)]*rand(1,10).^4;x=x(:).';
%  * >> [y,y2]=envdet(.4,x);
%  * >> plot(1:length(x),db(x)/2,'.',1:length(x),db(y)/2,'-r',1:length(x),db(y2)/2,'-g')
%  *
%  * ******************************** *
%  * John Flynn 11/24/2009, based on 1995 aud2midi algorithm.
%  *========================================================*/
