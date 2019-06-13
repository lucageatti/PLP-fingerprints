/**
 * Fingerprint recognition encoding in LPAD.
 * 
 * Francesco Fabiano (University of Udine),
 * Luca Geatti (University of Udine, FBK-Fondazione Bruno Kessler).
*/

:- use_module(library(pita)).
:- pita.

:- begin_lpad.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUX PREDICATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Symmetry breaking:
% we take only left-to-right edges.
:- edge(X1,Y1,X2,Y2),
   minutia(X1,Y1,_,_),
   minutia(X2,Y2,_,_),
   X1>X2.

:- edge(X1,Y1,X2,Y2),
   minutia(X1,Y1,_,_),
   minutia(X2,Y2,_,_),
   X1=X2, Y1>=Y2.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STRUCTURE CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constraint 1:
% Each B-minutia has exactly 3 incident edges.
:- aggregate_all(count, edge(X,Y,X1,Y1), Count_1),
   aggregate_all(count, edge(X2,Y2,X,Y), Count_2),
   Sum is Count_1+Count_2, Sum =\= 3,
   minutia(X,Y,_,b),
   minutia(X1,Y1,_,_),
   minutia(X2,Y2,_,_).

% Constraint 2:
% Each E-minutia has exactly 2 incident edges.
:- aggregate_all(count, edge(X,Y,X1,Y1), Count_1),
   aggregate_all(count, edge(X1,Y1,X,Y), Count_2),
   Sum is Count_1+Count_2, Sum =\= 1,
   minutia(X,Y,_,e),
   minutia(X1,Y1,_,_).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROBABILITY RULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%DIRECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constraint TA1:
%edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T),
%               				 minutia(X2,Y2,D2,T),
%    						 Weight is float(rdiv((pi - (abs(abs(D1-D2) - pi))),pi)).

archs :- type_fingerprint(tented_archs).
archs :- type_fingerprint(plain_archs).
weight_rule(1) :- archs.
weight_rule(W,Y1,Y2) :- type_fingerprint(left_loop),
    			  		minutia(_,Y1,_,_),
    			  		minutia(_,Y2,_,_),
    			  		W is float(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY)).
weight_rule(W,Y1,Y2) :- type_fingerprint(right_loop),
    			  		minutia(_,Y1,_,_),
    			  		minutia(_,Y2,_,_),
    			  		W is float(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY)).

%Constraint TA1-MOD:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T),
    						 minutia(X2,Y2,D2,T),
    						 Alpha is float(atan2(abs(X1-X2),abs(Y1-Y2))),
    						 Diff1 is float(abs(D1-Alpha)),
    						 Diff2 is float(abs(D2-Alpha)),
    						 weight_rule(W,Y1,Y2),
    						 Weight is float(W*(1 - rdiv((Diff1+Diff2),(2*pi)))).

%%Calculate the line angle and then make D1,D2 relative to it
% Constraint TA2:
%edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T1),
%               				 minutia(X2,Y2,D2,T2),
%    						 T1 \== T2,
%    						 Half_pi is rdiv(pi,2),
%    						 Weight is float(1-rdiv((Half_pi - abs((abs(abs(abs(D1-D2) - Half_pi)-Half_pi)-Half_pi))),Half_pi)).


%Constraint TA2-MOD:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T1),
    						 minutia(X2,Y2,D2,T2),
    						 T1 \== T2,
    						 Alpha is float(atan2(abs(X1-X2),abs(Y1-Y2))),
    						 Diff1 is float(abs(D1-Alpha)),
    						 Diff2 is float(abs(D2-Alpha)),
    						 weight_rule(W,Y1,Y2),
    						 Weight is float(W*(1 - rdiv((Diff1+Diff2),(2*pi)))).


%%%%%%%%%%%%%%%%%%%%%%DISTANCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constraint TA4:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,_),
               				 minutia(X2,Y2,_,_),
    						 max_Y(MY),
    						 archs,
    						 Weight is float(1 - rdiv(abs(Y1-Y2),MY)).

%Constraint TA5:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,b),
               				 minutia(X2,Y2,_,b),
        					 max_X(MX),
                             archs,
    						 Weight is float(1 - rdiv(abs(X1-X2),MX)).

%Constraint TA6:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
        					 max_X(MX),
    						 archs,
    						 Weight is float(1 - rdiv(abs(rdiv(abs(X1-X2),2) - rdiv(MX,2)), MX)).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEFT LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TODO:
%  1. right loop
%  2. contraint on X for loops
%  3. TA2 correct angle
%%%
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
    						 type_fingerprint(left_loop),
        				     max_Y(MY),
          				     max_X(MX),
        					 Weight_ruleY is float(1-(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY))),
                             Weight_rule is float(Weight_ruleY*(((1 - rdiv(X1,MX)))*abs(1 - rdiv(X2,MX)))),
    						 Weight is float(Weight_rule*(1 - rdiv(rdiv(abs(Y1-Y2),2) - rdiv(MY,4), MY))).
    
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
    						 type_fingerprint(left_loop),
        				     max_Y(MY),
        					 Weight_rule is float(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY)),
    						 Weight is float(Weight_rule*(1- rdiv(abs(Y1-Y2),MY))).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KNOWLEDGE BASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_Y(450).
max_X(338).
type_fingerprint(tended_archs).
minutia(60,84,0.14888994760949725,e).
minutia(47,166,5.734635904729058,e).
minutia(39,302,0.8478169733934057,e).
minutia(58,256,5.31793364398966,b).
minutia(94,104,0.09966865249116204,b).
minutia(62,272,5.135242906517631,b).
minutia(88,292,5.11091142605111,e).
minutia(116,292,5.9917285127017195,b).
minutia(133,346,0.3028848683749714,e).
minutia(158,256,3.2904826011992903,b).
minutia(163,270,0.049958395721942765,b).
minutia(193,122,3.3865713167166573,b).
minutia(191,150,3.141592653589793,e).
minutia(171,254,3.5644465797227336,b).
minutia(193,170,6.134295359570089,b).
minutia(191,222,0.5404195002705842,e).
minutia(196,252,3.2904826011992903,b).
minutia(218,232,3.7670776938290222,e).
minutia(214,280,2.7187387274568526,b).
minutia(62,35,2.6387494426619322,e).
minutia(82,355,0.5880026035475675,b).
minutia(169,117,3.3389882134396736,e).
minutia(190,59,3.682012153860377,e).
minutia(165,245,0.33667481938672716,b).
minutia(42,198,2.4186133001883023,b).
/* */
minutia(72,58,0.14888994760949725,e).
minutia(45,232,5.526350801620898,e).
/* */
minutia(97,152,2.761086276477428,e).
minutia(79,276,1.9513027039072617,e).
minutia(91,252,2.473303234759209,b).
minutia(90,262,4.9786410295356145,e).
minutia(101,256,2.2531128816696446,e).
minutia(111,292,5.465540261346884,b).
minutia(179,236,3.5644465797227336,b).
minutia(203,180,0.4228539261329407,e).
minutia(226,126,3.690142056040321,e).
minutia(205,304,2.7187387274568526,b).
minutia(201,344,2.8809902618424523,e).
minutia(216,290,2.743070207923373,b).
minutia(49,99,5.407127256581393,b).
minutia(49,107,4.2487413713838835,e).
minutia(107,311,0.049958395721942765,e).
minutia(223,199,3.7523186179790016,e).
minutia(58,76,0.14888994760949725,e).
minutia(38,292,0.9380474917927134,e).
minutia(86,348,3.47826747297652,e).
minutia(166,124,6.233226911457644,e).
minutia(169,150,0.0,e).
minutia(170,180,3.141592653589793,b).
minutia(175,282,3.644435864517654,b).
minutia(163,346,3.43304944806766,e).
minutia(190,212,3.5644465797227336,b).
minutia(190,236,3.283489708193957,e).
minutia(230,172,0.6350267353903137,e).
minutia(216,320,5.780342096251726,e).
minutia(133,313,0.14888994760949725,b).
minutia(165,185,0.2914567944778671,b).
minutia(167,333,0.0,b).
minutia(225,137,3.4474715249946453,e).
minutia(30,138,2.9147938055359073,b).
minutia(62,50,6.134295359570089,b).
minutia(41,156,3.0916342578678506,b).
minutia(51,240,5.407127256581393,e).
minutia(48,262,3.776619388980107,b).
minutia(104,54,2.992702705980296,b).
minutia(96,110,3.0916342578678506,b).
minutia(64,318,3.7295952571373605,b).
minutia(123,64,0.14888994760949725,b).
minutia(123,80,3.1915510493117356,b).
minutia(158,146,3.0916342578678506,b).
minutia(148,270,0.4636476090008061,b).
minutia(148,286,3.141592653589793,b).
minutia(221,164,3.246469592320027,e).
minutia(212,270,2.677945044588987,b).
minutia(232,230,5.135242906517631,e).
minutia(100,241,5.742765806909002,e).
minutia(165,55,3.3865713167166573,b).
minutia(205,319,5.672459342790377,e).

:- end_lpad.