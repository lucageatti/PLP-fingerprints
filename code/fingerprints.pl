/*
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
%%% STRUCTURE PREDICATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
archs:- type_fingerprint(tented_archs).
archs:- type_fingerprint(plain_archs).
loops:- type_fingerprint(left_loop).
loops:- type_fingerprint(right_loop).

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

%%%%DIRECTION

%%%Weights
weight_rule(W,_,_) :- archs,
    				  W is 1.
weight_rule(W,Y1,Y2) :- loops,
    			  		minutia(_,Y1,_,_),
    			  		minutia(_,Y2,_,_),
    			  		W is float(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY)).

%Constraint POS_1:
%edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T),
%    						 minutia(X2,Y2,D2,T),
% 							%Alpha1 is float(atan2(abs(Y1-Y2),abs(X1-X2))),
%            			     Alpha1 is float(atan2(Y1-Y2,X1-X2)),
%        					 Alpha2 is float(atan2(Y2-Y1,X2-X1)),
%    						 %Diff1 is float(abs(D1-Alpha1)),
%    						 %Diff2 is float(abs(D2-Alpha1)),
%    						 Diff1 is float(pi - (abs(abs(D1-Alpha1) - pi))),
%        					 Diff2 is float(pi - (abs(abs(D2-Alpha2) - pi))),
%    						 weight_rule(W,Y1,Y2),
%    						 Weight is float(W*(1 -rdiv((Diff1+Diff2),(2*pi)))).
%OLD
%edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T),
%               			  minutia(X2,Y2,D2,T),
%    						  Weight is float(rdiv((pi - (abs(abs(D1-D2) - pi))),pi)).

% Constraint POS_2:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,_),
    						 minutia(X2,Y2,D2,_),
    						 %Alpha1 is float(atan2(abs(Y1-Y2),abs(X1-X2))),
            			     Alpha1 is float(atan2(Y1-Y2,X1-X2)),
        					 Alpha2 is float(atan2(Y2-Y1,X2-X1)),
    						 %Diff1 is float(abs(D1-Alpha1)),
    						 %Diff2 is float(abs(D2-Alpha1)),
    						 Diff1 is float(abs(D1-Alpha1)),
        					 Diff2 is float(abs(D2-Alpha2)),
    						 weight_rule(W,Y1,Y2),
    						 Weight is float(W*(1 -rdiv((Diff1+Diff2),(2*pi)))).
%OLD
%edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,T1),
%               			  minutia(X2,Y2,D2,T2),
%    						  T1 \== T2,
%    						  Half_pi is rdiv(pi,2),
%    						  Weight is float(1-rdiv((Half_pi - abs((abs(abs(abs(D1-D2) - Half_pi)-Half_pi)-Half_pi))),Half_pi)).


%%%%DISTANCE

%%%ARCHS
%Constraint DIST_ARCH_1:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,_),
               				 minutia(X2,Y2,_,_),
    						 max_Y(MY),
    						 archs,
    						 Weight is float(1 - rdiv(abs(Y1-Y2),MY)).

%Constraint DIST_ARCH_2:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,b),
               				 minutia(X2,Y2,_,b),
        					 max_X(MX),
                             archs,
    						 Weight is float(1 - rdiv(abs(X1-X2),MX)).

%Constraint DIST_ARCH_3:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
        					 max_X(MX),
    						 archs,
    						 Weight is float(1 - rdiv(abs(rdiv(abs(X1-X2),2) - rdiv(MX,2)), MX)).

%%%LOOPS

%%%Left loop
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
    						 type_fingerprint(left_loop),
        				     max_Y(MY),
          				     max_X(MX),
        					 Weight_ruleY is float(1-(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY))),
                             Weight_rule is float(Weight_ruleY*(((1 - rdiv(X1,MX)))*abs(1 - rdiv(X2,MX)))),
    						 Weight is float(Weight_rule*(1 - rdiv(rdiv(abs(Y1-Y2),2) - rdiv(MY,4), MY))).
    
%%Right loop
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
    						 type_fingerprint(right_loop),
        				     max_Y(MY),
          				     max_X(MX),
        					 Weight_ruleY is float(1-(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY))),
                             Weight_rule is float(Weight_ruleY*((rdiv(X1,MX))*abs(rdiv(X2,MX)))),
    						 Weight is float(Weight_rule*(1 - rdiv(rdiv(abs(Y1-Y2),2) - rdiv(MY,4), MY))).

%loops
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
    						 loops,
        				     max_Y(MY),
        					 Weight_rule is float(abs(1 - rdiv(Y1,MY)))*abs(1 - rdiv(Y2,MY)),
    						 Weight is float(Weight_rule*(1- rdiv(abs(Y1-Y2),MY))).


:- mostProb(X1,Y1,X2,Y2),
   edge(X1,Y1,X2,Y2),
   prob(edge(X1,Y1,X2,Y2),P1),
   prob(edge(_,_,_,_),P2),
   P1 < P2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KNOWLEDGE BASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_Y(450).
max_X(338).
type_fingerprint(tented_archs).
minutia(338,84,3.14888994760949725,e).
minutia(47,450,0.1334635904729058,e).
minutia(39,302,0.8478169733934057,e).
minutia(58,256,5.31793364398966,b).
minutia(94,104,0.09966865249116204,b).
minutia(62,272,5.135242906517631,b).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QUERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prob edges
%prob(edge(X1,Y1,X2,Y2),P).

:- end_lpad.