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
good_edge(X1,Y1,X2,Y2) :-
   minutia(X1,Y1,_,_),
   minutia(X2,Y2,_,_),
   X1<X2.
good_edge(X1,Y1,X2,Y2) :-
   minutia(X1,Y1,_,_),
   minutia(X2,Y2,_,_),
   X1=X2,Y1<Y2.

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

% SA1:
% Each B-minutia has exactly 3 incident edges.
b_inc_edge :-  aggregate_all(count, edge(X,Y,X1,Y1), Count1),
               aggregate_all(count, edge(X2,Y2,X,Y), Count2),
               Sum is Count1+Count2, Sum == 3,
               minutia(X,Y,_,b),
               minutia(X1,Y1,_,_),
               minutia(X2,Y2,_,_).

% SA2:
% Each E-minutia has exactly 1 incident edge.
e_inc_edge :-  aggregate_all(count, edge(X,Y,X1,Y1), Count1),
               aggregate_all(count, edge(X1,Y1,X,Y), Count2),
               Sum is Count1+Count2, Sum == 1,
               minutia(X,Y,_,e),
               minutia(X1,Y1,_,_).

valid_graph :- b_inc_edge, e_inc_edge.
valid_graph.

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
    					max_Y(MY),
    			  		W is float(abs(1 - rdiv(Y1,MY))*abs(1 - rdiv(Y2,MY))).


% Constraint POS:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,D1,_),
    						 minutia(X2,Y2,D2,_),
    						 good_edge(X1,Y1,X2,Y2),
    						 %Alpha1 is float(atan2(abs(Y1-Y2),abs(X1-X2))),
            			     Alpha1 is float(atan2(Y1-Y2,X1-X2)),
        					 Alpha2 is float(atan2(Y2-Y1,X2-X1)),
    						 %Diff1 is float(abs(D1-Alpha1)),
    						 %Diff2 is float(abs(D2-Alpha1)),
    						 Diff1 is float(abs(D1-Alpha1)),
        					 Diff2 is float(abs(D2-Alpha2)),
    						 weight_rule(W,Y1,Y2),
    						 Weight is float(W*(1 -rdiv((Diff1+Diff2),(2*pi)))).

%%%%DISTANCE

%%%ARCHS
%Constraint DIST_ARCH_1:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,_),
               				 minutia(X2,Y2,_,_),
    						 good_edge(X1,Y1,X2,Y2),
    						 max_Y(MY),
    						 archs,
    						 Weight is float(1 - rdiv(abs(Y1-Y2),MY)).

%Constraint DIST_ARCH_2:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,b),
               				 minutia(X2,Y2,_,b),
                             good_edge(X1,Y1,X2,Y2),
        					 max_X(MX),
                             archs,
    						 Weight is float(1 - rdiv(abs(X1-X2),MX)).

%Constraint DIST_ARCH_3:
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
                             good_edge(X1,Y1,X2,Y2),
        					 max_X(MX),
    						 archs,
    						 Weight is float(1 - rdiv(abs(rdiv(abs(X1-X2),2) - rdiv(MX,2)), MX)).

%%%LOOPS

%%% CONSTRAINT LL 
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
                             good_edge(X1,Y1,X2,Y2),
    						 type_fingerprint(left_loop),
        				     max_Y(MY),
          				     max_X(MX),
        					 weight_rule(Weight_ruleY,Y1,Y2),
                             Weight_rule is float((1-Weight_ruleY)*(((1 - rdiv(X1,MX)))*abs(1 - rdiv(X2,MX)))),
    						 Weight is float(Weight_rule*(1 - rdiv(abs(rdiv(abs(Y1-Y2),2) - 3*rdiv(MY,4)), MY))).
    
%%% CONSTRAINT RL 
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
                             good_edge(X1,Y1,X2,Y2),
    						 type_fingerprint(right_loop),
        				     max_Y(MY),
          				     max_X(MX),
        					 weight_rule(Weight_ruleY,Y1,Y2),
                             Weight_rule is float((1-Weight_ruleY)*((rdiv(X1,MX))*abs(rdiv(X2,MX)))),
    						 Weight is float(Weight_rule*(1 - rdiv(abs(rdiv(abs(Y1-Y2),2) - 3*rdiv(MY,4)), MY))).

%%% CONSTRAINT LOOP
edge(X1,Y1,X2,Y2): Weight :- minutia(X1,Y1,_,e),
               				 minutia(X2,Y2,_,e),
                             good_edge(X1,Y1,X2,Y2),
    						 loops,
        				     max_X(MX),
        					 weight_rule(Weight_rule,Y1,Y2),
    						 Weight is float((1-Weight_rule)*(1- rdiv(abs(X1-X2),MX))).




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KNOWLEDGE BASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_X(504).
max_Y(480).
type_fingerprint(plain_archs).
minutia(126,224,3.38657,e).
minutia(229,162,3.89997,e).
minutia(184,436,4.37571,e).
minutia(180,319,0.896055,e).
minutia(130,366,2.65711,e).
minutia(168,242,0.576375,e).
minutia(132,458,1.00748,e).
minutia(141,416,1.03038,e).
minutia(162,388,3.92699,b).
minutia(164,426,4.13601,e).
minutia(223,204,3.81633,e).
minutia(226,454,4.13601,e).
minutia(50,335,2.07364,b).
minutia(67,284,6.08579,b).
minutia(89,318,2.56522,e).
minutia(115,386,5.99173,b).
minutia(122,356,3.04192,e).
minutia(145,382,2.9927,e).
minutia(177,289,3.70491,b).
minutia(188,432,1.10715,e).
minutia(90,214,0.197396,b).
minutia(91,320,6.03821,e).
minutia(108,344,0,b).
minutia(121,308,3.09163,e).
minutia(130,419,1.32582,b).
minutia(193,164,3.71797,e).
minutia(162,352,3.75232,b).
minutia(209,263,3.81633,b).
minutia(135,387,0.380506,e).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QUERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prob edges
%prob(edge(X1,Y1,X2,Y2),P).
%
%Find all Edges with related probability
%findall([X1,Y1,X2,Y2,P],prob(edge(X1,Y1,X2,Y2),P),Results).


edge(162, 352, 162, 388).
edge(91, 320, 180, 319).
edge(115, 386, 162, 388).
edge(89, 318, 180, 319).
edge(108, 344, 115, 386).
edge(193, 164, 229, 162).
edge(50, 335, 108, 344).
edge(108, 344, 162, 352).
edge(108, 344, 130, 419).
edge(135, 387, 162, 388).
edge(184, 436, 188, 432).
edge(115, 386, 130, 419).
edge(130, 419, 141, 416).
edge(115, 386, 135, 387).
edge(177, 289, 209, 263).
edge(184, 436, 188, 432).
edge(89, 318, 91, 320).
edge(132, 458, 226, 454).
edge(130, 419, 162, 388).
edge(50, 335, 67, 284).
edge(162, 352, 177, 289).
edge(67, 284, 177, 289).
edge(135, 387, 145, 382).
edge(122, 356, 162, 352).
edge(90, 214, 108, 344).
edge(50, 335, 90, 214).
edge(89, 318, 121, 308).
edge(67, 284, 90, 214).
edge(162, 388, 177, 289).
edge(164, 426, 188, 432).
edge(130, 419, 162, 352).
edge(115, 386, 145, 382).
edge(115, 386, 162, 352).
edge(108, 344, 162, 388).
edge(50, 335, 180, 319).

:- end_lpad.