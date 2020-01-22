(* ::Package:: *)

(* ::Input:: *)
(*k= 1; a=1;*)
(*keff = 1+1I;*)
(*\[Theta]in = 3/10;*)
(*A[n_,x_] = I^n E^(-I n \[Theta]in) E^(I x keff);*)
(**)
(*F[n_,X_,Y_]= (-1)^n Exp[I n ArcTan[X,Y] ]HankelH1[n,Sqrt[X^2+Y^2]];*)
(**)
(*ClearAll[integrateA]*)
(*integrateA[m_,x_,\[Theta]in_:0, M_:0]:=Sum[*)
(* NIntegrate[*)
(*Boole[(x2- x)^2+ y2^2>Abs[k]^2a^2 ]A[n,x2]Exp[I y2 Sin[\[Theta]in]]F[n-m, x2- x, y2]*)
(*,{y2,-\[Infinity],\[Infinity]},{x2,0,\[Infinity]}*)
(*,(*MaxRecursion\[Rule]12,*) PrecisionGoal->6]*)
(*,{n,-M,M}]*)
(*x1s = {0,1,2}*)


(* ::Input:: *)
(*data = Outer[(d=integrateA[#2,#1,\[Theta]in,3];Print["{x1s,\[Theta]} = ",{#1,#2},", result = ",d]; d)&,x1s,Range[-3,3,1]]     *)
(**)
(*Export[Directory[]<>"/Integrated_As.csv",data]*)


im=I;
 {
{65.7258-147.667im, 24.1063-8.37573im, -5.18432+3.75268im, 0.926351-3.3608im, 7.32809+0.00509993im, 22.9749+38.4621im, -50.276+91.6842im}
 ,{-32.8167+1.54958im, -12.1046-16.272im, 5.91279-10.8138im, 8.56817+3.75119im, -4.88707+9.31373im, -12.1282-5.56668im, 3.68169-13.925im}
 ,{-24.2364-34.8114im, 10.5494-15.7369im, 12.0539+6.19178im, -3.0705+11.7496im, -12.7322-1.11301im, 3.08584-14.5042im, 21.713+4.88521im}
}
