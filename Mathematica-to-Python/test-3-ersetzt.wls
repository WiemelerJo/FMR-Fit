#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
target=Rest@$ScriptCommandLine;
splitTarget=StringSplit[target,"|"];
{inputstringF,inputstringFilePath,inputstringParams,inputstringStartVal,inputstringRanges,inputstringSteps,inputstringfixedParams,inputstringfixedValues,inputstringanglestep,inputstringIterations,inputstringOutputPath,inputstringBresColumn,inputstringWinkelColumn,inputstringShift}=Flatten@Table[s,{s,splitTarget}];
(*F[B_,thetaB_,phiB_,M_,theta_,phi_]=M*B(Sin[theta]*Sin[thetaB]*Cos[phi-phiB]+Cos[theta]*Cos[thetaB])-(1/2*\[Mu]0*M^2-K2s)Sin[theta]^2-K2p*Sin[theta]^2*Cos[phi-phiu]^2-1/2*K4s*Cos[theta]^4-1/8*K4p(3+Cos[4*phi])Sin[theta]^4;*)
ToExpression[inputstringF];\[AliasDelimiter]
hl=Simplify[1/(M^2 Sin[theta]^2) (D[F[B,thetaB,phiB,M,theta,phi],{theta,2}]D[F[B,thetaB,phiB,M,theta,phi],{phi,2}]-(D[F[B,thetaB,phiB,M,theta,phi],theta,phi])^2)];
ResField[thetaB_,phiB_,M_,theta_,phi_]=(B/.Solve[hl==(omega/gamma)^2,B])[[2]];(*Not sure if Number 2 is automatically the physically valid solution*)

SolveAngles[B_,thetaB_,phiB_,rule_]:=(Sort[Table[Quiet@FindMinimum[{F[B,thetaB,phiB,M,theta,phi]/.rule},{{phi,phiB-n \[Pi]/8.},{theta,thetaB-m \[Pi]/8.}}],{n,-1,1},{m,-1,1}],#1[[1]]<#2[[1]]&][[1,2]]);
angleFunc[B_,thetaB_,phiB_,rule_]:={theta,phi}/.(SolveAngles[B,thetaB,phiB,rule][[2]]);

ResFieldNumInp[thetaB_,phiB_,rule_]:=(ResField[thetaB,phiB,M,#[[1]],#[[2]]]/.rule)&[angleFunc[BinterAligned[phiB],thetaB,phiB,rule]];
para = Import[inputstringFilePath,"Table"];
WinkelString = ToExpression[inputstringWinkelColumn];
BresString = ToExpression[inputstringBresColumn];
data=Transpose[{Table[para[[i,WinkelString]],{i,1,Length[#]}],#}&@Map[#[[BresString]]&,para]][[18;;-20]];

(*data=Sort[Transpose[{Mod[Transpose[data][[1,1;;Length[Transpose[data][[1]]]-1]]-20,180],
Transpose[data][[2,1;;Length[Transpose[data][[2]]]-1]]}]];*)(*Readjust the position of the Anisotropy peaks (might be complicated, depending on how the data looks like)*)
Shift=ToExpression[inputstringShift];
interTab=Sort@DeleteDuplicates@Transpose[{\[Pi]/180 Mod[#[[1]],360],#[[2]]}&@(Transpose[data]+{Shift,0})];
phimin=Min[Transpose[interTab][[1]]];
phimax=Max[Transpose[interTab][[1]]];
BinterAligned=Check[Interpolation[interTab],0];
ToExpression[inputstringParams](*e.g. fitParameters={K2p,K2s,K4p,phiu}*);

(*startvalues={startk2p,startk2s,startk4p,startk4s};*)
ToExpression[inputstringStartVal](*{863.25,261345,13720.6,5.07568359375`}*);
ranges=ToExpression[inputstringRanges](*{0.5,0.5,0.5,0.5}*);(*The ranges as a fraction of the starting Parameters eg.: 0.5 \[Rule] Try fitting from startvalue - 0.5*starvalue to startvalue + 0.5*startvalue*)
steps=ToExpression[inputstringSteps](*{0.5,0.5,0.5,0.5}*);(*The steps to take while fitting, as a fraction of the ranges*)

fixedParameters=ToExpression[inputstringfixedParams](*{omega,g,M,K4s}*);
fixedParameterValues=ToExpression[inputstringfixedValues](*{2\[Pi]*9.8782*10^9,2.05,1.53*10^6,0}*);

constants={mu0,\[Mu]b,\[HBar]};
constantValues={4*\[Pi]*10^-7,9.27*10^-24,(6.62606957*10^-34)/(2\[Pi])};

angleStep=ToExpression[inputstringanglestep](*\[Pi]/45*);

gamma=(g*\[Mu]b)/\[HBar];
parameterRules=Table[fixedParameters[[i]]->fixedParameterValues[[i]],{i,1,Length[fixedParameters]}];
constantRules=Table[constants[[i]]->constantValues[[i]],{i,1,Length[constants]}];
rules=Flatten[{parameterRules,constantRules}];



(*
pre=Table[{phiB,ResFieldNumInp[\[Pi]/2,phiB,Flatten[{Table[fitParameters[[l]]\[Rule]startvalues[[l]],{l,1,Length[startvalues]}],rules}]]},{phiB,phimin,phimax,angleStep}];
Show[Plot[BinterAligned[phi],{phi,phimin,phimax},GridLines\[Rule]{Table[i,{i,-2\[Pi],2\[Pi],\[Pi]/32}],Automatic},Frame\[Rule]True,PlotRange\[Rule]All],ListPlot[interTab,PlotStyle\[Rule]Black],ListLinePlot[pre,PlotStyle\[Rule]Red,PlotRange\[Rule]All]];*)
iterations=inputstringIterations(*5*);
fitParamsEvolution={};
run=True;
If[Length[constants]!=Length[constantValues],Print["ConstatsTables are of unequal Length. Skipping Evaluation"];run=False]
If[Length[fixedParameters]!=Length[fixedParameterValues],Print["ParameterTables are of unequal Length. Skipping Evaluation"];run=False]
If[Length[fitParameters]!=Length[ranges],Print["Ranges are of unequal Length. Skipping Evaluation"];run=False]
If[Length[fitParameters]!=Length[steps],Print["Steps are of unequal Length. Skipping Evaluation"];run=False]
If[Length[fitParameters]!=Length[startvalues],Print["Startvalues are of unequal Length. Skipping Evaluation"];run=False]
If[run,
For[i=1,i<=iterations,i++,
iterationInstructions=Tuples[Table[Table[fitParameters[[ii]]->l,{l,startvalues[[ii]](1-ranges[[ii]]),startvalues[[ii]](1+ranges[[ii]]),startvalues[[ii]]*ranges[[ii]]*steps[[ii]]}],{ii,1,Length[fitParameters]}]];
difftab=ParallelTable[Sum[Abs[ResFieldNumInp[\[Pi]/2,phiB,Flatten[{it,rules}]]-BinterAligned[phiB]],{phiB,phimin,phimax,angleStep}],{it,iterationInstructions}];
pos=Position[difftab,Min[difftab]][[1,1]];
AppendTo[fitParamsEvolution,iterationInstructions[[pos]]];
startvalues=fitParameters/.iterationInstructions[[pos]];
ranges=0.5*ranges;(*The ranges as a fraction of the starting Parameters eg.: 0.5 \[Rule] Try fitting from startvalue - 0.5*starvalue to startvalue + 0.5*startvalue*)
];
fittedParams=iterationInstructions[[pos]]
]
Min[difftab]/Length[Table[1,{phiB,phimin,phimax,angleStep}]]*(M/.rules)
tmp=ParallelTable[{phiB,ResFieldNumInp[\[Pi]/2,phiB,Flatten[{fittedParams,rules}]]},{phiB,phimin,phimax,angleStep}];
(*Show[ListLinePlot[tmp,PlotStyle\[Rule]Red],Plot[BinterAligned[phiB],{phiB,phimin,phimax}],PlotRange\[Rule]All]*)
tmp2=Table[ParallelTable[{phiB,ResFieldNumInp[\[Pi]/2,phiB,Flatten[{i,rules}]]-BinterAligned[phiB]},{phiB,phimin,phimax,angleStep}],{i,fitParamsEvolution}];
(*ListPlot[tmp2,PlotLegends\[Rule]Table[i,{i,1,Length[fitParamsEvolution]}],PlotRange\[Rule]All]*)
Konstanten = {K2p /. fittedParams,K2s /. fittedParams,K4p /. fittedParams,0,phiu /.fittedParams};
polarplot = Show[PolarPlot[BinterAligned[phiB],{phiB,phimin,phimax}],ListPolarPlot[BinterAligned,PlotRange->All,PlotStyle->Black],ListPolarPlot[tmp,Joined->True,PlotStyle->Red]];
plot = Show[ListPlot[BinterAligned,PlotStyle->Black],Plot[BinterAligned[phiB],{phiB,phimin,phimax}],ListLinePlot[tmp,PlotStyle->Red],PlotRange->All];
outputPath = ToExpression[inputstringOutputPath];
Export[outputPath<>"outputPolarPlot.pdf",polarplot];
Export[outputPath<>"outputPlot.pdf",plot];
Export[outputPath<>"outputAnisotropie-Konstanten.dat",Konstanten];
Export[outputPath<>"outputFunktion.dat",tmp];
