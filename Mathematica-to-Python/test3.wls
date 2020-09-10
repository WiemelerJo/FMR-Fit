#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
target=Rest@$ScriptCommandLine;
target=Rest@$ScriptCommandLine;
splitTarget=StringSplit[target,"|"];
{inputstringF,inputstringFilePath,inputstringParams,inputstringStartVal,inputstringRanges,inputstringSteps,inputstringfixedParams,inputstringfixedValues,inputstringanglestep}=Flatten@Table[s,{s,splitTarget}];
Print[inputstringF]
