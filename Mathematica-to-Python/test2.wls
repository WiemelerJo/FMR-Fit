#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
target=Rest@$ScriptCommandLine;
spitTarget=StringSplit[target,"|"];
{inputstringA,inputstringB,inputstringC}=Flatten@Table[s,{s,spitTarget}];
{inputstringA,inputstringB,inputstringC}={inputstringA,ToExpression[inputstringB],ToExpression[inputstringC]};
Print[a]

