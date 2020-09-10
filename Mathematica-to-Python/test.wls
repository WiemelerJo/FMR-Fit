#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
results=Table[RandomReal[],{10}];
target=Rest@$ScriptCommandLine;
(*target="/home/benjamin/Documents/test/output/results.csv";*)
Export[target,results];
Print["target = "  <> target]
Print[results]
