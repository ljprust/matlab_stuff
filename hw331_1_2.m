% Logan Prust - Aer E 331 - Homework 1 - Problem 2

clear,clc;

A=[.0997 -.9933;.5101 -.2452];
B=[.0997;.5101];
C=eye(2);
D=[0;0];
[NUM,DEN]=ss2tf(A,B,C,D)