clear;
syms u v dudy nu

A = [u 0  0; 1 0 0 ; 0 0 0];
B = [v 0 -nu; 0 1 0; 1 0 0];
[V D] = eig(inv(B)*A)