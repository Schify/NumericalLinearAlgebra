clear all
addpath("../..")
addpath("../../matlab2tikz/src/")
addpath("../NLAHW1_Lanczos/")
addpath("../NLAHW1_AFM/")
close all
load("../NLAHW1_AFM/HW1.mat")

spy(max(abs(full(M'-M))))