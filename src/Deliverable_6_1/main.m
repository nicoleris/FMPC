clc
clear all
close all

quad = Quad();
CTRL = ctrl_NMPC(quad);
sim = quad.sim(CTRL)
quad.plot(sim)
