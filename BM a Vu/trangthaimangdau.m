function trangthaimangdau
clc;
clf;

C   = 0;
D   = 0;
e   = 0;
Fr  = 0;
Ft  = 0;
h   = 0;
M   = 0;
p   = 0;
R   = 0;    %Ban kinh ngong truc
T   = 0;    % Thoi gian vo huong
t   = 0;    % Thoi gian
W   = 0;    % Tai trong
W0  = 0;    % Tai trong tinh
epxilon     = 0;
neta        = 0;
phi         =0;
phihoa1     = 0;
phihoa2     =0;
teta1       = pi() + phihoa1;
teta2       = pi + phihoa2;
Om          = OmegaP/ Omega;
Omega       = 0;    % Van toc goc c?a Truc.
OmegaP      = 0;

function [epxilon, phi] = tinh(
 
