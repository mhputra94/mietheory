% Fungsi Riccati Bessel Tipe 1
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 21 Februari 2015

function [psi, pso] = ricbestu(n,z)
% Output
% psi merupakan fungsi Riccati Bessel Tipe 1
% pso merupakan turunan fungsi RIccati Bessel Tipe 1

% Input
% n orde fungsi Riccati Bessel
% z variabel bebas

% Definisikan Fungsi Bessel Bola Tipe 1
kons = 0.5*pi./z;
jn = ((kons).^(0.5)).*besselj(n+0.5,z);

% Fungsi Riccati Bessel Tipe 1
psi = z.* jn;

% Mencari Nilai Turunan Fungsi Riccati Bessel Tipe 1
jnl = (sqrt(0.5*pi./z)).*besselj(n-0.5,z); 
pso = z.*jnl - n.*jn;




