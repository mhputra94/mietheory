% Fungsi Riccati Bessel Tipe 2
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 21 Februari 2015

function [psi, pso] = ricbesdua(n,z)
% Output
% psi merupakan fungsi Riccati Bessel Tipe 2
% pso merupakan turunan fungsi RIccati Bessel Tipe 2

% Input
% n orde fungsi Riccati Bessel
% z variabel bebas

% Definisikan Fungsi Bessel Bola Tipe 2
kons = 0.5*pi./z;
yn = ((kons).^(0.5)).*bessely(n+0.5,z);

% Fungsi Riccati Bessel Tipe 2
psi = z.* yn;

% Mencari Nilai Turunan Fungsi Riccati Bessel Tipe 2
ynl = (sqrt(0.5*pi./z)).*bessely(n-0.5,z); 
pso = z.*ynl - n.*yn;




