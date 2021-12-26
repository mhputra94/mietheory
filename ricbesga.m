% Fungsi Riccati Bessel Tipe 3
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 21 Februari 2015

function [psit, psot] = ricbesga(n,z)
% Output
% psi merupakan fungsi Riccati Bessel Tipe 3
% pso merupakan turunan fungsi Riccati Bessel Tipe 3

% Input
% n orde fungsi Riccati Bessel
% z variabel bebas

% Definisikan Fungsi Bessel Bola Tipe 3
i = sqrt(-1);
kons = 0.5*pi./z;
jn = ((kons).^(0.5)).*besselj(n+0.5,z);
yn = ((kons).^(0.5)).*bessely(n+0.5,z);
hn1 = jn + i*yn;
% Fungsi Riccati Bessel Tipe 3
psit = z.* hn1;

% Mencari Nilai Turunan Fungsi Riccati Bessel Tipe 3
jnl = (sqrt(0.5*pi./z)).*besselj(n-0.5,z); 
ynl = (sqrt(0.5*pi./z)).*bessely(n-0.5,z);
hnl = jnl + i*ynl;
psot = z.*hnl - n.*hn1;




