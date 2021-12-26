% Program Mie Theory untuk Partikel Bola Tunggal
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 21 Februari 2015

%clear all
%clc

mat = 'silverpalik.dat'; % Input Indeks Material Nanopartikel
next = 1.3334; % Indeks Bias Eksternal atau Lingkungan
erei = 280:1:800; % Input Panjang Gelombang
r = 50; % jari - jari
l = 3; % pole mode (dipole, quadrupole or multipole)

[ene, n, k] = textread(mat, '%f %f %f','commentstyle','matlab');

units

enei = eV2nm ./ ene;
ni = spline( enei, n, erei );
ki = spline( enei, k, erei );
i = sqrt(-1);
nmat = ni + i*ki; % Indeks Bias Material

% Memulai menghitung Mie Theory
x = (2*pi*next*r)./erei;
m = nmat./next;
k = (2*pi*next)./erei;

o = length(erei);
sca = zeros(o,1);
ext = zeros(o,1);

% Menghitung Output
for j = 1:o
   %for t = 1:l
       % Menentukkan Koefisien Hamburan
       
       [psi,psit] = ricbestu(l,(m(j)*x(j)));
       [pso,psot] = ricbestu(l,x(j));
       [zet,zett] = ricbesga(l,x(j));
       anom = (m(j)*psi*psot) - (pso*psit);
       adenom = (m(j)*psi*zett) - (zet*psit);
       an = anom/adenom;
       bnom = (psi*psot) - (m(j)*pso*psit);
       bdenom = (psi*zett) - (m(j)*zet*psit);
       bn = bnom/bdenom;
       
       insca = (2*pi/(k(j))^2)*(2*l+1)*((abs(an))^2 + (abs(bn))^2);
       sca(j) = sca(j) + insca;
       inext = (2*pi/(k(j))^2)*(2*l+1)*((real(an + bn)));
       ext(j) = ext(j) + inext;
   %end  
end

% Menghitung Serapan
abso = ext - sca;

% Menghitung nilai effisiensi
luas = pi*(r^2);
scaef = sca/luas;
extef = ext/luas;
absef = abso/luas;

% Plotting Grafik 
plot(erei',scaef,'b')
hold on
plot(erei',absef,'g')
plot(erei',extef,'r')
legend('Scattering','Absorption','Extinction')
grid on;