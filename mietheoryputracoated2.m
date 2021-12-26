% Program Mie Theory untuk Partikel Bola Tunggal Coated
% Bagian Dalam Non Logam
% Bagian Luar Logam
% Dibuat oleh Miftahussurur Hamidi Putra
% Tanggal 22 Februari 2015

clear all
clc
% Input Material
matlu = 'silverpalik.dat'; % Material Luar
nmatdal = 1.7; % Indeks Bias Material Dalam
nling = 1.3334; % Indeks Bias Lingkungan
erei = 300:1:700; % Input Panjang Gelombang
rdal = 30; % jari - jari dalam
rlua = 31; % jari - jari luar
l = 50; % banyak iterasi
[ene, n, k] = textread(matlu, '%f %f %f','commentstyle','matlab');

units

enei = eV2nm ./ ene;
ni = spline( enei, n, erei );
ki = spline( enei, k, erei );
i = sqrt(-1);
nmat = ni + i*ki; % Indeks Bias Material Luar

% Menghitung Mie Theory

x = (2*pi*nling*rdal)./erei;
y = (2*pi*nling*rlua)./erei;
m1 = nmatdal/nling;
m2 = nmat./nling;
k = (2*pi*nling)./erei;

o = length(erei);
sca = zeros(o,1);
ext = zeros(o,1);

% Menghitung Output
for j = 1:o
   for t = 1:l
       % Menentukkan Koefisien Hamburan
       
       % Menghitung Fungsi Riccati Bessel
      [psi, psit] = ricbestu(t,m2(j)*x(j));
      [pso, psot] = ricbestu(t,m1*x(j));
      [pse, pset] = ricbestu(t,y(j));
      [psa, psat] = ricbestu(t,m2(j)*y(j));
      [cha, chat] = ricbesdua(t,m2(j)*x(j));
      [che, chet] = ricbesdua(t,m2(j)*y(j));
      [zet, zett] = ricbesga(t,y(j));
      chi = -1*cha; chit = -1*chat;
      cho = -1*che; chot = -1*chet;
      
      % Koefisien An => Ane dan Bn => Bne
      Anum = m2(j)*psi*psot - m1*psit*pso;
      Adenum = m2(j)*chi*psot - m1*chit*pso;
      Ane = Anum/Adenum;
      
      Bnum = m2(j)*pso*psit - m1*psi*psot;
      Bdenum = m2(j)*chit*pso - m1*psot*chi;
      Bne = Bnum/Bdenum;
      
      % Koefisien an dan bn
      anom = pse*(psat - (Ane*chot)) - ((m2(j)*pset)*(psa - (Ane*cho)));
      adenom = zet*(psat - (Ane*chot)) - ((m2(j)*zett)*(psa - (Ane*cho)));
      an = anom/adenom;
      
      bnom = ((m2(j)*pse)*(psat - (Bne*chot))) - pset*(psa - (Bne*cho));
      bdenom = ((m2(j)*zet)*(psat - (Bne*chot))) - zett*(psa - (Bne*cho));
      bn = bnom/bdenom;
      
      % Menghitung Hamburan dan Ekstinsi
      insca = (2*pi/(k(j))^2)*(2*t+1)*((abs(an))^2 + (abs(bn))^2);
       sca(j) = sca(j) + insca;
       inext = (2*pi/(k(j))^2)*(2*t+1)*((real(an + bn)));
       ext(j) = ext(j) + inext;
             
   end  
end

% Menghitung Serapan
abso = ext - sca;
plot(erei',abso)
grid on



