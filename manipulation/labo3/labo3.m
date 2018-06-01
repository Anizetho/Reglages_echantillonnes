% Pour Matlab 2014a
close all
clear;
%description du systeme continu
nH=1;
dH=conv([1 0],[1 0.5]);
H=tf(nH,dH);% Construction de l'objet du systeme continu

%Cahier des charges du systeme en BF
T95=1;
T=T95/3;
pcr=-1/T;
amo=0.69;
W0=-pcr/amo;
Wd=sqrt(W0^2-pcr^2);
pcr1=pcr+1i*Wd;
pcr2=pcr-1i*Wd;
dfc=poly([pcr1;pcr2]);
nfc=dfc(3);
Fc=tf(nfc,dfc);

% poles discrets imposes en Bf
h=0.1;%periode d'echantillonnage du regulateur discret
pdr1=exp(pcr1*h);%poles discrets du modele de reference
pdr2=exp(pcr2*h);
Fd=tf([1 0],poly([pdr1;pdr2]),h);
Fd=Fd/dcgain(Fd);
[nFd,dFd]=tfdata(Fd, 'v');

Hd=c2d(H,h);%Equivalent echantillonne bloque du systeme continu a la periode h
Hd=zpk(Hd);%ecriture du system sous la forme poles/zeros

%sortir la partie inversible du systeme discret
[zd,pd,kd]=zpkdata(Hd,'v');%Extraire les poles et zeroes du systeme
pdd=pd(abs(pd-1)>eps*100);%On sort tous les poles qu'on veut simplifier, sauf l'integrateur

%calcul du regulateur discret
zrd=pdd;%les zeros du regulateur discret sont les poles
prd=zd;%les poles du regulateurs sont les zeroes du systemes
Rdsi=zpk(zrd,prd,1/kd,h);%partie inversible du regulateur, on converse l'integrateur du systeme
id=tf(1,[1 -1],h);% on construit l'intregrateur discret de gain d'Evans unite
Rd=series(Rdsi,id);%on ajoute l'integrateur dans le regulateur


ke=(2-pdr1-pdr2);%calcul du gain d'Evans de la boucle ouverte
zra=(1-pdr1*pdr2)/ke;%calcul du zero additionnel
Hza=zpk(zra,[],ke,h);%zero additionnel et gain
Rdt=Rd*Hza;
Bo=minreal(series(Rdt,Hd));% La boucle ouverte est un double integrateur de gain d'Evans unite
%zrt=[zrd zra];
[nzd,dzd]=tfdata(Rdt,'v');
F=feedback(Bo,1);
[nf,df]=tfdata(F,'v');

%ouvrir le modele simulink
slabo3
sim('slabo3',4)
t1=S1.time;
D1=[S1.signals.values];
t2=S2.time;
D2=[S2.signals.values];
subplot(211);
plot(t1,D1(:,1));grid;
A=axis;
px=A(1)+(A(2)-A(1))*0.4;
py=A(3)+(A(4)-A(3))*0.3;
dpy=(A(4)-A(3))/10;
nFtxt=sprintf('%3.7g   ',nzd);
dFtxt=sprintf('%3.7g   ',dzd);
text(px,py,['nF(Z)= ' nFtxt]);
text(px,py-dpy,['dF(Z)= ' dFtxt]);
subplot(212);
plot(t2,D2(:,2),t2,D2(:,3),t1,D1(:,2));
grid
A=axis;
px=A(1)+(A(2)-A(1))*0.4;
py=A(3)+(A(4)-A(3))*0.5;
dpy=(A(4)-A(3))/10;
nzdtxt=sprintf('%3.7g  ',nzd);
dzdtxt=sprintf('%3.7g  ',dzd);
pstxt=sprintf('%3.7g  ',zra);
prdtxt=sprintf('%3.7g ',prd);
htxt=sprintf('%3.7g  ',h);
text(px,py,['NRd(Z)= ' nzdtxt]);
text(px,py-dpy,['DRd(Z)= ' dzdtxt]);
text(px,py-2*dpy,['pole ajoute en  ' pstxt]);
text(px,py-3*dpy,['pole simplifie en  ' prdtxt]);
text(px,py-4*dpy,['Periode d''echantillonnage  ' htxt]);
clc
