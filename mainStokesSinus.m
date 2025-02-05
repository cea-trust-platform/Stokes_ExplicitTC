%{
/****************************************************************************
* Copyright (c) 2024, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% mainStokesSinus.m:
%
% Resolution du probleme de Stokes "Sinus" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(b-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                               sin(npi*x)*{(b+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=sin(npi*x)*sin(npi*y);
%
% Programme principal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES MODIFIABLES
%%%%%%%%%%%%%%%%%%%%%%%%
m0=1; % CHOICE : 1 to 4
m1=4; % CHOICE : m0 to 4
%
global eps nitMAX
eps=1.e-12;  % precision GCP
nitMAX=1000; % nb iterations max GCP
global ordre % CHOICE : 1 or 2
ordre=1;
% Basic TC :    lambda>0.25*(Cdiv)^2, Cdiv<~sqrt(2) dans le carre
% Explicit TC : lambda>0
global nu lambda nmax epsTC
% Remarque : on garde lambda=1 lorsque la pression initiale est assez precise
%nu=1.e+00; lambda=1.e+00; epsTC=eps; nmax=8;
nu=1.e-06; lambda=1.e+01; epsTC=eps; nmax=8;
%
global gama bp
gama=1;    % Longueur d'onde. Choisir un nombre entier.
bp=1;      % coefficient devant la pression. Choisir 1 ou 0.
%%%%%%%%%%%%%%%%%
% NON MODIFIABLES
%%%%%%%%%%%%%%%%%
global Vol=1;
nmesh=m1-m0+1;
meshstep=[0.1,0.05,0.025,0.0125,0.00625];
TabNbpt  =zeros(1,nmesh); %[142,  568, 2212,  8558,  3mi4239];
TabNbtri =zeros(1,nmesh); %[242, 1054, 4262, 16794,  68476];
TabNbedg =zeros(1,nmesh); %[383, 1621, 6473, 25351, 102074];
nEup=7;
Eu0=zeros(nmesh,nEup); Eu1=zeros(nmesh,nEup); Ep0=zeros(nmesh,nEup);
Eud=zeros(nmesh,nEup); 
nit=zeros(nmesh,nEup); tps=zeros(nmesh,nEup);
TabLG=zeros(1,nmesh); TabNC=zeros(1,nmesh);
%
% Image finale
fig=-nmesh+1;
%%%%%%%%%
% CALCULS
%%%%%%%%%
im=0;
for mii=m0:m1
  im=im+1;
  meshname='Square_h';
  filename=sprintf('%s%i',meshname,mii);
  fprintf(' Maillage %s\n',filename);
  global Nbpt CoorNeu CoorNeu2 RefNeu
  global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
  global Nbedg NumEdg CoorMil RefEdg LgEdg2 EdgNorm EdgTri SomOpp
  global mi pas
  global npi=gama*2*pi npi2=npi^2
  mi=mii; pas=meshstep(mii);
  %
  id = tic;
  [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil,...
         RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,LgEdg2,EdgNorm,Aires]=readmeshfiles(filename);
  elapsed_time=toc(id);
  fprintf('Temps de lecture des fichiers de maillage = %7.2e s\n',elapsed_time);
  fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
  %
  Nbpt = size(CoorNeu,1); TabNbpt(im) = Nbpt;
  Nbtri= size(NumTri,1) ; TabNbtri(im)= Nbtri;
  Nbedg= size(NumEdg,1) ; TabNbedg(im)= Nbedg;
  TabLG(1,im)=2*Nbpt+Nbtri;
  TabNC(1,im)=2*Nbedg+Nbtri;
  %
  id = tic;
  [Ku,Mu,Mp,invMp,Bx,By] = MatStokesSV();
  elapsed_time=toc(id);
  fprintf('Temps d assemblage des matrices LG = %7.2e s\n',elapsed_time);
  id = tic;
  [KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC] = MatStokesNC();
  elapsed_time=toc(id);
  fprintf('Temps d assemblage des matrices NC = %7.2e s\n',elapsed_time);
  %
  if (fig<=0)
    fig=fig+1;
  end
  [Eu0(im,:),Eu1(im,:),Eud(im,:),Ep0(im,:),nit(im,:),tps(im,:),fig]=...
  SolveStokesSinus(fig,Ku,Mu,Mp,invMp,Bx,By,KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC);
end % end mi
%%%%%%%%%%%%%
% CONVERGENCE
%%%%%%%%%%%%%
if (nmesh>1)
  % Maillage
  logH =log10(meshstep)';
  dlogH =logH (1:nmesh-1)-logH (2:nmesh);
  invdlogH=dlogH.^(-1); 
  % Nb DoF LG
  logLG=0.5*log10(TabLG)'; 
  dlogLG=logLG(1:nmesh-1)-logLG(2:nmesh);
  invdlogLG=dlogLG.^(-1); 
  % Nb DoF NC
  logNC=0.5*log10(TabNC);
  dlogNC=logNC(1:nmesh-1)-logNC(2:nmesh);
  invdlogNC=dlogNC.^(-1);
  % Erreurs vitesse pression
  logEu0=log10(Eu0); dlogEu0=logEu0(1:nmesh-1,:)-logEu0(2:nmesh,:);
  logEu1=log10(Eu1); dlogEu1=logEu1(1:nmesh-1,:)-logEu1(2:nmesh,:);
  logEud=log10(Eud); dlogEud=logEud(1:nmesh-1,:)-logEud(2:nmesh,:);
  logEp0=log10(Ep0); dlogEp0=logEp0(1:nmesh-1,:)-logEp0(2:nmesh,:);
  %  
  % Convergence en pas du maillage
  vtauEu0=dlogEu0.*invdlogH; tauEu0=sum(vtauEu0,1)/(nmesh-1);
  vtauEu1=dlogEu1.*invdlogH; tauEu1=sum(vtauEu1,1)/(nmesh-1);
  vtauEud=dlogEud.*invdlogH; tauEud=sum(vtauEud,1)/(nmesh-1);
  vtauEp0=dlogEp0.*invdlogH; tauEp0=sum(vtauEp0,1)/(nmesh-1);
  % Convergence en nombre de ddl
  vtauEu0n=-dlogEu0.*invdlogLG; tauEu0n=sum(vtauEu0n,1)/(nmesh-1);
  vtauEu1n=-dlogEu1.*invdlogLG; tauEu1n=sum(vtauEu1n,1)/(nmesh-1);
  vtauEudn=-dlogEud.*invdlogLG; tauEudn=sum(vtauEudn,1)/(nmesh-1);
  vtauEp0n=-dlogEp0.*invdlogLG; tauEp0n=sum(vtauEp0n,1)/(nmesh-1);
  %
  if (fig>1)
    titre=sprintf('Cv P%i-P%i,nu=%7.2e',ordre,ordre-1,nu);
    figure(fig);
    plot(logLG,logEu0(:,2),'-b;Eu0 CR;',...
         logLG,logEu1(:,2),'-r;Eu1 CR;',...
         logLG,logEp0(:,2),'-g;Ep0 CR;',...
         logLG,logEu0(:,3),'--b;Eu0 TS;',...
         logLG,logEu1(:,3),'--r;Eu1 TS;',...
         logLG,logEp0(:,3),'--g;Ep0 TS;',...
         logLG,logEu0(:,4),'-.b;Eu0 NS;',...
         logLG,logEu1(:,4),'-.r;Eu1 NS;',...
         logLG,logEp0(:,4),'-.g;Ep0 NS;');
    xlabel ('log10(Ndof)');
    ylabel ('log10(Erreur)');
  endif % fig>1
  %
  sortie=sprintf('StokesSinus-P%iP%i',ordre,ordre-1);
  sortie=sprintf('%s-nu=%5.0e-lambda=%5.0e-eps=%5.0e-nmax=%i',sortie,nu,lambda,epsTC,nmax);
  %
  tps(:,3)=tps(:,2)+tps(:,3); 
  tps(:,4)=tps(:,2)+tps(:,4);
  titre1='eU0EP    eU1EP    eUDEP    eP0EP    tpsEP    eU0NC    eU1NC    eUDNC    eP0NC    tpsNC    nNC ';
  titre2='eU0TS    eU1TS    eUDTS    eP0TS    tpsTS    eU0NS    eU1NS    eUDNS    eP0NS    tpsNS    nNS h\n';
  entete=sprintf('%s%s',titre1,titre2);
  %
  fid=fopen(sortie,'w');    
  fprintf(fid,entete);
  for m=m0:m1
    fprintf(fid,'%7.2e %7.2e %7.2e %7.2e %7.2e ',     Eu0(m,1),Eu1(m,1),Eud(m,1),Ep0(m,1),tps(m,1));
    fprintf(fid,'%7.2e %7.2e %7.2e %7.2e %7.2e %3i ', Eu0(m,2),Eu1(m,2),Eud(m,2),Ep0(m,2),tps(m,2),nit(m,2));
    fprintf(fid,'%7.2e %7.2e %7.2e %7.2e %7.2e ',     Eu0(m,3),Eu1(m,3),Eud(m,3),Ep0(m,3),tps(m,3));
    fprintf(fid,'%7.2e %7.2e %7.2e %7.2e %7.2e %3i ', Eu0(m,4),Eu1(m,4),Eud(m,4),Ep0(m,4),tps(m,4),nit(m,4));
    fprintf(fid,'%7.2e \n',meshstep(1,m));
  end
  %
  toprint='%7.2e %7.2e %7.2e %7.2e';
  space1='          '; 
  space2='              ';
  last=' \n';
  toprintfull=sprintf('%s%s%s%s%s%s%s%s%s%s',toprint,space1,toprint,space2,toprint,space1,toprint,last);
  
  for j=1:nmesh-1
    fprintf(fid,toprintfull,...
             vtauEu0(j,1),vtauEu1(j,1),vtauEud(j,1),vtauEp0(j,1),...
             vtauEu0(j,2),vtauEu1(j,2),vtauEud(j,2),vtauEp0(j,2),...
             vtauEu0(j,3),vtauEu1(j,3),vtauEud(j,3),vtauEp0(j,3),...
             vtauEu0(j,4),vtauEu1(j,4),vtauEud(j,4),vtauEp0(j,4));
  end
  %
  fprintf(fid,toprintfull,...
             tauEu0(1),tauEu1(1),tauEud(1),tauEp0(1), ...
             tauEu0(2),tauEu1(2),tauEud(2),tauEp0(2), ...
             tauEu0(3),tauEu1(3),tauEud(3),tauEp0(3), ...
             tauEu0(4),tauEu1(4),tauEud(4),tauEp0(4));
  fclose(fid);
  %
  %%%%%%%%%%%%%%
  %
  sortieRT=sprintf('StokesSinusRT-P%iP%i',ordre,ordre-1);
  sortieRT=sprintf('%s-nu=%5.0e-lambda=%5.0e-eps=%5.0e-nmax=%i',sortieRT,nu,lambda,epsTC,nmax);
  %
  tps(:,6)=tps(:,5)+tps(:,6); 
  tps(:,7)=tps(:,5)+tps(:,7);
  titre1RT='eU0RT    eU1RT    eUDRT    eP0RT    tpsRT    nRT ';
  titre2RT='eU0TS    eU1TS    eUDTS    eP0TS    tpsTS    eU0NS    eU1NS    eUDNS    eP0NS    tpsNS    nNS h\n';
  enteteRT=sprintf('%s%s',titre1RT,titre2RT);
  %
  fidRT=fopen(sortieRT,'w');    
  fprintf(fidRT,enteteRT);
  for m=m0:m1
    fprintf(fidRT,'%7.2e %7.2e %7.2e %7.2e %7.2e %3i ', Eu0(m,5),Eu1(m,5),Eud(m,5),Ep0(m,5),tps(m,5),nit(m,5));
    fprintf(fidRT,'%7.2e %7.2e %7.2e %7.2e %7.2e ',     Eu0(m,6),Eu1(m,6),Eud(m,6),Ep0(m,6),tps(m,6));
    fprintf(fidRT,'%7.2e %7.2e %7.2e %7.2e %7.2e %3i ', Eu0(m,7),Eu1(m,7),Eud(m,7),Ep0(m,7),tps(m,7),nit(m,7));
    fprintf(fidRT,'%7.2e \n',meshstep(1,m));
  end
  %
  toprint='%7.2e %7.2e %7.2e %7.2e';
  space1='          '; 
  space2='              ';
  last=' \n';
  toprintfull=sprintf('%s%s%s%s%s%s%s%s%s%s',toprint,space1,toprint,space2,toprint,space1,last);
  
  for j=1:nmesh-1
    fprintf(fidRT,toprintfull,...
             vtauEu0(j,5),vtauEu1(j,5),vtauEud(j,5),vtauEp0(j,5),...
             vtauEu0(j,6),vtauEu1(j,6),vtauEud(j,6),vtauEp0(j,6),...
             vtauEu0(j,7),vtauEu1(j,7),vtauEud(j,7),vtauEp0(j,7));
  end
  %
  fprintf(fidRT,toprintfull,...
             tauEu0(5),tauEu1(5),tauEud(5),tauEp0(5), ...
             tauEu0(6),tauEu1(6),tauEud(6),tauEp0(6), ...
             tauEu0(7),tauEu1(7),tauEud(7),tauEp0(7));
  fclose(fidRT);
  %
 end
