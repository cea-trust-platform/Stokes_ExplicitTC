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
% SolveStokesSinus.m:
%
% Resolution du probleme de Stokes "Sinus" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(1-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                               sin(npi*x)*{(1+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ;
%         (cos(npi*y)-1)*sin(npi*x)];
% p(x,y)=sin(npi*x)*sin(npi*y);
%
% Construction et resolution du systeme lineaire, calcul des erreurs
%
% SYNOPSIS [Eu0,Eu1,Ep0,nit,tps,fig]=SolveStokesSinusSV(fig,Ku,Mu,Sp,Mp,invMp,Bx,By,...
%                                  KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordonnees (x, y) des noeuds P2
%        - RefNeu(Nbpt,1) : reference des sommets
%        - CoorBary(Nbtri,3) :coordonnees (x, y) des barycentres des triangles
%        - CoorMil(Nbedg,2)   : Coordonnees des milieux d'aretes
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets) 
%        - NumTri2(4*Nbtri,3) : liste de triangles du maillage P2
%                   (3 numeros de sommets)
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - LgEdg2(Nbedg,1) : longueurs des aretes au carre
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1) : aires des triangles
%        - fig : numero de figure si visualisation
%        - mi  : numero du maillage
%        - ordre    : ordre d'approximation
% OUTPUT - Eu0 : Erreur L2 normalisee de la vitesse, calcul decompose
%        - Eu1 : Erreur H1 normaliseede la vitesse, calcul decompose
%        - Ep0 : Erreur L2 normalisee de la pression, calcul decompose
%        - fig = numero de la derniere figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Eu0,Eu1,Eud,Ep0,nit,tps,fig]=SolveStokesSinus(fig,Ku,Mu,Mp,invMp,Bx,By,KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC)
%
global Nbpt CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
global Nbedg CoorMil RefEdg LgEdg2 EdgNorm EdgTri
global ordre mi nu nitMAX stab Vol bp nmax lambda
global npi npi2
%
nu2=nu*nu;
normPex2=0.25*bp^2;
normGUex2=2*npi2;
normUP2=normGUex2+normPex2/nu2;
nu2normUP2=nu2*normGUex2+normPex2;
%
nEup=7;
Eu0=zeros(1,nEup); Eu1=zeros(1,nEup); Ep0=zeros(1,nEup); nit=ones(1,nEup); tps=zeros(1,nEup);
Eud=zeros(1,nEup);
%
oU=ordre; oP=oU-1;
%
Np=size(Mp,1); UnP=ones(Np,1);
Vol=UnP'*Mp*UnP;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traitement des CL de Dirichlet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DofU_S=find(RefNeu==0); DofUB_S=find(RefNeu~=0);
DofU_A=find(RefEdg==0); DofUB_A=find(RefEdg~=0);
if (ordre==1)
  DofU=DofU_S; DofUB=DofUB_S;
  DofUNC=DofU_A; DofUNCB=DofUB_A; 
  Nu=Nbpt; NuNC=Nbedg;
end
if (ordre==2)
  NP2=Nbpt+Nbedg;
  DofU=[DofU_S;Nbpt+DofU_A];
  DofU_T=linspace(1,Nbtri,Nbtri)';
  DofUNC=[DofU_S;Nbpt+DofU_A;NP2+DofU_T];
  Nu=Nbpt+Nbedg; NuNC=Nu+Nbtri;
  DofUB_A=find(RefEdg~=0);
  DofUB  =[DofUB_S;Nbpt+DofUB_A];
  DofUNCB=DofUB;
end
%
NuB=size(DofUB,1); NuNCB=size(DofUNCB,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximations de la solutions exacte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PexP,moyPexP,PexLG]=PexStokesSinus(Mu,Mp);
%
npiXYS=npi*CoorNeu;
cosXYS=cos(npiXYS); sinXYS=sin(npiXYS);
UxS=((1-cosXYS(:,1)).*sinXYS(:,2));
UyS=((cosXYS(:,2)-1).*sinXYS(:,1));
%
npiXYA=npi*CoorMil;
cosXYA=cos(npiXYA); sinXYA=sin(npiXYA); 
UxA=((1-cosXYA(:,1)).*sinXYA(:,2));
UyA=((cosXYA(:,2)-1).*sinXYA(:,1));
%
if (ordre==1)
  UxI=UxS;   UyI=UyS;
  UxINC=UxA; UyINC=UyA;
  npiXYT=npi*CoorBary;
  PexI=sin(npiXYT(:,1)).*sin(npiXYT(:,2));
end
%
if (ordre==2)
  UxI=[UxS;UxA]; UyI=[UyS;UyA];
  UxINC=[UxI;zeros(Nbtri,1)]; UyINC=[UyI;zeros(Nbtri,1)];
  PexS=sin(npiXYS(:,1)).*sin(npiXYS(:,2));
  PexI=zeros(3*Nbtri,1);
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    fin=3*t;
    PexI(fin-2:fin,1)=PexS(IGLO,1);
  endfor
endif
%
PexI=PexI-UnP'*Mp*PexI/Vol;
%
UxI(DofUB,1)=zeros(NuB,1); UyI(DofUB,1)=zeros(NuB,1);
UxINC(DofUNCB,1)=zeros(NuNCB,1); UyINC(DofUNCB,1)=zeros(NuNCB,1);
%
UxhB=zeros(Nu,1); UyhB=zeros(Nu,1);
UxhNCB=zeros(NuNC,1); UyhNCB=zeros(NuNC,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREMIER CALCUL : Pi_{dg}(p) (presque) exacte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 1: P%i-P%idg, mesh %i\n',oU,oP,mi);
[RHS_Ux,RHS_Uy]=RHSStokesSinusSV();
lambda0=lambda;
lambda=1;
id = tic;
[Uxh,Uyh,Ph,Eud(1)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PexP,Vol);
tps(1)=toc(id);
%
dUxhI=UxI-Uxh; dUyhI=UyI-Uyh; dPhP=PexP-Ph; 
Eu0(1)=sqrt((dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI)/normUP2);
Eu1(1)=sqrt((dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI)/normUP2);
Eud(1)/=sqrt(normUP2);
Ep0(1)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
%
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(1));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(1));
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1(1));
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',Eud(1));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(1));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%
% SECOND CALCUL : EF NC
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 2: P%inc-P%idg, mesh %i\n',oU,oP,mi);
[RHS_UxNC,RHS_UyNC]=RHSStokesSinusNC();
id = tic;
[UxhNC,UyhNC,PhNC,nit(2)]=UzawaNC(KuNC,BxNC,ByNC,MpNC,invMpNC,UnP,RHS_UxNC,RHS_UyNC,UxhNCB,UyhNCB,DofUNC,DofUNCB,Vol);
tps(2)=toc(id);
%
dUxhINC=UxINC-UxhNC; dUyhINC=UyINC-UyhNC; dPhP=PexP-PhNC;
Eu0(2)=sqrt((dUxhINC'*MuNC*dUxhINC+dUyhINC'*MuNC*dUyhINC)/normUP2);
Eu1(2)=sqrt((dUxhINC'*KuNC*dUxhINC+dUyhINC'*KuNC*dUyhINC)/normUP2);
Eud(2)=ErrorDivUNC(UxhNC,UyhNC)/sqrt(normUP2);
Ep0(2)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
%
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(2),nit(2))
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(2));
fprintf('         ||Uex-Uh||_h/||(U,P)||_X = %7.2e\n',Eu1(2));
fprintf('             ||Uh||_J/||(U,P)||_X = %7.2e\n',Eud(2));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(2));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TROISIEME CALCUL : PRESSION NC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 3:  P%i-P%idg, mesh %i\n',oU,oP,mi);
lambda=lambda0;
id = tic;
[Uxh,Uyh,Ph,Eud(3)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PhNC,Vol);
tps(3)=toc(id);
%
dUxhI=UxI-Uxh; dUyhI=UyI-Uyh; dPhP=PexP-Ph; 
Eu0(3)=sqrt((dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI)/normUP2);
Eu1(3)=sqrt((dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI)/normUP2);
Eud(3)/=sqrt(normUP2);
Ep0(3)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
%
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(3));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(3));
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1(3));
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',Eud(3));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(3));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUATRIEME CALCUL : PRESSION NC + iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 4: P%i-P%idg, mesh %i\n',oU,oP,mi);
id = tic;
dPh0=Ph-PhNC;
[Uxh,Uyh,Ph,Eud(4),nit(4)]=SolveStokesExplicitTC_it(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Uxh,Uyh,DofU,DofUB,Ph,dPh0,Vol);
tps(4)=toc(id);
%
dUxhI=UxI-Uxh; dUyhI=UyI-Uyh; dPhP=PexP-Ph;  
Eu0(4)=sqrt((dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI)/normUP2);
Eu1(4)=sqrt((dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI)/normUP2);
Eud(4)/=sqrt(normUP2);
Ep0(4)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
%
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(4),nit(4));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(4));
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1(4));
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',Eud(4));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(4));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CINQUIEME CALCUL : EF NC + PROJECTION SUR LES RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RHS_UxNCRT,RHS_UyNCRT]=RHSStokesSinusNCRT();
id = tic;
[UxhNCRT,UyhNCRT,PhNCRT,nit(5)]=UzawaNC(KuNC,BxNC,ByNC,MpNC,invMpNC,UnP,RHS_UxNCRT,RHS_UyNCRT,UxhNCB,UyhNCB,DofUNC,DofUNCB,Vol);
tps(5)=toc(id);
%
dUxhINCRT=UxINC-UxhNCRT; dUyhINCRT=UyINC-UyhNCRT; dPhP=PexP-PhNCRT;
Eu0(5)=sqrt((dUxhINCRT'*MuNC*dUxhINCRT+dUyhINCRT'*MuNC*dUyhINCRT)/normUP2);
Eu1(5)=sqrt((dUxhINCRT'*KuNC*dUxhINCRT+dUyhINCRT'*KuNC*dUyhINCRT)/normUP2);
Eud(5)=ErrorDivUNC(UxhNCRT,UyhNCRT)/sqrt(normUP2);
Ep0(5)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
%
fprintf('Calcul 5: P%inc-P%idg+RT, mesh %i\n',oU,oP,mi);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(5),nit(5));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(5));
fprintf('         ||Uex-Uh||_h/||(U,P)||_X = %7.2e\n',Eu1(5));
fprintf('             ||Uh||_J/||(U,P)||_X = %7.2e\n',Eud(5));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(5));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection UNCRT sur EF Pk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Projection du calcul 5, mesh %i\n',mi);
id = tic;
[Uxh_6LG,Uyh_6LG]=NCtoLG(UxhNCRT,UyhNCRT,Mu);
tps_LG=toc(id);
dUxI=Uxh_6LG-UxI; dUyI=Uyh_6LG-UyI;
Eu0I=sqrt((dUxI'*Mu*dUxI+dUyI'*Mu*dUyI)/normUP2);
Eu1I=sqrt((dUxI'*Ku*dUxI+dUyI'*Ku*dUyI)/normUP2);
%
divUI=Bx*Uxh_6LG+By*Uyh_6LG; EudI=sqrt(divUI'*invMp*divUI/normUP2);
%
fprintf('-------\n');
fprintf('Temps de projection (P%inc vers P%i) = %7.2e s\n',oU,oU,tps_LG);
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0I);
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1I);
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',EudI);
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIXIEME CALCUL : PRESSION NC-RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 6:  P%i-P%idg, mesh %i\n',oU,oP,mi);
lambda0=lambda;
lambda=1;
id = tic;
[Uxh_6,Uyh_6,Ph_6,Eud(6)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PhNCRT,Vol);
tps(6)=toc(id);
%
dUxhI=UxI-Uxh_6; dUyhI=UyI-Uyh_6; dPhP=PexP-Ph_6;  
Eu0(6)=sqrt((dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI)/normUP2);
Eu1(6)=sqrt((dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI)/normUP2);
Eud(6)/=sqrt(normUP2);
Ep0(6)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(6));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(6));
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1(6));
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',Eud(6));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(6));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUATRIEME CALCUL : PRESSION NC + iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 7: P%i-P%idg, mesh %i\n',oU,oP,mi);
id = tic;
dPh0=Ph_6-PhNCRT;
[Uxh_7,Uyh_7,Ph_7,Eud(7),nit(7)]=SolveStokesExplicitTC_it(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Uxh_6,Uyh_6,DofU,DofUB,Ph_6,dPh0,Vol);
tps(7)=toc(id);
%
dUxhI=UxI-Uxh_7; dUyhI=UyI-Uyh_7; dPhP=PexP-Ph_7;  
Eu0(7)=sqrt((dUxhI'*Mu*dUxhI+dUyhI'*Mu*dUyhI)/normUP2);
Eu1(7)=sqrt((dUxhI'*Ku*dUxhI+dUyhI'*Ku*dUyhI)/normUP2);
Eud(7)/=sqrt(normUP2);
Ep0(7)=sqrt(dPhP'*Mp*dPhP/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(7),nit(7));
fprintf('         ||Uex-Uh||_0/||(U,P)||_X = %7.2e\n',Eu0(7));
fprintf('           |Uex-Uh|_1/||(U,P)||_X = %7.2e\n',Eu1(7));
fprintf('        ||div(Uh)||_0/||(U,P)||_X = %7.2e\n',Eud(7));
fprintf('(1/nu)||P(Pex)-Ph||_0/||(U,P)||_X = %7.2e\n',Ep0(7));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
lambda=lambda0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fig>0)
  %
  oPT=oP*ones(Nbtri,1);
  %
  PhLG=DGPktoLG(Mu,Ph,oPT,oU);
  PhNCLG=DGPktoLG(Mu,PhNC,oPT,oU);
  %
  if (ordre==1)
    NT=NumTri; CN=CoorNeu;
  end
  if (ordre==2)
    NT=NumTri2; CN=CoorNeu2;
  end
  %
  texX=sprintf('LG P%i, mesh%i, Ux',oU,mi);
  tSVX=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%i, mesh%i, Uxh',oU,oP,nu,lambda,mi);
  %
  texY=sprintf('LG P%i, mesh%i, Uy',oU,mi);
  tSVY=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%i, mesh%i, Uyh',oU,oP,nu,lambda,mi);
  %
  texP=sprintf('LG P%i, mesh%i, Pr',oU,mi);
  tSVP=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%i, mesh%i, Prh',oU,oP,nu,lambda,mi);
  %
  tNCP=sprintf('NC P%i-P%i, nu=%7.2e, lambda=%i, mesh%i, Prh',oU,oP,nu,lambda,mi);
  %
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),UxI);
  view(2);
  shading interp
  title(texX)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Uxh);
  view(2);
  shading interp
  title(tSVX)
  colorbar;
  %
  fig+=1;
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),UyI);
  view(2);
  shading interp
  title(texY)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Uyh);
  view(2);
  shading interp
  title(tSVY)
  colorbar;
  fig=fig+1;
  %
  figure(fig)
  colormap ("jet");
  subplot(1,2,1)
  trisurf(NT,CN(:,1),CN(:,2),PexLG);
  view(2);
  shading interp
  title(texP)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),PhLG);
  view(2);
  shading interp
  title(tSVP)
  colorbar;
  fig=fig+1;%
  figure(fig)
  colormap ("jet");
  subplot(1,2,1)
  trisurf(NT,CN(:,1),CN(:,2),PexLG);
  view(2);
  shading interp
  title(texP)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),PhNCLG);
  view(2);
  shading interp
  title(tNCP)
  colorbar;
  fig=fig+1;
end
