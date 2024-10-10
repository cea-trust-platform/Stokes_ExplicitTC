%{
/****************************************************************************
* Copyright (c) 2023, CEA
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
% SolveStokesXY.m:
%
% Resolution du probleme de Stokes "U=(-Y,X)" dans un carre [0,1]*[0,1]
%   avec les elements finis P1-P0 ou P2-P1dg
%   -nu*Delta U + grad p=2*[(X-X_0,Y-Y_0]
%    div U = 0;
%
% U(x,y)=(-Y,X);
% p(x,y)=X^3+Y^3-1/6;
%
% Construction et resolution du systeme lineaire, calcul des erreurs
%
% SYNOPSIS [Eu0,Eu1,Ep0,nit,tps,fig]=SolveStokesXY(fig,Ku,Mu,Mp,invMp,Bx,By,...
%                                            KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC)
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

function [Eu0,Eu1,Eud,Ep0,nit,tps,fig]=SolveStokesXY(fig,Ku,Mu,Mp,invMp,Bx,By,...
                                            KuNC,MuNC,MpNC,invMpNC,BxNC,ByNC)
%
global Nbpt CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
global Nbedg CoorMil RefEdg LgEdg2 EdgNorm EdgTri
global ordre mi nu Vol lambda
%
nEup=7;
%
Eu0=zeros(1,nEup); Eu1=zeros(1,nEup); Ep0=zeros(1,nEup); nit=ones(1,nEup); tps=zeros(1,nEup);
Eud=zeros(1,nEup);
%
oU=ordre; oP=oU-1;
%
%%%%%%%%%%%%%%%%%
% Solution exacte
%%%%%%%%%%%%%%%%%
Xnod=CoorNeu(:,1); Ynod=CoorNeu(:,2);
Xedg=CoorMil(:,1); Yedg=CoorMil(:,2);
%
if (ordre==1)
  Uxex  =-Ynod; Uyex  =Xnod;
  UxexNC=-Yedg; UyexNC=Xedg;
end
%
if (ordre==2)
  Uxex=-[Ynod;Yedg]; Uyex=[Xnod;Xedg];
  UxexNC=[Uxex;zeros(Nbtri,1)];
  UyexNC=[Uyex;zeros(Nbtri,1)];
end
%
[PexDG,moyPexDG,PexLG,Pex2h]=PexStokesXY(Mu,Mp);
%
 Uex2h=Uxex'*Mu*Uxex+Uyex'*Mu*Uyex;
GUex2h=Uxex'*Ku*Uxex+Uyex'*Ku*Uyex;
%
 Uex2hNC=UxexNC'*MuNC*UxexNC+UyexNC'*MuNC*UyexNC;
GUex2hNC=UxexNC'*KuNC*UxexNC+UyexNC'*KuNC*UyexNC;
% NORME ET QUALITE APPROXIMATION
Eud_exNC=ErrorDivUNC(UxexNC,UyexNC);
div=Bx*Uxex+By*Uyex;
Eud_ex=sqrt(div'*invMp*div);
%fprintf('     ||UexLG||_0 = %7.2e, ||UexNC||_0 = %7.2e\n',sqrt(Uex2h),sqrt(Uex2hNC));
%fprintf('       |UexLG|_1 = %7.2e, ||UexNC||_h = %7.2e\n',sqrt(GUex2h),sqrt(GUex2hNC));
%fprintf('||div(UexLG)||_0 = %7.2e, ||UexNC||_J = %7.2e\n',Eud_ex,Eud_exNC);
%fprintf('       ||Pex||_0 = %7.2e\n',sqrt(Pex2h));
%
Np=size(Mp,1); UnP=ones(Np,1);
Vol=UnP'*Mp*UnP;
Uex2=2/3; GUex2=2;
Pex2=9/56;
nu2=nu*nu; normUP2=GUex2+Pex2/nu2; nu2normUP2=nu2*2+Pex2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traitement des CL de Dirichlet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DofU_S=find(RefNeu==0); DofUB_S=find(RefNeu~=0);
DofU_A=find(RefEdg==0); DofUB_A=find(RefEdg~=0);
if (ordre==1)
  DofU=DofU_S; DofUB=DofUB_S;
  DofUNC=DofU_A; DofUNCB=DofUB_A;
  Nu=Nbpt;
end
if (ordre==2)
  NP2=Nbpt+Nbedg;
  DofU=[DofU_S;Nbpt+DofU_A];
  DofU_T=linspace(1,Nbtri,Nbtri)';
  DofUNC=[DofU_S;Nbpt+DofU_A;NP2+DofU_T];
  Nu=Nbpt+Nbedg;
  DofUB_A=find(RefEdg~=0);
  DofUB  =[DofUB_S;Nbpt+DofUB_A];
  DofUNCB=DofUB;
end
UxhB=zeros(Nu,1); UxhB(DofUB,1)=Uxex(DofUB,1);
UyhB=zeros(Nu,1); UyhB(DofUB,1)=Uyex(DofUB,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREMIER CALCUL : PRESSION EXACTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 1: P%i-P%idg, mesh %i\n',oU,oP,mi);
lambda0=lambda;
lambda=1;
[RHS_Ux,RHS_Uy]=RHSStokesXYSV();
id = tic;
[Uxh,Uyh,Ph,Eud(1)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PexDG,Vol);
tps(1)=toc(id);
%
dP=Ph-PexDG; dUx=Uxh-Uxex; dUy=Uyh-Uyex;
Eu0(1)=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/normUP2);
Eu1(1)=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/normUP2);
Eud(1)/=sqrt(normUP2);
Ep0(1)=sqrt((dP'*Mp*dP)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(1));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,Pex)||_X = %7.2e\n',Eu0(1));
fprintf('      ||Ph(Uex)-U_h||_1/||(Uex,Pex)||_X = %7.2e\n',Eu1(1));
fprintf('         ||div(U_h)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud(1));
fprintf('         ||div(Uex)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud_ex);
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,Pex)||_X = %7.2e\n',Ep0(1));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%
% SECOND CALCUL : EF NC
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 2: P%i-P%idg, mesh %i\n',oU,oP,mi);
%divNC=BxNC*UxexNC+ByNC*UyexNC; Eud_ex=sqrt(divNC'*Mp*divNC/normUP2);
Eud_ex=ErrorDivUNC(UxexNC,UyexNC)/sqrt(normUP2);
%
[RHS_UxNC,RHS_UyNC]=RHSStokesXYNC();
%
NuNC=size(RHS_UxNC,1);
UxhNCB=zeros(NuNC,1); UxhNCB(DofUNCB,1)=UxexNC(DofUNCB,1);
UyhNCB=zeros(NuNC,1); UyhNCB(DofUNCB,1)=UyexNC(DofUNCB,1);
%
##DiffX=nu*KuNC(DofUNC,DofUNC)*UxexNC(DofUNC,1)+nu*KuNC(DofUNC,DofUNCB)*UxhNCB(DofUNCB,1)+BxNC'(DofUNC,:)*PexNC-RHS_UxNC(DofUNC,1);
##norm(DiffX,2)
##DiffY=nu*KuNC(DofUNC,DofUNC)*UyexNC(DofUNC,1)+nu*KuNC(DofUNC,DofUNCB)*UyhNCB(DofUNCB,1)+ByNC'(DofUNC,:)*PexNC-RHS_UyNC(DofUNC,1);
##norm(DiffY,2)
%
id = tic;
[UxhNC,UyhNC,PhNC,nit(2)]=UzawaNC(KuNC,BxNC,ByNC,MpNC,invMpNC,UnP,RHS_UxNC,RHS_UyNC,UxhNCB,UyhNCB,DofUNC,DofUNCB,Vol);
tps(2)=toc(id);
dUxNC=UxhNC-UxexNC; dUyNC=UyhNC-UyexNC; dPNC=PhNC-PexDG;
Eu0(2)=sqrt((dUxNC'*MuNC*dUxNC+dUyNC'*MuNC*dUyNC)/normUP2);
Eu1(2)=sqrt((dUxNC'*KuNC*dUxNC+dUyNC'*KuNC*dUyNC)/normUP2);
Eud(2)=ErrorDivUNC(UxhNC,UyhNC)/sqrt(normUP2);
Ep0(2)=sqrt((dPNC'*MpNC*dPNC)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(2),nit(2));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,Pex)||_X = %7.2e\n',Eu0(2)); 
fprintf('      ||Ph(Uex)-U_h||_h/||(Uex,Pex)||_X = %7.2e\n',Eu1(2));
fprintf('              ||U_h||_J/||(Uex,Pex)||_X = %7.2e\n',Eud(2));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,Pex)||_X = %7.2e\n',Ep0(2));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TROISIEME CALCUL : PRESSION NC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 3: P%i-P%idg, mesh %i\n',oU,oP,mi);
lambda=lambda0;
id = tic;
[Uxh,Uyh,Ph,Eud(3)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PhNC,Vol);
tps(3)=toc(id);
dUx=Uxh-Uxex; dUy=Uyh-Uyex; dP=Ph-PexDG;
Eu0(3)=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/normUP2);
Eu1(3)=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/normUP2);
Eud(3)/=sqrt(normUP2);
Ep0(3)=sqrt((dP'*Mp*dP)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(3));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,pex)||_X = %7.2e\n',Eu0(3)); 
fprintf('        |Ph(Uex)-U_h|_1/||(Uex,pex)||_X = %7.2e\n',Eu1(3));
fprintf('         ||div(U_h)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud(3));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,pex)||_X = %7.2e\n',Ep0(3));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUATRIEME CALCUL : PRESSION NC + iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 4: P%i-P%idg, mesh %i\n',oU,oP,mi);
id = tic;
dPh0=Ph-PhNC;
[Uxh,Uyh,Ph,Eud(4),nit(4)]=SolveStokesExplicitTC_it(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Uxh,Uyh,DofU,DofUB,Ph,dPh0,Vol);
tps(4)=toc(id);
dUx=Uxh-Uxex; dUy=Uyh-Uyex; dP=Ph-PexDG;
Eu0(4)=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/normUP2);
Eu1(4)=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/normUP2);
Eud(4)/=sqrt(normUP2);
Ep0(4)=sqrt((dP'*Mp*dP)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(4),nit(4));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,pex)||_X = %7.2e\n',Eu0(4)); 
fprintf('        |Ph(Uex)-U_h|_1/||(Uex,pex)||_X = %7.2e\n',Eu1(4));
fprintf('         ||div(U_h)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud(4));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,pex)||_X = %7.2e\n',Ep0(4));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CINQUEME CALCUL : EF NC + PROJECTION SUR LES RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 5: P%inc-P%idg+RT, mesh %i\n',oU,oP,mi);
[RHS_UxNCRT,RHS_UyNCRT]=RHSStokesXYNCRT();
id = tic;
[UxhNCRT,UyhNCRT,PhNCRT,nit(5)]=UzawaNC(KuNC,BxNC,ByNC,MpNC,invMpNC,UnP,RHS_UxNCRT,RHS_UyNCRT,UxhNCB,UyhNCB,DofUNC,DofUNCB,Vol);
tps(5)=toc(id);
dUxNC=UxhNCRT-UxexNC; dUyNC=UyhNCRT-UyexNC; dPNC=PhNCRT-PexDG;
Eu0(5)=sqrt((dUxNC'*MuNC*dUxNC+dUyNC'*MuNC*dUyNC)/normUP2);
Eu1(5)=sqrt((dUxNC'*KuNC*dUxNC+dUyNC'*KuNC*dUyNC)/normUP2);
Eud(5)=ErrorDivUNC(UxhNCRT,UyhNCRT)/sqrt(normUP2);
Ep0(5)=sqrt((dPNC'*MpNC*dPNC)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(5),nit(5));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,Pex)||_X = %7.2e\n',Eu0(5)); 
fprintf('      ||Ph(Uex)-U_h||_h/||(Uex,Pex)||_X = %7.2e\n',Eu1(5));
fprintf('              ||U_h||_J/||(Uex,Pex)||_X = %7.2e\n',Eud(5));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,Pex)||_X = %7.2e\n',Ep0(5));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIXIEME CALCUL : PRESSION NC-RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 6:  P%i-P%idg, mesh %i\n',oU,oP,mi);
lambda0=lambda;
lambda=1;
id = tic;
[Uxh_6,Uyh_6,Ph_6,Eud(6)]=SolveStokesExplicitTC(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,UxhB,UyhB,DofU,DofUB,PhNCRT,Vol);
tps(6)=toc(id);
dUx=Uxh_6-Uxex; dUy=Uyh_6-Uyex; dP=Ph_6-PexDG;
Eu0(6)=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/normUP2);
Eu1(6)=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/normUP2);
Eud(6)/=sqrt(normUP2);
Ep0(6)=sqrt((dP'*Mp*dP)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, 1 it\n',tps(6));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,pex)||_X = %7.2e\n',Eu0(6)); 
fprintf('        |Ph(Uex)-U_h|_1/||(Uex,pex)||_X = %7.2e\n',Eu1(6));
fprintf('         ||div(U_h)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud(6));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,pex)||_X = %7.2e\n',Ep0(6));
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUATRIEME CALCUL : PRESSION NC + iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calcul 7: P%i-P%idg, mesh %i\n',oU,oP,mi);
id = tic;
dPh0=Ph_6-PhNCRT;
[Uxh_7,Uyh_7,Ph_7,Eud(7),nit(7)]=SolveStokesExplicitTC_it(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,Uxh,Uyh,DofU,DofUB,Ph_6,dPh0,Vol);
tps(7)=toc(id);
dUx=Uxh_7-Uxex; dUy=Uyh_7-Uyex; dP=Ph_7-PexDG;
Eu0(7)=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/normUP2);
Eu1(7)=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/normUP2);
Eud(7)/=sqrt(normUP2);
Ep0(7)=sqrt((dP'*Mp*dP)/nu2normUP2);
fprintf('-------\n');
fprintf('Temps de resolution (Ax=b) = %7.2e s, %i it\n',tps(7),nit(7));
fprintf('      ||Ph(Uex)-U_h||_0/||(Uex,pex)||_X = %7.2e\n',Eu0(7)); 
fprintf('        |Ph(Uex)-U_h|_1/||(Uex,pex)||_X = %7.2e\n',Eu1(7));
fprintf('         ||div(U_h)||_0/||(Uex,Pex)||_X = %7.2e\n',Eud(7));
fprintf('(1/nu)||Ph(Pex)-p_h||_0/||(Uex,pex)||_X = %7.2e\n',Ep0(7));
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
  tSVX=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%7.2e, mesh%i, Uxh',oU,oP,nu,lambda,mi);
  %
  texY=sprintf('LG P%i, mesh%i, Uy',oU,mi);
  tSVY=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%7.2e, mesh%i, Uyh',oU,oP,nu,lambda,mi);
  %
  texP=sprintf('LG P%i, mesh%i, Pr',oU,mi);
  tSVP=sprintf('SV P%i-P%i, nu=%7.2e, lambda=%7.2e, mesh%i, Prh',oU,oP,nu,lambda,mi);
  %
  tNCP=sprintf('NC P%i-P%i, nu=%7.2e, lambda=%7.2e, mesh%i, Prh',oU,oP,nu,lambda,mi);
  %
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Uxex);
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
  trisurf(NT,CN(:,1),CN(:,2),Uyex);
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
