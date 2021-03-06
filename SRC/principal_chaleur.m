% =====================================================
%
% principal_chaleur;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour 
% 1) l'equation de la chaleur suivante stationnaire, avec condition de
% Dirichlet non homogene
%
% | \alpha T - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2
% |         T = T_\Gamma,   sur le bord
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% 2) l'equation de la chaleur dependant du temps avec condition de 
% Dirichlet non homogene
%
% | dT/dt - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2 et pour tout t< t_max
% |         T = T_\Gamma,   sur le bord et pour tout t< t_max
% |         T = T_0       dans \Omega et pour t=0  
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure,
% T_0 est la valeur initiale de la temp?rature
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
% =====================================================
% Donnees du probleme
% ---------------------------------

% constantes:
h = 0.05;
alpha = 1;
lambda = 1;
T_Gamma = 290;

% génération et chargement du maillage et constantes
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) ' geomChaleur.geo']);
nom_maillage = 'geomChaleur.msh' ;
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% pour choisir quelle type de problème résoudre (en choisir uniquement un)
pb_stationnaire_Dirichlet   = 'non'; % exercice 1   (ne pas oublier de changer
pb_stationnaire_Fourier = 'non'; % exercice 2     les fonctions appropriées
pb_temporel               = 'oui'; % exercice 3        f, Tc, sigma_1 et sigma_2) 

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
SS = sparse(Nbpt,Nbpt); % matrice de surface
FF = zeros(Nbpt,1);     % vecteur approché pour f
TT_c = zeros(Nbpt,1);   % vecteur approché pour Tc

% construction des matrices MM et KK (commun aux 3 problèmes)
% ------------------------
% On itère sur tout les triangles
for l=1:Nbtri 
   % calcul des matrices elementaires du triangle l 
   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:),Reftri(l));
   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:));   
   % On fait l'assemblage des matrices globales
   for i = 1:3
      I = Numtri(l, i);
      for j = 1:3
         J = Numtri(l, j);
         MM(I,J) += Mel(i,j);
         KK(I,J) += Kel(i,j);
      end;
   end;
end

% calcul de f approché
for I = 1:Nbpt
    FF(I) = f(Coorneu(I,1), Coorneu(I,2));
end

% code spécifique pour résoudre le problème de Dirichlet (exercice 1)
if strcmp(pb_stationnaire_Dirichlet, 'oui')
  AA = alpha*MM+KK;
  LL = MM * FF;
  % on réalise la pseudo élimination seulement dans l'exercice 1
  [AA, LL] = elimine(AA, LL, Refneu);
  UU = AA\LL;
  UU = UU + T_Gamma;
end

% code spécifique pour résoudre le problème de Fourier (exercice 2)
if strcmp(pb_stationnaire_Fourier, 'oui')
  % calcul de Tc approché
  for I = 1:Nbpt
    TT_c(I) = Tc(Coorneu(I,1), Coorneu(I,2));
  end
  % boucle sur les arrêtes pour la construction de SS
  % ------------------------
  for l=1:Nbaretes
     [Sel]=mat_elem_surface(Coorneu(Numaretes(l,1),:), Coorneu(Numaretes(l,2),:));
     i = Numaretes(l,1);
     j = Numaretes(l,2);
     SS(i,i) += Sel(1,1);
     SS(i,j) += Sel(1,2);
     SS(j,i) += Sel(2,1);
     SS(j,j) += Sel(2,2);
  end;
  AA = alpha*MM+KK+lambda*SS;
  LL = MM * FF + lambda * SS * TT_c;
  UU = AA\LL;
end

##% code pour construire le vecteur solution exacte (pour vérification)
##for I = 1:Nbpt
##  x = Coorneu(I, 1);
##  y = Coorneu(I, 2);
##  TT(I) = cos(pi*x)*cos(pi*y);
##  %TT(I) = sigma_2(x,y);
##end;

% Si l'on souhaite afficher la valeur max de la solution
% disp(max(TT));

if strcmp(pb_stationnaire_Dirichlet, 'oui')
  % Calcul des erreurs L^1 et H^1 (valable pour le problème 1)
  % ----------
  UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
  % Calcul de l erreur L2
  sqrt( transpose(UU_exact)*MM*UU_exact+transpose(UU)*MM*UU-2*transpose(UU_exact)*MM*UU); 
  % Calcul de l erreur H1
  sqrt( transpose(UU_exact)*KK*UU_exact+transpose(UU)*KK*UU-2*transpose(UU_exact)*KK*UU);
  % attention de bien changer le terme source (dans FF)
end

% visualisation pour les problèmes stationnaires exercices 1 et 2
% -------------
if strcmp(pb_stationnaire_Dirichlet, 'oui') || strcmp(pb_stationnaire_Fourier, 'oui')
  affiche(UU, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
end
% =====================================================
% =====================================================
% Pour le probleme temporel
% ---------------------------------
if strcmp(pb_temporel,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    % on initialise la condition initiale
    % -----------------------------------
    T_initial = zeros(Nbpt,1);
    for I = 1:Nbpt
      T_initial(I) = condition_initiale(Coorneu(I,1), Coorneu(I,2));
    end
    
    
	  % solution a t=0
	  % --------------
    UU = T_initial - T_Gamma ;
    TT = UU + T_Gamma;


    % visualisation
    % -------------
    figure;
    hold on;
    affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(0)]);
    axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
        290,330,290,300]);
    hold off;

	% Boucle sur les pas de temps
	% ---------------------------
    for k = 1:N_t
        LL_k = zeros(Nbpt,1);
        % Calcul du second membre F a l instant k*delta t
        % -----------------------------------------------
        for I = 1:Nbpt
          FF(I) = f_t(Coorneu(I,1), Coorneu(I,2), k * delta_t);
        end
		    LL_k = MM*(FF + (1/delta_t)*UU);
        AA = alpha*MM + KK;
        [tilde_AA,tilde_LL_k] = elimine(AA, LL_k, Refneu);
		    % inversion
		    % ----------
        UU = tilde_AA\tilde_LL_k;
        TT = UU + T_Gamma;

        % visualisation 
		    % -------------
        pause(0.05)
        affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
        axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
            290,330,290,320]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020

