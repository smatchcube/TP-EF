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
h = 0.05;
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) ' geomChaleur.geo']);
nom_maillage = 'geomChaleur.msh' ;

validation = 'oui';
pb_stationnaire = 'non';
pb_temporel = 'non';

if strcmp(validation,'oui')
    alpha = 1;
    T_Gamma = 0;
end

if strcmp(pb_stationnaire,'oui')
    alpha = 1;
    T_Gamma = 290;
end

if strcmp(pb_temporel,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 290;
end

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  
  % calcul des matrices elementaires du triangle l 
  
   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:),Reftri(l));
   % LA ROUTINE matK_elem.m DOIT ETRE MODIFIEE

   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:));
    
    % On fait l'assemblage de la matrice globale
    % A COMPLETER
end % for l

% Matrice EF
% -------------------------
AA = alpha*MM+KK;

% =====================================================
% =====================================================
% Pour le probleme stationnaire et la validation
% ---------------------------------

% Calcul du second membre F
% -------------------------
% A COMPLETER EN UTILISANT LA ROUTINE f.m
FF = f(...);
LL = ...;

% inversion
% ----------
% tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
% APRES PSEUDO_ELIMINATION 
% ECRIRE LA ROUTINE elimine.m ET INSERER L APPEL A CETTE ROUTINE
% A UN ENDROIT APPROPRIE
UU = tilde_AA\tilde_LL;
TT = ...;

% validation
% ----------
if strcmp(validation,'oui')
    UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
	% Calcul de l erreur L2
	% A COMPLETER
	% Calcul de l erreur H1
	% A COMPLETER
	% attention de bien changer le terme source (dans FF)
end

% visualisation
% -------------
if ( strcmp(validation,'oui') || strcmp(pb_stationnaire,'oui') )
    affiche(TT, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
end

% =====================================================
% =====================================================
% Pour le probleme temporel
% ---------------------------------
if strcmp(pb_temporel,'oui')

    % on initialise la condition initiale
    % -----------------------------------
    T_initial = condition_initiale(Coorneu(:,1),Coorneu(:,2));

	% solution a t=0
	% --------------
    UU = ...;
    TT = ...;

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
		% A COMPLETER EN UTILISANT LA ROUTINE f_t.m et le terme precedent (donne par UU)
		LL_k = ...;

		% inversion
		% ----------
		% tilde_AA ET tilde_LL_k SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
		% APRES PSEUDO_ELIMINATION 
		% ECRIRE LA ROUTINE elimine.m ET INSERER L APPEL A CETTE ROUTINE
		% A UN ENDROIT APPROPRIE
        UU = tilde_AA\tilde_LL_k;
        TT = ...;

        % visualisation 
		& -------------
        pause(0.05)
        affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
        axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
            290,330,290,320]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020

