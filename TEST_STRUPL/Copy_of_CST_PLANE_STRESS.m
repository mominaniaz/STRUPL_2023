% THIS PROGRAM USES AN 3-NODE LINEAR TRIANGULAR ELEMENT FOR THE
% LINEAR ELASTIC STATIC ANALYSIS OF A TWO DIMENSIONAL PROBLEM
%
clear all
clc
tic
%
% Make these variables global so they can be shared by other functions
%
global nnd nel nne nodof eldof n
global geom connec dee nf Nodal_loads
%
format long g
%
% ALTER NEXT LINES TO CHOOSE THE NAME OF THE OUTPUT FILE
%
fid =fopen('CST_COARSE_MESH_RESULTS.txt','w');
%
% To change the size of the problem or change elastic properties
% supply another input file
%
CST_COARSE_MESH_DATA;
%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assemble the global force vector
%
fg=zeros(n,1);
for i=1: nnd
if nf(i,1) ~= 0
fg(nf(i,1))= Nodal_loads(i,1);
end
if nf(i,2) ~= 0
fg(nf(i,2))= Nodal_loads(i,2);
end
end
%
% Assembly of the global stiffness matrix
%
% initialize the global stiffness matrix to zero
%
KK = zeros(n, n);
%
for i=1:nel
[bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
ke=thick*A*bee'*dee*bee; % Compute stiffness matrix
KK=form_KK(KK,ke, g); % assemble global stiffness matrix
end
%
%
%%%%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Compute global displacement vector and reaction%%
delta = KK\fg ; % solve for unknown displacements
Reaction=KK*delta-fg;
%
node_disp=zeros(nnd,2);
%
for i=1: nnd %
if nf(i,1) == 0 %
x_disp =0.; %
else
x_disp = delta(nf(i,1)); %
end
%
if nf(i,2) == 0 %
y_disp = 0.; %
else
y_disp = delta(nf(i,2)); %
end
node_disp(i,:) =[x_disp y_disp];
end
%
% Retrieve the x_coord and y_disp of the nodes located on the neutral axis
%
k = 0;
for i=1:nnd;
if geom(i,2)== 0.
k=k+1;
x_coord(k) = geom(i,1);
vertical_disp(k)=node_disp(i,2);
end
end
%% Compute stress and strain %%
%
for i=1:nel
[bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
eld=zeros(eldof,1); % Initialize element displacement to zero
for m=1:eldof
if g(m)==0 
    eld(m)=0.;
else %
eld(m)=delta(g(m)); % Retrieve element displacement
end
end
%
eps=bee*eld; % Compute strains
EPS(i,:)=eps ; % Store strains for all elements
sigma=dee*eps; % Compute stresses
SIGMA(i,:)=sigma ; % Store stress for all elements
end
%
%% Compute for each element the plastic multiplier that activates the plastic mode This part is a function related to yield function drucker pragher [alpha_b,Normal_stress]
%Macro-cracks
sigma_t = 0.1198; % Limit tension [MPa]
tau = 0.1198; % Limit shear [MPa]
tau_m = tau*3 ; %Valore di tau per macrofess
Gf_t = 0.035; % Fracture Energy traction [N/mm]
Gf_s = Gf_t ; % Fracture Energy shear[N/mm]
w_t = 2 * Gf_t/sigma_t ;
w_s = 2 * Gf_s/(tau) ;
h_s_t = -sigma_t /w_t  ; %softening parameter
h_s_s = -tau /w_s  ; %softening parameter
%Micro-cracks (Drucker-Prager parameters)
Plastic_mode = [1:1:nel]'; %Plastic mode
d_u = 0.019 ; %Ultimate deformation
coeff_h = 0.01 ; %percentuale di resistanza in corrispondenza del quale si attiva l'hardning
% h_t = ((1-coeff_h)*sigma_t /d_u) ; %hardening parameter
% h_s = ((1-coeff_h)*tau/d_u) ; %hardening parameter
h_t = (0.9*sigma_t /d_u) ; %hardening parameter


sigma_t_r = coeff_h * sigma_t ; % Limit traction
sigma_c = sigma_t_r  ; % Limit compression
%-------------------------
% sigma_t_tol = sigma_t_r + eps*sigma_t_r ;
% sigma_c_tol = sigma_c + eps*sigma_c ;
Aa = (2/sqrt(3))*((sigma_t_r*sigma_c)/(sigma_t_r+sigma_c)) ;
Bb = (1/sqrt(3))*((sigma_t_r-sigma_c)/(sigma_t_r+sigma_c)) ;
% Aa_tol = (2/sqrt(3))*((sigma_t_tol*sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;
% Bb_tol = (1/sqrt(3))*((sigma_t_tol-sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;


% %Volume of elements
% Volume = [];
% for i = 1 : Plastic_mode
%     Volume = [Volume; polyarea(ex(Plastic_mode(i),:),ey(Plastic_mode(i),:))*thick] ;
% end


for j = 1 : size(Plastic_mode,1)
    i = Plastic_mode(j) ;
    %Compute the principal stress 
    SIGMA_principal(i,1) = ((SIGMA(i,1)+SIGMA(i,2))/2) + (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %sigmaI
    SIGMA_principal(i,2) = ((SIGMA(i,1)+SIGMA(i,2))/2) - (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %sigmaII
    SIGMA_principal(i,3) = (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %taumax
    
    %Compute the first stress invariant
    I1(j,1) = SIGMA(i,1) + SIGMA(i,2) ;
    
    %Compute the  second invariant of the deviatoric stress it describes a state of pure shear.
    J2(j,1) = 1/6*(SIGMA(i,1)-SIGMA(i,2))^2 + 1/6*SIGMA(i,1)^2 + 1/6*SIGMA(i,2)^2 + SIGMA(i,3)^2 ;
    alpha_b(j,1) = Aa/(sqrt(J2(j))-Bb*I1(j)) ;
end
clear i j 

    % Solve the elastic problem amplifying the load  by multiplication of
    % fg *alpha to computes the normal stress
    %     [u_alpha]=solve(Kk,alpha_b(j)*F,bc) ;
   

eld_alpha = []; %element stress matrix
Normal_stress = [] ;
SIGMA_alpha = [] ;
for j = 1 : size(Plastic_mode,1)
    delta_alpha = KK\fg*alpha_b(j) ; % solve for unknown displacements
    Reaction_alpha=KK*delta_alpha-fg*alpha_b(j);
    %
    node_disp_alpha=zeros(nnd,2);
    %
    for i=1: nnd %
        if nf(i,1) == 0 %
            x_disp =0.; %
        else
            x_disp = delta_alpha(nf(i,1)); %
        end
        %
        if nf(i,2) == 0 %
            y_disp = 0.; %
        else
            y_disp = delta_alpha(nf(i,2)); %
        end
        node_disp(i,:) =[x_disp y_disp];
    end
    %
    % Retrieve the x_coord and y_disp of the nodes located on the neutral axis
    %
    k = 0;
    for i=1:nnd
        if geom(i,2)== 0.
            k=k+1;
            x_coord(k) = geom(i,1);
            vertical_disp(k)=node_disp(i,2);
        end
    end
    %% Compute stress and strain%% 
    %
    for i=1:nel
        [bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
        eld_alpha=zeros(eldof,1); % Initialize element displacement to zero
        for m=1:eldof
            if g(m)==0
                eld_alpha(m)=0.;
            else %
                eld_alpha(m)=delta_alpha(g(m)); % Retrieve element displacement
            end
        end
        %
        eps_alpha=bee*eld_alpha; % Compute strains
        EPS_alpha(i,:)=eps_alpha ; % Store strains for all elements
        sigma_alpha=dee*eps_alpha; % Compute stresses
        SIGMA_alpha(i,:)=sigma_alpha ; % Store stress for all elements
    end
    clear i

    i = Plastic_mode(j) ;

    %Compute the  second invariant of the deviatoric stress it describes a state of pure shear.
    J2_alpha = 1/6*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))^2 + 1/6*SIGMA_alpha(i,1)^2 + 1/6*SIGMA_alpha(i,2)^2 + SIGMA_alpha(i,3)^2 ;

    %Compute normal at a surface point 
    dphi_dsigmax = -Bb + 0.5*J2_alpha^(-0.5) * (1/3*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))+1/3*SIGMA_alpha(i,1)) ;
    dphi_dsigmay = -Bb + 0.5*J2_alpha^(-0.5) * (-1/3*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))+1/3*SIGMA_alpha(i,2)) ;
    dphi_dsigmaxy =  0.5*J2_alpha^(-0.5) * 2*SIGMA_alpha(i,3) ;
    %
    nnn = [dphi_dsigmax dphi_dsigmay dphi_dsigmaxy]  ;
    nnn_n = nnn/norm(nnn) ;
    Normal_stress = [Normal_stress; nnn_n] ;
end
clear j
%%  Computes Resistance vector Function%% 
crack_path=[];% this vector is the former IL vector from manuela and franchi codes is formed by 1st column prone crack nodes and 2dn column is the length a pre check loop is placed here where it is checked the limit tensio
% exceeded in the SIGMA_alpha array
%
 
RR = zeros(size(crack_path)) ;
RR(1:2:end) = crack_path(1:2:end)*thick*sigma_t;% find the resistance of nodes that failure in tension
RR(2:2:end) = crack_path(2:2:end)*thick*tau_m;%nodes that failure in shear

% for crack nodes the resistance is the multplication of stress by normal
% resistenze relative ai modi di fessurazione localizzata (prodotto della
% % tensione  limite per l'area di influenza di ciascun nodo, calcolata a

R_f = zeros(2*2*size(Cp,1),1) ;
for rr = 1 : size(RR,1)
    R_f(2*rr-1) = RR(rr) ;
end
clear rr

R_p = [] ;
for j = 1 : size(Plastic_mode,1)
    i = Plastic_mode(j) ;
    RR_p = SIGMA_alpha(j,:)* Normal_stress(j,:)' ;
    R_p = [R_p; RR_p] ;
end
clear j
%Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
R_p = R_p.*Vol ;

Resistance_vector = [R_f; R_p];
%% Assemble vector b 

% Print results to file
%
print_CST_results;
%
% Plot the stresses in the x_direction
%
x_stress = SIGMA(:,1);
cmin = min(x_stress);
cmax = max(x_stress);
clim([cmin cmax])
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',x_stress, ...
'Facecolor','flat','Marker','o')
colorbar;

toc