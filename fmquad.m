function[der,fun] = fmquad(samp,Element_Type,dim, ig,jg)%,lg
%
% This function returns the vector of the shape function and their
% derivatives with respect to xi and eta at the gauss points for
% an 8-nodded quadrilateral
%
% Element_type=8;
% Element_type=3;
% Element_type=4;
% dim=2;
% dim=3;
xi=samp(ig,1);
eta=samp(jg,1);
% zeta=samp(lg,1);
etam=(1.-eta);
etap=(1.+eta);
xim=(1.-xi);
xip=(1.+xi);
% zetam=(1.-zeta);
% zetap=(1.+zeta);

% 
 if Element_Type==8 && dim==2
% shape functions
fun(1) = -0.25*xim*etam*(1.+ xi + eta);
fun(2) = 0.5*(1.- xi^2)*etam;
fun(3) = -0.25*xip*etam*(1. - xi + eta);
fun(4) = 0.5*xip*(1. - eta^2);
fun(5) = -0.25*xip*etap*(1. - xi - eta);
fun(6) = 0.5*(1. - xi^2)*etap;
fun(7) = -0.25*xim*etap*(1. + xi - eta);
fun(8) = 0.5*xim*(1. - eta^2);
% derivatives
der(1,1)=0.25*etam*(2.*xi + eta); der(1,2)=-1.*etam*xi;
der(1,3)=0.25*etam*(2.*xi-eta); der(1,4)=0.5*(1-eta^2);
der(1,5)=0.25*etap*(2.*xi+eta); der(1,6)=-1.*etap*xi;
der(1,7)=0.25*etap*(2.*xi-eta); der(1,8)=-0.5*(1.-eta^2);
%
der(2,1)=0.25*xim*(2.*eta+xi); der(2,2)=-0.5*(1. - xi^2);
der(2,3)=-0.25*xip*(xi-2.*eta); der(2,4)=-1.*xip*eta;
der(2,5)=0.25*xip*(xi+2.*eta); der(2,6)=0.5*(1.-xi^2);
der(2,7)=-0.25*xim*(xi-2.*eta); der(2,8)=-1.*xim*eta;
%
 elseif  Element_Type==4 && dim==2
% shape functions
 fun(1)=0.25*xim*etam;
 fun(2)=0.25*xip*etam;
 fun(3)=0.25*xip*etap;
 fun(4)=0.25*xim*etap;
% derivatives
 der(1,1)=-0.25*etam;  der(1,2)=-0.25*xim;
 der(2,1)=0.25*etam;  der(2,2)=-0.25*xip;
 der(3,1)=0.25*etap;  der(3,2)=0.25*xip;
 der(4,1)=-0.25*etap; der(4,2)=0.25*xim;
 end

 if Element_Type==3 && dim==2
% shape functions
    fun(1) = 1. - xi - eta;
    fun(2) =  xi;
    fun(3) =  eta; 
% derivatives
    der(1,1)= -1.; der(1,2)=-1.;
    der(2,1)= 1.; der(2,2)= 0;
    der(3,1)=0; der(3,2)=1.;

 end
 
  if Element_Type==8 && dim==3  %brick element
% shape functions
fun(1) = 0.125*xim*etam*zetam;
fun(2) = 0.125*xip*etam*zetam;
fun(3) = 0.125*xip*etap*zetam;
fun(4) = 0.125*xim*etap*zetam;
fun(5) = 0.125*xim*etam*zetap;
fun(6) = 0.125*xip*etam*zetap;
fun(7) = 0.125*xip*etap*zetap;
fun(8) = 0.125*xim*etap*zetap;
% derivatives
der(1,1)=0.125*etam*zetam; der(1,2)=0.125*etam*zetam; 
der(1,3)=0.125*etap*zetam; der(1,4)=0.125*etap*zetam;
der(1,5)=0.125*etam*zetap; der(1,6)=0.125*etam*zetap;
der(1,7)=0.125*etap*zetap; der(1,8)=0.125*etap*zetap;
%
der(2,1)=0.125*xim*zetam; der(2,2)=0.125*xip*zetam;
der(2,3)=0.125*xip*zetam; der(2,4)=0.125*xim*zetam;
der(2,5)=0.125*xim*zetap; der(2,6)=0.125*xip*zetap;
der(2,7)=0.125*xip*zetap; der(2,8)=0.125*xim*zetap;
%
der(3,1)=0.125*xim; der(3,2)=0.125*xip;
der(3,3)=0.125*xip; der(3,4)=0.125*xim;
der(3,5)=0.125*xim; der(3,6)=0.125*xip;
der(3,7)=0.125*xip; der(3,8)=0.125*xim;
  end
% 
der'
% 

% end function fmquad