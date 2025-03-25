k:=GF(47);
P<x,y,z>:=PolynomialRing(k,3);
PP:=ProjectiveSpace(P);
a:=Random(k);
b:=Random(k);
f:=x*(y^3+z^3)+y^2*z^2+a*x^2*y*z+b*x^4;
C:=Curve(PP,f);
if IsSingular(C) eq false then
print #GeometricAutomorphismGroup(C);
end if;
Twists(C);
