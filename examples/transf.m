QQ := Rationals();
R<x,y,z> := PolynomialRing(QQ, 3);
PP2 := ProjectiveSpace(R);

T := Matrix(QQ, [[1,1,0],[0,1,0],[0,0,1]]);
P := [2,3,4];
TP := [5,3,4];
f := z - 2*x;
f := y - (3/2)*x;
fT := f^(T^(-1));

print Evaluate(fT, TP);
print x^T;
print y^T;
print z^T;

// So this is straight up the transformation needed, in the following sense:
// If f2 = f1^(T^(-1)), then T is the transformation on the ambient
// So f1 = f2^T
