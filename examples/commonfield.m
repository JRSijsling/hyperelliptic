import "magma/TernaryForms.m": Dehomogenization,Homogenization;

SetSeed(1);
SetVerbose("IsGL2Equiv", 0);

F := FiniteField(NextPrime(10^6));
F := FiniteField(NextPrime(10^6), 3);
F := Rationals();
F := NumberField(x^2+x+1);
_<x> := PolynomialRing(F);

deg := 50; B := 100; D := [-B..B];
while true do
    repeat
        f := &+[ Random(D)*x^i : i in [0..deg] ];
    until Degree(f) eq deg and Discriminant(f) ne 0;
    repeat
        T := Matrix(F, 2, 2, [ Random(D) : i in [1..4] ]);
    until Determinant(T) ne 0;
    g := f^T;
    print "Covariant:";
    time _, IsoLst := IsGL2GeometricEquivalent(f, g, deg : geometric := false);
    print "Direct:";
    time _, IsoLst := IsGL2GeometricEquivalent(f, g, deg : geometric := false, covariant := false);
    print Universe(IsoLst[1]);
end while;

f := x^5 + x; g := f;
f := 4*x^6 + 1; g := f;
f := x^6 + x^5 + 2*x^3 + 5*x^2 + x + 1; g := f;

while true do
    _, IsoLst := IsGL2GeometricEquivalent(f, f, deg : commonfield := false);
    _, IsoLst := IsGL2GeometricEquivalent(f, f, deg : commonfield := false, covariant := false);
    print "Covariant:";
    time _, IsoLst := IsGL2GeometricEquivalent(f, f, deg : commonfield := true);
    print "Direct:";
    time _, IsoLst := IsGL2GeometricEquivalent(f, f, deg : commonfield := true, covariant := false);
    print #IsoLst;
    print Universe(IsoLst[1]);
end while;
