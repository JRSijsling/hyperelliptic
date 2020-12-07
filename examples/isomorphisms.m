SetVerbose("IsGL2Equiv", 0);

P<x> := PolynomialRing(Rationals());
f1 := x^12 - 1;
C1 := HyperellipticCurve(f1); 
f2 := x^12 + 1;
C2 := HyperellipticCurve(f2);

time test, Ts := IsIsomorphicHyperelliptic(C1, C2);
time test, Ts := IsIsomorphicHyperelliptic(C1, C2 : geometric := true);

time test, Ts := IsIsomorphicHyperelliptic(C1, C2);
time test, Ts := IsIsomorphicHyperelliptic(C1, C2 : geometric := true, commonfield := true);

