SetVerbose("IsGL2Equiv", 0);
//SetSeed(1340053158);
SetSeed(2698839856);

counter := 0;
while counter lt 10^3 do
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until not p in [2,3,5,7];
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    repeat
        f := &+[ Random(FF)*R.1^i : i in [0..8] ];
        d := Degree(f);
    until Discriminant(f) ne 0 and d in [7, 8];
    f /:= LeadingCoefficient(f);

    C := HyperellipticCurve(f);

    S, W := ShiodaInvariants(C);
    // Transvectant problem
    //M := MaedaInvariants(C);

    print "";
    print S;
    print DiscriminantFromShiodaInvariants(S);

    D := HyperellipticCurveFromShiodaInvariants(S);
    T := ShiodaInvariants(D);
    assert WPSEqual(W, S, T);
    assert IsIsomorphicHyperellipticCurves(C, D : geometric := true);

end while;


