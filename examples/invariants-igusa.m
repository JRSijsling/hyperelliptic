SetVerbose("IsGL2Equiv", 0);
//SetSeed(1340053158);
//SetSeed(960798546);

counter := 0;
while counter lt 10^3 do
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until not p eq 2;
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    repeat
        f := &+[ Random(FF)*R.1^i : i in [0..6] ];
        d := Degree(f);
        test := f ne 0;
        if test then
            test := Discriminant(f) ne 0 and d in [5, 6];
        end if;
    until test;

    C := HyperellipticCurve(f);
    I, W := IgusaInvariants(C);

    printf ".";
    //print I;
    //print DiscriminantFromIgusaInvariants(I);
    // Empty
    //print IgusaAlgebraicRelations(I);

    // No version for f
    g2s := G2Invariants(C);
    assert g2s eq IgusaToG2Invariants(I);
    J := G2ToIgusaInvariants(g2s);
    assert WPSEqual(W, I, J);

    D := HyperellipticCurveFromIgusaInvariants(I);
    J := IgusaInvariants(D);
    assert WPSEqual(W, I, J);
    assert IsIsomorphicHyperellipticCurves(C, D : geometric := true);

end while;


