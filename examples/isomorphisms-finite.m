SetVerbose("IsGL2Equiv", 0);
//SetSeed(1340053158);
SetSeed(1662962995);

counter := 0;
while counter lt 10^3 do
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until not p eq 2;
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    repeat
        f1 := &+[ Random(FF)*R.1^i : i in [0..12] ];
        d := Degree(f1);
    until Discriminant(f1) ne 0 and d in [11,12];
    f1 /:= LeadingCoefficient(f1);

    repeat
        T := Matrix(FF, 2, 2, [ Random(FF) : i in [1..4] ]);
    until Determinant(T) ne 0;

    f2 := Evaluate(f1, (T[1,1]*R.1 + T[1,2])/(T[2,1]*R.1 + T[2,2]));
    f2 := (T[2,1]*R.1 + T[2,2])^12 * f2;
    f2 := R ! f2;
    f2 /:= LeadingCoefficient(f2);

    C1 := HyperellipticCurve(f1);
    C2 := HyperellipticCurve(f2);

    //printf ".";
    test, Ts := IsIsomorphicHyperellipticCurves(C1, C2);
    test, Ts := IsIsomorphicHyperellipticCurves(C1, C2 : geometric := true);

    found := false;
    for tup in Ts do
        U := tup[1];
        U := U^(-1);
        if WPSEqual([1,1,1,1], Eltseq(ChangeRing(T, BaseRing(U))), Eltseq(U)) then
            found := true;
        end if;
    end for;
    assert found;

    // Behavior of this function is a bit bizarre
    test, Ts := IsReducedIsomorphicHyperellipticCurves(C1, C2);
    //assert test;

end while;

// Matrix gives action on PP^1 covariantly
