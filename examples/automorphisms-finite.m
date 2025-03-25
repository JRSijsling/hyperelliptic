SetVerbose("IsGL2Equiv", 0);

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
    C1 := HyperellipticCurve(f1);

    printf ".";
    Ts := AutomorphismsOfHyperellipticCurve(C1);
    Ts := AutomorphismsOfHyperellipticCurve(C1 : geometric := true);

    for tup in Ts do
        T := tup[1];
        f2 := Evaluate(f1, (T[1,1]*R.1 + T[1,2])/(T[2,1]*R.1 + T[2,2]));
        f2 := (T[2,1]*R.1 + T[2,2])^12 * f2;
        f2 := R ! f2;

        g1 := f1/LeadingCoefficient(f1);
        g2 := f2/LeadingCoefficient(f2);

        assert g2 eq (Parent(g2) ! g1);
    end for;

    G := AutomorphismGroupOfHyperellipticCurve(C1);
    G := GeometricAutomorphismGroup(C1);

    Ts := ReducedAutomorphismsOfHyperellipticCurve(C1);
    Ts := ReducedAutomorphismsOfHyperellipticCurve(C1 : geometric := true);
    Ts := ReducedAutomorphismGroupOfHyperellipticCurve(C1);
    Ts := ReducedAutomorphismGroupOfHyperellipticCurve(C1 : geometric := true);
end while;
