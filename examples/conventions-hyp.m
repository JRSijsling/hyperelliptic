SetVerbose("IsGL2Equiv", 0);

counter := 0;
while counter lt 10^3 do
    printf ".";
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until p ne 2;
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    repeat
        f1 := R.1^12 + &+[ Random(FF)*R.1^i : i in [0..11] ];
        test := false;
        if f1 ne 0 then
            D := Discriminant(f1);
            test := D ne 0;
        end if;
    until test;

    repeat
        T := Matrix(FF, 2, 2, [ Random(FF) : i in [1..4] ]);
        U := Matrix(FF, 2, 2, [ Random(FF) : i in [1..4] ]);
    until Determinant(T)*Determinant(U) ne 0;

    TU := T*U;
    f2 := f1^(T*U);
    g2 := (f1^T)^U;
    if not (Degree(f1^T) eq 12 and Degree(f2) eq 12 and Degree(g2) eq 12) then continue; end if;
    assert f2 eq g2;

    X1 := HyperellipticCurve(f1);
    X2 := HyperellipticCurve(f2);

    test, Ts := IsIsomorphicHyperellipticCurves(X1, X2);
    test := false;
    for tup in Ts do
        T := tup[1];
        T := T^(-1);
        if WPSEqual([ 1 : i in [1..4] ], Eltseq(TU), Eltseq(T)) then
            test := true;
        end if;
    end for;
    assert test;
end while;
// Returns T such that f2^T = f1, so that in fact the curves transform via the left action
