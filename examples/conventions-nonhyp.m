SetVerbose("IsGL2Equiv", 0);

counter := 0;
while counter lt 10^3 do
    printf ".";
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until p ne 3;
    FF := FiniteField(p);
    R := PolynomialRing(FF, 3);
    PP2 := ProjectiveSpace(R);

    repeat
        F1 := &+[ Random(FF)*mon : mon in MonomialsOfDegree(R, 4) ];
        test := false;
        if F1 ne 0 then
            D := DiscriminantOfTernaryQuartic(F1);
            test := D ne 0;
        end if;
    until test;

    repeat
        T := Matrix(FF, 3, 3, [ Random(FF) : i in [1..9] ]);
        U := Matrix(FF, 3, 3, [ Random(FF) : i in [1..9] ]);
    until Determinant(T)*Determinant(U) ne 0;

    TU := T*U;
    F2 := F1^(T*U);
    G2 := (F1^T)^U;
    assert Curve(PP2, F2) eq Curve(PP2, G2);

    X1 := Curve(PP2, F1);
    X2 := Curve(PP2, F2);

    test, Ts := IsIsomorphicPlaneQuartics(X1, X2);
    test := false;
    for T in Ts do
        Ti := T^(-1);
        if WPSEqual([ 1 : i in [1..9] ], Eltseq(TU), Eltseq(Ti)) then
            test := true;
        end if;
    end for;
    assert test;
end while;
// Returns T such that F2^T = F1, so that in fact the curves transform via the left action
