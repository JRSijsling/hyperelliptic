SetVerbose("IsGL2Equiv", 0);

counter := 0;
while counter lt 10^3 do
    printf ".";
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until not p in [2,5];
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    f := R.1^5 - 1;
    I := IgusaInvariants(f);

    Cs := TwistsFromIgusaInvariants(I);
    for i := 1 to #Cs do
        for j := i + 1 to #Cs do
            C1 := Cs[i];
            C2 := Cs[j];
            assert not IsIsomorphicHyperellipticCurves(C1, C2);
        end for;
    end for;

    fs := TwistsOfHyperellipticPolynomials(f);
    for i := 1 to #fs do
        for j := i + 1 to #fs do
            C1 := HyperellipticCurve(fs[i]);
            C2 := HyperellipticCurve(fs[j]);
            assert not IsIsomorphicHyperellipticCurves(C1, C2);
        end for;
    end for;

    repeat
        p := RandomPrime(5);
    until not p in [2,3,5,7];
    FF := FiniteField(p);
    R := PolynomialRing(FF);

    f := R.1^7 - 1;
    I := ShiodaInvariants(f);

    Cs := TwistsFromShiodaInvariants(I);
    for i := 1 to #Cs do
        for j := i + 1 to #Cs do
            C1 := Cs[i];
            C2 := Cs[j];
            assert not IsIsomorphicHyperellipticCurves(C1, C2);
        end for;
    end for;

    fs := TwistsOfHyperellipticPolynomials(f);
    for i := 1 to #fs do
        for j := i + 1 to #fs do
            C1 := HyperellipticCurve(fs[i]);
            C2 := HyperellipticCurve(fs[j]);
            assert not IsIsomorphicHyperellipticCurves(C1, C2);
        end for;
    end for;

    // Char 2?
    fs := TwistsOfHyperellipticPolynomials([ f, 0 ]);
    for i := 1 to #fs do
        for j := i + 1 to #fs do
            C1 := HyperellipticCurve(fs[i]);
            C2 := HyperellipticCurve(fs[j]);
            assert not IsIsomorphicHyperellipticCurves(C1, C2);
        end for;
    end for;

end while;
