SetVerbose("IsGL2Equiv", 0);

counter := 0;
while counter lt 10^3 do
    counter +:= 1;

    repeat
        p := RandomPrime(5);
    until not p eq 3;
    //p := 2;

    FF := FiniteField(p);
    R := PolynomialRing(FF, 3);
    repeat
        F := &+[ Random(FF)*mon : mon in MonomialsOfDegree(R, 4) ];
        test := false;
        if F ne 0 then
            D := DiscriminantOfTernaryQuartic(F);
            test := D ne 0;
        end if;
    until test;

    printf ".";
    C := Curve(ProjectiveSpace(R), F);
    Cs := Twists(C);
    for i := 1 to #Cs do
        for j := i + 1 to #Cs do
            C1 := Cs[i];
            C2 := Cs[j];
            assert not IsIsomorphicPlaneQuartics(C1, C2);
        end for;
    end for;

end while;
