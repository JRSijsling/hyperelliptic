//freeze;

/***
 *  Mini Toolboox for reconstructing genus 2 hyperelliptic curves with the
 *  conic and cubic method.
 *
 *  Distributed under the terms of the GNU Lesser General Public License (L-GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  Copyright 2011-2020, R. Lercier & C. Ritzenthaler
 */

import "g2conic_123.m" : Genus2ConicAndCubic123;

import "g2conic_uv1245.m" : Genus2ConicAndCubicUV1245;
import "g2conic_uv1246.m" : Genus2ConicAndCubicUV1246;


function FindPointOnConic(L : RationalPoint := true, RandomLine := true, Legendre := false, B := 100)
    /* B is the maximal height of the integral coefficients of the intersecting line. */

    K := BaseRing(Parent(L));
    P := ProjectiveSpace(K, 2); x := P.1; y := P.2; z := P.3;
    C := Conic(P, L);

    /* Can we find a rational point on this conic ? */
    if (Type(K) eq FldFin) or (RationalPoint
	and ((Type(K) in {FldRat, FldFin, RngInt}) or
	     (Type(K) eq FldAlg and (
	            Degree(K) eq 1 or IsSimple(K))) or
	     (Type(K) eq FldFunRat and (
		    Type(BaseField(K)) eq FldRat or
		    ISA(Type(BaseField(K)),FldNum) or
		    (IsFinite(BaseField(K)) and Characteristic(BaseField(K)) ne 2 ))) or

	     (ISA(Type(K), FldFunG) and Characteristic(K) ne 2)))
	then
	HasRatPoint, Pt := HasRationalPoint(C);
	if HasRatPoint then
            vprintf G2Twists, 1 : "Conic has a rational point\n";
	    return Parametrization(C, Place(Pt), Curve(ProjectiveSpace(K, 1)));
	end if;
	vprintf G2Twists, 1 : "Conic has no rational point\n";
    end if;

    /* Since we have no rational point on it, let us construct a quadratic extension that contains one */
    if Legendre then
        LC, LMmap := LegendreModel(C); LL := DefiningPolynomial(LC);
    else
        LC := C; LL := DefiningPolynomial(LC);
    end if;

    if RandomLine then
        D := [-B..B];
        repeat
            c1 := Random(D);
            c2 := Random(D);
            c3 := Random(D);
        until c2 ne 0;
        R<t> := PolynomialRing(K);
        h := hom<Parent(LL) -> R | [R.1, -(c1/c2)*R.1 - 1, 1]>;
        S := ExtensionField < K, t | h(LL) >;
        Pt := [ S | S.1, (-c1/c2)*S.1 - 1, 1];
    else
        a := K ! MonomialCoefficient(LL, x^2); b := K ! MonomialCoefficient(LL, x*z); c := K ! MonomialCoefficient(LL, z^2);
        S := ExtensionField < K, x | a*x^2 + b*x + c >;
        Pt := [ S | S.1, 0, 1];
    end if;

    CS := Conic(ProjectiveSpace(S, 2), L);
    if Legendre then
        Pt := [ Evaluate(p, Pt) : p in DefiningPolynomials(Inverse(LMmap)) ];
    end if;
    return Parametrization( CS, Place(CS!Pt), Curve(ProjectiveSpace(S, 1)) );

end function;

/* The affine case */
function GbToPoints(gb, tail)
    local up, spez, i, j, rt, f, res;

    if #gb eq 0 then return [tail]; end if;
    up := gb[#gb]; /* This polynomial is univariate */
    suc,up := IsUnivariate(gb[#gb], #gb);
    assert suc;  /* ideal is triangular and zero-dimensional */
    rt := Roots(up);
    res := [];
    for i := 1 to #rt do
        spez := [Evaluate(gb[j], #gb, rt[i][1]) : j in [1..#gb-1]];
        res := res cat $$(spez,[rt[i][1]] cat tail);
    end for;
    return res;
end function;

/*
Routines due to Stoll et al. for minimizing point clusters
apply in order to minimize the (conic, cubic) pair obtained in our algorithms.
*/

/* All points defined over nf will be returned */
function PointsOfIdeal(id, nf)
    local tri, erg, i, gb, up, r_nf;

    erg := [];
    if 1 in id then return erg; end if;

    tri := TriangularDecomposition(Radical(id));
    for i := 1 to #tri do
        gb := GroebnerBasis(tri[i]);
        r_nf := PolynomialRing(nf,Rank(Parent(gb[1])));
        erg := erg cat GbToPoints([r_nf!f : f in gb], []);
    end for;
    return erg;
end function;

function ReduceMestreConicAndCubic(f, Q)
    local hes,id,i,j,pl, CC, q, mt, subs, res, prec, mul;

    if false and Max([AbsoluteValue(a) : a in Coefficients(f)]) eq 1 then
        return f,[Parent(f).i : i in [1..3]];
    end if;

    mul := 3;
    prec := Round(100 + 4 * Log(1+Max([AbsoluteValue(i) : i in Coefficients(f)])));
    repeat
        mul := mul^2;
        CC := ComplexField(prec);

        id := ideal<PolynomialRing(RationalField(),3) | Q, f, Parent(f).1 - 1 >;
        pl := PointsOfIdeal(id,CC);

        /* Next two are usually not needed: */
        id := ideal<PolynomialRing(RationalField(),3) | Q, f, Parent(f).1, Parent(f).2-1 >;
        pl := pl cat PointsOfIdeal(id,CC);

        id := ideal<PolynomialRing(RationalField(),3) | Q, f, Parent(f).1, Parent(f).2, Parent(f).3-1 >;
        pl := pl cat PointsOfIdeal(id,CC);

        pts := [Vector([ChangePrecision(pl[i][j], prec div 2) : j in [1..3]]) : i in [1..#pl]];
        ptsnew,_,Trr := ReduceCluster(pts : eps := 10^(-prec div 20));

        subs := [&+[Parent(f).i * Trr[i,j] : i in [1..3]] : j in [1..3]];
        res :=  Evaluate(f,subs);

        prec := prec * 2;
    until Max([AbsoluteValue(a) : a in Coefficients(res)]) le mul * Max([AbsoluteValue(a) : a in Coefficients(f)]);

    if Max([AbsoluteValue(a) : a in Coefficients(res)]) ge Max([AbsoluteValue(a) : a in Coefficients(f)]) then
        return f,[Parent(f).i : i in [1..3]];
    end if;

    return Evaluate(f,subs),subs;
end function;

function MinimizeLinearEquationOverRationals(LE)

    u := Parent(LE).1; v := Parent(LE).2;

    an := Numerator(MonomialCoefficient(LE, u));
    ad := Denominator(MonomialCoefficient(LE, u));

    bn := Numerator(MonomialCoefficient(LE, v));
    bd := Denominator(MonomialCoefficient(LE, v));

    ct := GCD([an*bd, bn*ad, ad*bd]);

    a := (an*bd) div ct;
    b := (bn*ad) div ct;
    c := (ad*bd) div ct;

    _, U, V := ExtendedGreatestCommonDivisor(a, b);

    return U*c, V*c;

end function;

function PartialFactorization(N : Proof := false, PollardRhoLimit := 10^6, ECMB1 := 10^5, ECMCurves := 20)

    F, Q0 := TrialDivision(N, 10^6 : Proof := Proof);

    Q1 := [];
    for q in Q0 do
        _F, Q := PollardRho(q, 1, 1, PollardRhoLimit : Proof := Proof);
        F *:= _F; Q1 cat:= Q;
    end for;
    Q0 := Q1;

    Q1 := [];
    for q in Q0 do
        n := q;
        repeat
            nc := 0; repeat nc +:= 1;
                f := ECM(n, ECMB1);
            until f ne 0 or nc ge ECMCurves;
            if f ne 0 then
                repeat
                    if IsPrime(f : Proof := Proof) then
                        F *:= SeqFact([ <f, 1> ]);
                    else
                        Q1 cat:= [ f ];
                    end if;
                    n := n div f;
                until (n mod f) ne 0;

                if IsPrime(n : Proof := Proof) then
                    F *:= SeqFact([ <n, 1> ]);
                    n := 1;
                end if;
            end if;
        until n eq 1 or f eq 0;
        if n gt 1 then
            Q1 cat:= [n];
        end if;
    end for;

    return F, Q1;
end function;

function WPSMinimizeQQ(W, I);

    dens := [ Integers() ! Denominator(i) : i in I ];

    lambda := LCM(dens);

    Inorm := [Integers() | lambda^(W[k]) * I[k] : k in [1..#I] ];

    primes := [ fac[1] : fac in PartialFactorization(GCD(Inorm)) ];

    Imin := Inorm;
    for p in primes do
	while Seqset([Valuation(Imin[k], p) ge W[k] : k in [1..#Imin] ]) eq {true} do
	    Imin := [ Imin[k] div p^W[k] : k in [1..#Imin] ];
	end while;
    end for;

    return Imin;
end function;

function Genus2ConicAndCubic(JI : models := true, RationalModel := true, Deterministic := false)

    FF := Universe(JI);

    R := FF!0;
    if true and Type(FF) eq FldRat and RationalModel eq true then

        Puv := PolynomialRing(FF, 2); u := Puv.1; v := Puv.2;
        R := Puv!0;

        ret := true; sqrQ := 1;
        if #JI eq 6 then
            JI := WPSMinimizeQQ([2, 4, 6, 8, 10, 15], JI);
        else
            /* Let us recover J15 */
            JI := WPSMinimizeQQ([1, 2, 3, 4, 5], JI);
            J2, J4, J6, J8, J10 := Explode(JI);

            J15_2 := 2^22 * (
                J2^6*J6^3-2*J2^5*J4^2*J6^2+J2^4*J4^4*J6-72*J10*J2^5*J4*J6+8*J10*J2^4*J4^3-72*J2^4*J4*J6^3+136*J2^3*J4^3*J6^2-64*J2^2*J4^5*J6-432*J10^2*J2^5-48*J10*J2^4*J6^2+4816*J10*J2^3*J4^2*J6-512*J10*J2^2*J4^4 +216*J2^3*J6^4+1080*J2^2*J4^2*J6^3-2304*J2*J4^4*J6^2+1024*J4^6*J6+28800*J10^2*J2^3*J4-12960*J10*J2^2*J4*J6^2-84480*J10*J2*J4^3*J6+8192*J10*J4^5-7776*J2*J4*J6^4+6912*J4^3*J6^3-96000*J10^2*J2^2*J6-512000*J10^2*J2*J4^2-129600*J10*J2*J6^3+691200*J10*J4^2*J6^2+11664*J6^5+11520000*J10^2*J4*J6+51200000*J10^3
                );

            if J15_2 lt 0 then
                J2  *:= -1; J6  *:= -1; J8  *:= -1; J10 *:= -1; J15_2 *:= -1;
            end if;

            vprintf G2Twists, 2 :  "J15_2 := %o;\n", J15_2;

            /* Let us make J15_ 2 be a square (if not too hard) */
            F, Q := TrialDivision(Numerator(J15_2), 10^6 : Proof := false);
            if #Q gt 0 then
                ret, sqrQ := IsSquare(&*Q);
            end if;

            if ret eq false then
                F, Q := PartialFactorization(Numerator(J15_2));
                if #Q gt 0 then
                    ret, sqrQ := IsSquare(&*Q);
                end if;
            end if;

            vprintf G2Twists, 2 :  "J15_2 := %o * %o;\n", F, Q;

            if #Q gt 0 then
                ret, sqrQ := IsSquare(&*Q);
            end if;
            if ret then
                for fct in F do
                    if fct[2] mod 2 eq 1 then
                        J2  *:= fct[1]; J4  *:= fct[1]^2; J6  *:= fct[1]^3;
                        J8  *:= fct[1]^4; J10 *:= fct[1]^5; J15_2 *:= fct[1]^15;
                    end if;
                end for;
                ret, J15 := IsSquare(AbsoluteValue(J15_2));
            end if;
            if ret then JI := [J2, J4, J6, J8, J10, J15]; end if;
        end if;

        if ret then
            vprintf G2Twists, 2 :  "[J2, J4, J6, J8, J10, J15] = %o\n\n", JI;

            Genus2ConicAndCubicUVFct := Genus2ConicAndCubicUV1245;
            R, _, _ := Genus2ConicAndCubicUVFct([u, v], JI : models := false);

            if R eq 0 or Degree(R,u)*Degree(R,v) eq 0 then
                Genus2ConicAndCubicUVFct := Genus2ConicAndCubicUV1246;
                R, _, _ := Genus2ConicAndCubicUVFct([u, v], JI : models := false);
            end if;

        end if;

        /* Let us find a conic with small discriminant */
	if R ne 0 and Degree(R,u)*Degree(R,v) ne 0 then

	    vprintf G2Twists, 2 :  "Let us minimize the discriminant of the conic to be used, i.e %o\n\n", R;

	    U, V := MinimizeLinearEquationOverRationals(R);

	    vprintf G2Twists, 2 :  "We set :";
	    vprintf G2Twists, 2 :  "  u = %o\n", U;
	    vprintf G2Twists, 2 :  "  v = %o\n", V;

	    R, C, M := Genus2ConicAndCubicUVFct([FF | U, V], JI : models := models);

	    vprintf G2Twists, 2 :  "So that :";
	    vprintf G2Twists, 2 :  "  R = %o\n", R;

	    /* Let us first remove the content of C and M */
	    ct := LCM([Denominator(c) : c in Coefficients(C)]) /
		  GCD([Numerator(c) : c in Coefficients(C)]);
	    C *:= ct;
	    ct := LCM([Denominator(c) : c in Coefficients(M)]) /
		GCD([Numerator(c) : c in Coefficients(M)]);
	    M *:= ct;

            vprintf G2Twists, 2 :
               "Factorization of conic discriminant before reduction: %o\n",
                Factorization(Integers() ! Discriminant(Conic(ProjectiveSpace(Parent(C)), C)));

            vprintf G2Twists, 2 : "Minimal model step...\n";
	    Cphi, phi := MinimalModel(Conic(ProjectiveSpace(Parent(C)), C));
	    C := DefiningPolynomial(Cphi);
	    M := Evaluate(M, DefiningPolynomials(phi));
	    ct := GCD([Denominator(c) : c in Coefficients(M)]) /
		  GCD([Numerator(c) : c in Coefficients(M)]);
	    M *:= ct;
	    vprintf G2Twists, 2 :  "Conic %o\n", C;
	    vprintf G2Twists, 2 :  "Cubic %o\n", M;

	else
	    R := FF!0;
	    vprintf G2Twists, 2 :  "None parametric conic works ! Let us do it as usual...\n\n", R;
	end if;
    end if;

    if R eq 0 then R, C, M := Genus2ConicAndCubic123(JI : models := models); end if;

    /* Hum, hum... No non-degenrate conic found, return immediatly */
    if R eq 0 then return false; end if;

    /* Computing conic and cubic */
    if models then

	phi := FindPointOnConic(C : RationalPoint := RationalModel, RandomLine := not Deterministic);

	f := Evaluate(M, DefiningPolynomials(phi));

        g := UnivariatePolynomial(Evaluate(f, Parent(f).2, 1));

        vprintf G2Twists, 1 :  "Hyperelliptic polynomial: %o\n", g;
	return g;

    end if;

    return true;

end function;
