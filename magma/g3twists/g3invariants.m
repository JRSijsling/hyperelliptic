//freeze;

/***
 *  Invariants of Genus 3 Curves.
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
 *  Copyright 2007-2020 R. Lercier & C. Ritzenthaler & J. Sijsling
 */

/***
 * Exported intrinsics.
 *
 * intrinsic ShiodaInvariants(f::RngUPolElt : normalize := false) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariants(fh::SeqEnum : normalize := false) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariants(C::CrvHyp : normalize := false) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariantsEqual(V1::SeqEnum, V2::SeqEnum) -> BoolElt
 * intrinsic ShiodaAlgebraicInvariants(PrimaryInvariants::SeqEnum : ratsolve := true) -> SeqEnum
 *
 * intrinsic MaedaInvariants(f::RngUPolElt) -> SeqEnum
 * intrinsic MaedaInvariants(C::CrvHyp) -> SeqEnum
 *
 ********************************************************************/

import "g3twists_char2.m" : ShiodaInvariantsChar2, ShiodaAlgebraicInvariantsChar2;
import "g3twists_char3.m" : ShiodaInvariantsChar3, ShiodaAlgebraicInvariantsChar3;
import "g3twists_char5.m" : ShiodaInvariantsChar5;
import "g3twists_char7.m" : ShiodaInvariantsChar7, ShiodaAlgebraicInvariantsChar7;
import "g3twists_charp.m" : ShiodaInvariantsCharp, ShiodaAlgebraicInvariantsCharp;

/***
 *
 * Shioda invariants, in fields of characteristic 0 or > 2,
 * see [Shioda1967].
 *
 ********************************************************************/
intrinsic ShiodaInvariants(f::RngUPolElt, p::RngIntElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum

    {Compute the Shioda invariants  'J2', 'J3', 'J4', 'J5', 'J6', 'J7',
    'J8', 'J9' and 'J10' of a polynomial of degree at most 8, considered as a
    binary form of degre 8 (see 1967 Shioda's  paper), and other invariants in
    characteristic 2, 3, 5 and 7. Weights of these
    invariants are returned too.
    }

    CR := CoefficientRing(Parent(f));
    require not Characteristic(CR) in {2} and not p in {2} :
        "2 must be invertible in the base ring.";
    require Degree(f) le 8:
        "Polynomial must have degree at most 8.";

    /*
    require not (p eq 5 and PrimaryOnly eq true) :
        "Primary invariants are not known in char 5";
    */

    /* Rings of small characteristic  */
    if p eq 3 then
	JI, Wght := ShiodaInvariantsChar3(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 5 then
	JI, Wght := ShiodaInvariantsChar5(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 7 then
	JI, Wght := ShiodaInvariantsChar7(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    else
        scaled := true; if IntegralNormalization then scaled := false; end if;

        JI, Wght := ShiodaInvariantsCharp(f : PrimaryOnly := PrimaryOnly, scaled := scaled, degmax := degmax, degmin := degmin);

        if IntegralNormalization then
            Scaling := [
                2^5 * 3^2,
                2^8 * 3^4 * 5,
                2^10 * 3^4 * 5^2,
                2^14 * 3^6 * 5^2,
                2^16 * 3^7 * 5^3,
                2^20 * 3^9 * 5^3,
                2^22 * 3^10 * 5^4,
                2^26 * 3^12 * 5^4,
                2^28 * 3^13 * 5^5
                ];
            for d := 2 to 10 do
                if d ge degmin and d le degmax then
                    JI[d-1] := Universe(JI)!(JI[d-1] / Scaling[d-1]);
                end if;
            end for;
        end if;

    end if;

    if normalize eq false then return JI, Wght; end if;

    w := GCD(Wght);
    return WPSNormalize([e div w : e in Wght], JI), Wght;

end intrinsic;

intrinsic ShiodaInvariants(f::RngMPolElt, p::RngIntElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum

    {Compute the Shioda invariants  'J2', 'J3', 'J4', 'J5', 'J6', 'J7',
    'J8', 'J9' and 'J10' of a polynomial of degree at most 8, considered as a
    binary form of degre 8 (see 1967 Shioda's  paper), and other invariants in
    characteristic 2, 3, 5 and 7. Weights of these
    invariants are returned too.
    }

    P := Parent(f);
    require
        Rank(P) eq 2 and {Degree(e) : e in Monomials(f)} eq {8} :
        "Input must be a binary homogeneous polynomial of degree 8";

    F := UnivariatePolynomial(Evaluate(f, 2, 1));

    return ShiodaInvariants(F, p :
        normalize := normalize,
        IntegralNormalization := IntegralNormalization,
        PrimaryOnly := PrimaryOnly,
        degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(f::RngUPolElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum

    {Compute the Shioda invariants  'J2', 'J3', 'J4', 'J5', 'J6', 'J7',
    'J8', 'J9' and 'J10' of a polynomial of degree at most 8, considered as a
    binary form of degre 8 (see 1967 Shioda's  paper), and other invariants in
    characteristic 2, 3, 5 and 7. Weights of these
    invariants are returned too.
    }

    require Degree(f) le 8: "Polynomial must have degree at most 8.";

    return ShiodaInvariants(f, Characteristic(Parent(f)) :
        normalize := normalize,
        IntegralNormalization := IntegralNormalization,
        PrimaryOnly := PrimaryOnly,
        degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(f::RngMPolElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum

    {Compute the Shioda invariants  'J2', 'J3', 'J4', 'J5', 'J6', 'J7',
    'J8', 'J9' and 'J10' of a polynomial of degree at most 8, considered as a
    binary form of degre 8 (see 1967 Shioda's  paper), and other invariants in
    characteristic 2, 3, 5 and 7. Weights of these
    invariants are returned too.
    }

    P := Parent(f);
    require
        Rank(P) eq 2 and {Degree(e) : e in Monomials(f)} eq {8} :
        "Input must be a binary homogeneous polynomial of degree 8";

    F := UnivariatePolynomial(Evaluate(f, 2, 1));

    return ShiodaInvariants(F :
        normalize := normalize,
        IntegralNormalization := IntegralNormalization,
        PrimaryOnly := PrimaryOnly,
        degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(fh::SeqEnum :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
    {Let fh be equal to [f(x), h(x)]. Compute the invariants of the hyperelliptic curve y^2 + h(x) * y = f(x)
    (see 1967 Shioda's  paper). Weights of these invariants are returned
    too. In rings of characteristic 0 or >= 11, these invariants are  'J2',
    'J3', 'J4', 'J5', 'J6', 'J7',  'J8', 'J9' and 'J10', of weight 2, 3.. 10.
    Other invariants are computed in characteristic 2, 3, 5 and 7.
    }

    f, h := Explode(fh);
    K := CoefficientRing(Parent(f));

    if Characteristic(K) eq 2 then
        JI, Wght := ShiodaInvariantsChar2(f, h :
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin);

        if normalize eq false then return JI, Wght; end if;

        w := GCD(Wght);
        return WPSNormalize([e div w : e in Wght], JI), Wght;

    end if;

    if h eq 0 then F := f; else F := 4*f+h^2; end if;

    return ShiodaInvariants(F :
            normalize := normalize,
            IntegralNormalization := IntegralNormalization,
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(C::CrvHyp :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false,
    degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
    {
    Compute geometric invariants of a genus 3 hyperelliptic curve, i.e. 'J2',
    'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9' and 'J10' in characteristic p = 0
    or p >= 11 (see 1967 Shioda's paper), and other invariants in
    characteristic 2, 3, 5 and 7.
    }

    require Genus(C) eq 3: "Curve must be of genus 3.";
    K := BaseField(C);
    R := PolynomialRing(K); x := R.1;
    f, h := HyperellipticPolynomials(C);
    require (Degree(h) le 4) and (Degree(f) le 8):
	"The curve must be of form y^2 + h(x) y = f(x), where f and h must have suitable degrees.";

    if Characteristic(K) eq 2 then
        JI, WG := ShiodaInvariants([f, h] : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

        if normalize eq false then return JI, WG; end if;

        w := GCD(WG);
        return WPSNormalize([e div w : e in WG], JI), WG;

    end if;

    return ShiodaInvariants([f, h] :
            normalize := normalize,
            IntegralNormalization := IntegralNormalization,
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin
            );

end intrinsic;

intrinsic ShiodaAlgebraicInvariants(PrimaryInvariants::SeqEnum : ratsolve := true) -> SeqEnum
    {
    This function returns the algebraic relations between the six primary
    invariants given in input and the other invariants.

    In Characteristic 0 or > 7, the primary invariants are J2, J3, J4, J5, J6
    and J7, in characteristic 3, they are J2, J4, J5, J7, J9, J12.

    By default (ratsolve := true), this function computes solutions to
    the system of relations and returns them as a list of invariants. Otherwise
    (ratsolve := false), the system is returned as is.
    }

    FF := Universe(PrimaryInvariants);
    p := Characteristic(FF);

    /* Not yet implemented */
    if p in {5} then
	error "[Hyperelliptic] currently, no computation done in fields of char. 5, sorry";
	return [];
    end if;

    /* Rings of small characteristic  */
    case p:

    when 2:
	require #PrimaryInvariants eq 6 or #PrimaryInvariants eq 5 or #PrimaryInvariants eq 4 or #PrimaryInvariants eq 3
	    : "Argument must be a sequence of size 3 (N7, N32, N40), 4 (j2, L, L', L''), 5 (j2<>0, K, K', K'', K'''') or (j2=0, M1, M3, M4, M5) or 6 (j2, j3, J5, J6, J8, J9) invariants";

	return ShiodaAlgebraicInvariantsChar2(PrimaryInvariants, ratsolve);

    when 3:
	require #PrimaryInvariants eq 6:
	    "Argument must be a sequence of six free Shioda invariants J2, J4, J5, J7, J9, J12.";
	return ShiodaAlgebraicInvariantsChar3(PrimaryInvariants, ratsolve);

    when 7:
	require #PrimaryInvariants eq 6:
	    "Argument must be a sequence of six free Shioda invariants J2, J4, J5, J7, J9, J12.";
	return ShiodaAlgebraicInvariantsChar7(PrimaryInvariants, ratsolve);

    end case;

    require #PrimaryInvariants eq 6:
	"Argument must be a sequence of six free Shioda invariants J2, J3, J4, J5, J6, J7.";

    if Universe(PrimaryInvariants) cmpeq Integers() then
	ChangeUniverse(~PrimaryInvariants, Rationals());
    end if;

    return ShiodaAlgebraicInvariantsCharp(PrimaryInvariants, ratsolve);

end intrinsic;


 /***
  *
  * Maeda invariants in characteristic 0 or > 7
  * see [Maeda1997].
  *
  ********************************************************************/

intrinsic MaedaInvariants(f::RngUPolElt) -> SeqEnum
    {Compute the Maeda field invariants  'I2', 'I3', 'I4', 'I4p', 'I8' and
    'I9' of a polynomial of degree at most 8, considered as a binary form of
    degre 8 (see 1990 Maeda's paper). It works only in characteristic p = 0
    or p >= 11.
    }

    require Degree(f) le 8: "Polynomial must have degree at most 8.";

    CR := CoefficientRing(Parent(f));
    require not Characteristic(CR) in {2, 3, 5, 7} : "2, 3, 5 and 7 must be invertible in the base ring.";

    f_cov := [* f, 1, 8 *];

    /* Q=(f,f)_6 : order 4, degree 2 */
    Q:=Transvectant(f_cov,f_cov,6);

    /* t=((Q,Q)_2,Q)_1 : order 6, degree 6 */
    t:=Transvectant(Transvectant(Q,Q,2),Q,1);

    /* theta=(f,t)_6 : order 2, degree 7 */
    theta:=Transvectant(f_cov,t,6);

    /* M=(f,f)_6 : order 0, degree 12 */
    M:=Transvectant(t,t,6);
    m:=M[1];

    /* j=((t,t)_2,t)_1 : order 12, degree 18 */
    j:=Transvectant(Transvectant(t,t,2),t,1);

    /* I2=(theta,theta)_2/M : order 0, degree 2  */
    I2n:=Transvectant(theta,theta,2);
    I2:=I2n[1]/m;

    /* I3=(theta^3,t)_6/M^2 : order 0, degree 3 */
    I3n:=Transvectant([* theta[1]^3, 3*theta[2], 3*theta[3] *], t, 6);
    I3:=I3n[1]/m^2;

    /* I4=(theta^4,(t,t)_2)_8/M^3 : order 0, degree 4 */
    I4n:=Transvectant([* theta[1]^4, 4*theta[2], 4*theta[3] *], Transvectant(t,t,2), 8);
    I4:=I4n[1]/m^3;

    /* J2=((theta,f)_1,(t,t)_2)_8*(theta^6,j)_12/M^6 : order 0, degree 8 : new notation I8 */
    I8n:=Transvectant(Transvectant(theta,f_cov,1),Transvectant(t,t,2),8)[1]
	*Transvectant([* theta[1]^6, 6*theta[2], 6*theta[3] *], j, 12)[1];
    I8:=I8n/m^6;

    /* J3=(36*(theta^2*f,j)_12/M^7+14*((theta^2,f)_3,t)_6/(9*M))*(theta^6,j)_12/M^5 :
       order 0, degree 3 */
    /* WARNING: it seems that there is a mistake and one should replace M^7 by M^2
       so J3 is of degree 9. New notation I9 */
    I9:=(36*Transvectant([* theta[1]^2*f, 2*theta[2]+1, 2*theta[3]+8*],j,12)[1]/m^2+
    14*Transvectant(Transvectant([* theta[1]^2, 2*theta[2], 2*theta[3]*],f_cov,3),t,6)[1]/(9*m))
    *Transvectant([* theta[1]^6, 6*theta[2], 6*theta[3] *],j,12)[1]/m^5;

    /* J4=-2(f*theta^3,t*(t,t)_2)_14/m^3+5*((f,theta^3)_1,j)_12/(21*M^3)+140*((f,theta^3)_4,t)_6/(297*M^2) :
        order 0, degree 4. New notation Ip4 */
    Ip4:=-2*Transvectant([* f*theta[1]^3, 3*theta[2]+1, 3*theta[3]+8 *],
	[* t[1]*Transvectant(t,t,2)[1], 3*t[2], 14 *],14)[1]/m^3+
    5*Transvectant(Transvectant(f_cov,[* theta[1]^3, 3*theta[2], 3*theta[3] *],1),j,12)[1]/(21*m^3)
    +140*Transvectant(Transvectant(f_cov,[* theta[1]^3, 3*theta[2], 3*theta[3] *],4),t,6)[1]/(297*m^2);

    return [ CR!I2 , CR!I3 , CR!I4, CR!Ip4, CR!I8, CR!I9 ];

end intrinsic;

intrinsic MaedaInvariants(C::CrvHyp) -> SeqEnum
    {Compute the Maeda field invariants  'I2', 'I3', 'I4', 'I4p', 'I8' and
    'I9' of a genus 3 hyperelliptic curve C. It works only in characteristic p = 0
    or p >= 11.}

    require Genus(C) eq 3: "Curve must be of genus 3.";
    K := BaseField(C);
    R := PolynomialRing(K); x := R.1;
    f, h := HyperellipticPolynomials(C);
    require (h eq 0) and (Degree(f) le 8):
    "The curve must be of form y^2 = f(x), where f has degree at most 8.";

    return MaedaInvariants(f);
end intrinsic;
