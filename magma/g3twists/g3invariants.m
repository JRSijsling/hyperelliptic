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
 * intrinsic ShiodaInvariants(f::RngUPolElt : normalize := false, scaled := true) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariants(fh::SeqEnum : normalize := false, scaled := true) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariants(C::CrvHyp : normalize := false) -> SeqEnum, SeqEnum
 * intrinsic ShiodaInvariantsEqual(V1::SeqEnum, V2::SeqEnum) -> BoolElt
 * intrinsic ShiodaAlgebraicInvariants(FreeJI::SeqEnum : ratsolve := true) -> SeqEnum
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
    IntegralNormalization := false, scaled := true,
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

    /* Rings of small characteristic  */
    if p eq 3 then
	JI, Wght := ShiodaInvariantsChar3(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 5 then
	JI, Wght := ShiodaInvariantsChar5(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    elif p eq 7 then
	JI, Wght := ShiodaInvariantsChar7(f : PrimaryOnly := PrimaryOnly, degmax := degmax, degmin := degmin);

    else
        JI, Wght := ShiodaInvariantsCharp(f : PrimaryOnly := PrimaryOnly, scaled := scaled, degmax := degmax, degmin := degmin);

    end if;

    if normalize eq false then return JI, Wght; end if;

    w := GCD(Wght);
    return WPSNormalize([e div w : e in Wght], JI), Wght;

end intrinsic;

intrinsic ShiodaInvariants(f::RngUPolElt :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false, scaled := true,
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
        IntegralNormalization := IntegralNormalization, scaled := scaled,
        PrimaryOnly := PrimaryOnly,
        degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(fh::SeqEnum :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false, scaled := true,
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
	return ShiodaInvariantsChar2(f, h :
            normalize := normalize,
            IntegralNormalization := IntegralNormalization, scaled := scaled,
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin);
    end if;

    if h eq 0 then F := f; else F := 4*f+h^2; end if;

    return ShiodaInvariants(F :
            normalize := normalize,
            IntegralNormalization := IntegralNormalization, scaled := scaled,
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin);

end intrinsic;

intrinsic ShiodaInvariants(C::CrvHyp :
    normalize := false,
    PrimaryOnly := false,
    IntegralNormalization := false, scaled := true,
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
            IntegralNormalization := IntegralNormalization, scaled := scaled,
            PrimaryOnly := PrimaryOnly,
            degmax := degmax, degmin := degmin
            );

end intrinsic;


intrinsic ShiodaInvariantsEqual(V1::SeqEnum, V2::SeqEnum) -> BoolElt
    {Check whether Shioda Invariants V1 en V2 of two genus 3 hyperelliptic curves or of
     two binary forms of degree 8 are equivalent.}

    CR := Universe(V1);

    /* Rings of small characteristic  */
    case Characteristic(CR):

    when 2:
	require (#V1 eq 10 and #V2 eq 10) or (#V1 eq 7 and #V2 eq 7)
		or (#V1 eq 6 and #V2 eq 6) or (#V1 eq 5 and #V2 eq 5)
		: "V1, V2 must be of size 5, 6, 7 or 10";

        case #V1:

	when 5:
	    if V1[1] eq 0 and V1[2] eq 0 and V2[1] eq 0 and V2[2] eq 0 then
		return WPSEqual([2, 3, 7, 32, 40], V1, V2);
	    elif V1[1] ne 0 and V1[2] eq 0 and V2[1] ne 0 and V2[2] eq 0 then
		return WPSEqual([2, 3, 2, 2, 2], V1, V2);
	    else
		"Type (3, 3) inconsistent with Type 7";
	    end if;

	when 6:
	    return WPSEqual([2, 3, 1, 3, 4, 5], V1, V2);

	when 7:
	    return WPSEqual([2, 3, 2, 2, 2, 2, 2], V1, V2);

	when 10:
	    return WPSEqual([2, 3, 2, 4, 5, 6, 8, 9, 11, 12], V1, V2);

	end case;

    when 3:
	require #V1 eq 10 and #V2 eq 10 : "V1, V2 must be of size 10 (J2, ..., J10, J12)";

	return WPSEqual([2, 3, 4, 5, 6, 7, 8, 9, 10, 12], V1, V2);

    when 5:
	require #V1 eq 62 and #V2 eq 62 : "V1, V2 must be of size 62";

	Wght := [ 1, 4, 6, 6, 8, 8, 9, 10, 10, 11, 12, 12, 12, 13,
		  14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17,
		  18, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21,
		  22, 22, 23, 23, 24, 24, 24, 25, 25, 25, 26, 27, 27,
		  28, 28, 29, 29, 30, 31, 32, 33, 33, 37 ];
	return WPSEqual(Wght, V1, V2);

    when 7:
	require #V1 eq 13 and #V2 eq 13 : "V1, V2 must be of size 13 (J2, ..., J10, J11, J13, J14, J15)";

	return WPSEqual([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15], V1, V2);

    end case;

    /* Other rings (p = 0 or p > 7) */
    require #V1 eq 9 and #V2 eq 9 : "V1, V2 must be of size 9 (J2, ..., J10)";

    return WPSEqual([2, 3, 4, 5, 6, 7, 8, 9, 10], V1, V2);

end intrinsic;



intrinsic ShiodaAlgebraicInvariants(FreeJI::SeqEnum : ratsolve := true) -> SeqEnum
    {
    This function returns the algebraic relations between the six primary
    invariants given in input and the other invariants.

    In Characteristic 0 or > 7, the primary invariants are J2, J3, J4, J5, J6
    and J7, in characteristic 3, they are J2, J4, J5, J7, J9, J12.

    By default (ratsolve := true), this function computes solutions to
    the system of relations and returns them as a list of invariants. Otherwise
    (ratsolve := false), the system is returned as is.
    }

    FF := Universe(FreeJI);
    p := Characteristic(FF);

    /* Not yet implemented */
    if p in {5} then
	error "[Hyperelliptic] currently, no computation done in fields of char. 5, sorry";
	return [];
    end if;

    /* Rings of small characteristic  */
    case p:

    when 2:
	require #FreeJI eq 6 or #FreeJI eq 5 or #FreeJI eq 4 or #FreeJI eq 3
	    : "Argument must be a sequence of size 3 (N7, N32, N40), 4 (j2, L, L', L''), 5 (j2<>0, K, K', K'', K'''') or (j2=0, M1, M3, M4, M5) or 6 (j2, j3, J5, J6, J8, J9) invariants";

	return ShiodaAlgebraicInvariantsChar2(FreeJI, ratsolve);

    when 3:
	require #FreeJI eq 6:
	    "Argument must be a sequence of six free Shioda invariants J2, J4, J5, J7, J9, J12.";
	return ShiodaAlgebraicInvariantsChar3(FreeJI, ratsolve);

    when 7:
	require #FreeJI eq 6:
	    "Argument must be a sequence of six free Shioda invariants J2, J4, J5, J7, J9, J12.";
	return ShiodaAlgebraicInvariantsChar7(FreeJI, ratsolve);

    end case;

    require #FreeJI eq 6:
	"Argument must be a sequence of six free Shioda invariants J2, J3, J4, J5, J6, J7.";

    if Universe(FreeJI) cmpeq Integers() then
	ChangeUniverse(~FreeJI, Rationals());
    end if;

    return ShiodaAlgebraicInvariantsCharp(FreeJI, ratsolve);

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
