//freeze;

/***
 *  Twists of curves.
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
 *  Copyright 2020, R. Lercier, C. Ritzenthaler & J. Sijsling
 */

 /***
  *
  * Given a finite field Fq and a curve C, either hyperelliptic or a plane smooth quartic
  * Twists(C) compute a set of representatives of the twists of C. This works
  * - for C hyperelliptic g=2 and 3 without restriction on q
  * - For C hyperelliptic g>3 with q odd
  * - For C plane smooth quartic and q a power of a prime>7.
  *
  * Note that the function TwistsOverFiniteField can also be useful if you already know
  * a subgroup H of the group of all automorphisms of C given as a list of
  * linear projective transformations (i.e. matrices) of the ambient projective space. In
  * that case you can compute twists with respect to these elements for any non-hyperelliptic
  * curve (not necessarily plane). Caution: all matrices in H
  * have been normalized so that they have a 1 as first non-zero coefficient runnning
  * the matrix in column.
  *
  *********************************************************************/

import "../g2twists/g2twists.m" : G2Models;
import "../g3twists/g3twists.m" : G3Models;

/* return f.M when f is seen as a form of degree n
   (M lives in an extension of the field of definition of f)
*/

function TransformForm(A,f,n)
  F2 := BaseRing(Parent(A));
  F := BaseRing(Parent(f));
  Embed(F,F2);
  f2 := ChangeRing(f,F2);
  P2 := PolynomialRing(F2,2); x := P2.1; z := P2.2;
  P := PolynomialRing(F2); u := P.1;
  phi := hom<P2 -> P | u,1>;
  f3 := P2!(Evaluate(f2,x/z)*z^n);
  f4 := Evaluate(f3,[A[1,1]*x+A[1,2]*z,A[2,1]*x+A[2,2]*z]);
  co := Coefficients(f4);
  return phi(f4)/co[1];
end function;


/* This function return a list [*m, lamda*] such that
   FrobeniusImage(M,e)^m-1 .. FrobeniusImage(M,e)M = lambda Id
   Where M is a NXN matrix defined over the finite field Fs, #Fs = q^r.
   F is the field other which my curve is defined. #F = p^e = q
   FrobeniusImage(M,e)^m = FrobeniusImage(..FrobeniusImage(M,e)..,e) m times.
*/
function FrobeniusOrderAutomorphism(M, F)

    e := Degree(F);
    P := M;
    m := 1;
    while IsScalar(P) eq false or P[1,1] notin F do
        P := FrobeniusImage(P,e)*M;
        m := m+1;
    end while;

    return [*m, P[1, 1]*];
end function;


/* This function return a Matrix with a 1 which is a multiple of M
*/

function NormalizedM(M)

    for j := 1 to Nrows(M) do
        for i := 1 to Nrows(M) do
            if M[i,j] ne 0 then return (1 / M[i,j]) * M; end if;
        end for;
    end for;

    return M;
end function;

/* check if the projective matrix M is defined over F
*/

function IsDefinedOver(M,F)

    if M eq 0 then return true, ChangeRing(M,F); end if;
    e:=Degree(F);
    if NormalizedM(M) eq FrobeniusImage(NormalizedM(M),e) then
        return true, ChangeRing(NormalizedM(M),F);
    end if;

    return false, M;
end function;

/* smallest field of definition
*/

function SmallestField(M,F)
    Fs := BaseRing(M);
    r := Degree(Fs,F);
    for d in Divisors(r) do
        Fd := ext<F | d>;
        if IsDefinedOver(M, Fd) then
            return ChangeRing(NormalizedM(M), Fd);
        end if;
    end for;

end function;

/* Given the (reduced) automorphism group H and the finite field F,
   this function compute it cohomology classes.
*/
function CohomologyClass(H, F)
    L := H;
    e := Degree(F);
    CohoClass := [* Identity(Parent(H[1])) *];
    while not IsEmpty(L) do
        for h in H do
            EqClassCoho := h^(-1)*CohoClass[#CohoClass]*FrobeniusImage(h,e);
            Exclude(~L, NormalizedM(EqClassCoho));
        end for;
        if IsEmpty(L) eq false then
            Append(~CohoClass, SmallestField(L[1],F));
        end if;
    end while;
    return CohoClass;
end function;



/* Given a finite field F of size q = p^e and a matrix M defined over an extension of F,
   this function compute an invertible matrix A such that FrobeniusImage(A,e)^(-1)*A=Mb.
   Where Mb = 1/X*M verify FrobeniusImage(Mb,e)^m-1 .. FrobeniusImage(Mb,e)Mb = Id
*/

function ComputationCobord(M, F)

    Fs := BaseRing(M);
    N := Nrows(M);
    r := Degree(Fs);
    e := Degree(F);
    L := FrobeniusOrderAutomorphism (M,F);
    m := L[1] ;
    K := ext<F | m>;
    _, X := NormEquation(K,F!L[2]);
    //  n:=Degree(Fs,F);
    //  K2:=ext<Fs | Lcm(m,n)>;

    Mb := ChangeRing(M,K)*1/X;
    // check S:=Mb;for i:=1 to m-1 do S:=FrobeniusImage(S,e)*Mb; end for;S;

    repeat
        Pb := RandomMatrix(K,N,N);
        A := Pb;
        for i := 1 to m-1 do
            Pb := FrobeniusImage(Pb,e)*Mb;
            A :=  A + Pb;
        end for;
    until IsUnit(A);

    // in this step we have computed A such that FrobeniusImage(A,e)*Mb-A =0
    // in this step we have computed A such that A*FrobeniusImage(A,e)^(-1)-Mb =0

    return A;
end function;

/* Given the reduced Automorphism group G of the hyperelliptic curve C/F,
   the cocycle 2x2 matrix A from above,
   this function compute the reduced automorphism  group defined over F of
   the twist associated to A.
*/

function ReducedAutomorphismGroupOfTwistDefinedOverF(Aut, f, A, g)

    F := BaseRing(f);
    K := BaseRing(A);
    P := Parent(f); x := P.1;
    Embed(F, K);
    L := []; // List of the reduced Automorphism on the twist
    eN := []; // List of the associate constant to find the Automorphisms
    ft:=P!TransformForm(A,f,2*g+2);
    for M in Aut do
        B := A^(-1)*M*A;
        boo, B := IsDefinedOver(B,F);
        if boo then
            Append(~L,B);
            // M; Evaluate(ft,(B[1,1]*x+B[1,2])/(B[2,1]*x+B[2,2]))*(B[2,1]*x+B[2,2])^(2*g+2)/ft;
            Append (~eN,F!(Evaluate(ft,(B[1,1]*x+B[1,2])/(B[2,1]*x+B[2,2]))*(B[2,1]*x+B[2,2])^(2*g+2)/ft));
        end if;
    end for;

    return [* L, eN *];
end function;

/* test of self duality for the twist of C  given by the cocycle A
*/

function IsSelfDual(Aut,f,A,g)
    F := BaseRing(f);
    RAGT := ReducedAutomorphismGroupOfTwistDefinedOverF(Aut,f,A,g);
    return 2*#[e : e in RAGT[2] | IsSquare(e)] eq #RAGT[2];
end function;



function ProjectiveMatrixGroup(L)

    n := Nrows(L[1]);
    FF := BaseRing(L[1]);
    prim := PrimitiveElement(FF);
    MM := MatrixAlgebra(FF,n);
    H :=  MatrixGroup< n, FF | L cat [prim*Identity(MM)]>;
    C :=  MatrixGroup< n, FF | [prim*Identity(MM)]>;
    _, I, _ := CosetAction(H,C);

    //phi, I, K := CosetAction(H,C);
    //psi := Inverse(phi);
    // return [MM!psi(h) : h in I], I;

    return I;
end function;

function Normalize22Column(T)
    col := Eltseq(Rows(Transpose(T))[1]);
    i0 := Minimum([ i : i in [1..#col] | col[i] ne 0 ]);
    return (1/col[i0]) * T;
end function;

/* Given a projective curve C and its Automorphism group Aut
   (given by a list of matrix N x N), this function compute the twists of C
*/

intrinsic TwistsOverFiniteField(C::Crv, Aut::.) -> SeqEnum[Crv]
    {Determines the twists of C using the known automorphisms in Aut.}

    F := BaseRing(C);

    require Type(F) eq FldFin :
	"Twist computations only available in finite fields";

    g := Genus(C);

    ishyper, H :=  IsHyperelliptic(C);
    if ishyper then
        t := PrimitiveElement(F);
        f := HyperellipticPolynomials(H);
        P := Parent(f);
    end if;

    Coh := CohomologyClass(Aut, F);
    T := [C];

    if ishyper then
        boo := IsSelfDual(Aut,f,Identity(Parent(Coh[1])), g);
        if boo eq false then
            Append(~T,HyperellipticCurve(t*f));
        end if;
    end if;

    for M in Coh do
        if M ne Identity(Parent(M)) then
            A := ComputationCobord(M^(-1),F);
            F2 := BaseRing(A);
            C2 := ChangeRing(C,F2);

            if ishyper then
                ft := P!TransformForm(A^(-1), f, 2*g+2);
                boo := IsSelfDual(Aut, f, A^(-1), g);
                Append(~T, HyperellipticCurve(ft));
                if boo eq false then
                    Append(~T,HyperellipticCurve(t*ft));
                end if;
            else
                X2 := Automorphism(AmbientSpace(C2), A^(-1));
                Append(~T, ChangeRing(X2(C2), F));
            end if;

        end if;
    end for;

    return T;
end intrinsic;

intrinsic Twists(H::CrvHyp : AutomorphismGroup := false) -> SeqEnum[CrvHyp], GrpPerm
    {Compute twisted  hyperelliptic curves and their automorphism groups}

    F := CoefficientRing(H);

    require Type(F) eq FldFin :
	"Twist computations only available in finite fields";

    g := Genus(H);

    if g eq 2 then
        twists, aut := G2Models(IgusaInvariants(H));
        if AutomorphismGroup then return twists, aut; end if;
        return twists;
    end if;

    if g eq 3 then
        twists, aut := G3Models(ShiodaInvariants(H));
        if AutomorphismGroup then return twists, aut; end if;
        return twists;
    end if;

    require Characteristic(F) ge 3 :
        "2 must be invertible in the base ring.";

    f, h := HyperellipticPolynomials(H);
    g := 4*f + h^2;
    d := 2*((Degree(g) + 1) div 2);

    _, Aut := IsGL2EquivalentExtended(f, f, d : geometric := true, commonfield := true);

    twists := TwistsOverFiniteField(HyperellipticCurve(g), [ Normalize22Column(A) : A in Aut ]);

    if AutomorphismGroup then return twists, ProjectiveMatrixGroup(Aut); end if;
    return twists;

end intrinsic;

intrinsic HyperellipticPolynomialTwists(f::RngUPolElt, n::RngIntElt) -> SeqEnum[RngUPolElt]
    {Found polynomials fp s.t. the curves y^2 = f(x) and y^2 = fp(x) are
    twisted each other.}

    _, Aut := IsGL2EquivalentExtended(f, f, n : geometric := true, commonfield := true);

    Twists := TwistsOverFiniteField(HyperellipticCurve(f), [ Normalize22Column(A) : A in Aut ]);

    return [HyperellipticPolynomials(g) : g in Twists];

end intrinsic;
