Exported intrinsics
--

### Calculate invariants

#### Genus 2

```
intrinsic IgusaInvariants(f::RngUPolElt, h::RngUPolElt :
    extend := false, normalize := false) -> SeqEnum, SeqEnum
intrinsic IgusaInvariants(f::RngUPolElt :
    Quick  := false, extend := false, normalize := false) -> SeqEnum, SeqEnum
intrinsic IgusaInvariants(f::RngMPolElt :
    Quick  := false, extend := false, normalize := false) -> SeqEnum, SeqEnum
intrinsic IgusaInvariants(f::RngUPolElt, p::RngIntElt :
    Quick  := false, extend := false, normalize := false) -> SeqEnum, SeqEnum
intrinsic IgusaInvariants(f::RngMPolElt, p::RngIntElt :
    Quick  := false, extend := false, normalize := false) -> SeqEnum, SeqEnum
intrinsic IgusaInvariants(C::CrvHyp :
    Quick := false, extend := false, normalize := false) -> SeqEnum, SeqEnum

intrinsic IgusaInvariantsEqual(V1::SeqEnum, V2::SeqEnum) -> BoolElt
intrinsic DiscriminantFromIgusaInvariants(JI::SeqEnum) -> .
intrinsic IgusaAlgebraicRelations(JI::SeqEnum) -> SeqEnum

intrinsic G2Invariants(H::CrvHyp) -> SeqEnum
intrinsic IgusaToG2Invariants(JI::SeqEnum) -> SeqEnum
intrinsic G2ToIgusaInvariants(GI::SeqEnum) -> SeqEnum
```

#### Genus 3

```
intrinsic ShiodaInvariants(f::RngUPolElt, p::RngIntElt :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
intrinsic ShiodaInvariants(f::RngMPolElt, p::RngIntElt :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
intrinsic ShiodaInvariants(f::RngUPolElt :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
intrinsic ShiodaInvariants(f::RngMPolElt :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
intrinsic ShiodaInvariants(fh::SeqEnum :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum
intrinsic ShiodaInvariants(C::CrvHyp :
    normalize := false, PrimaryOnly := false, IntegralNormalization := false, degmax := Infinity(), degmin := 1) -> SeqEnum, SeqEnum

intrinsic ShiodaInvariantsEqual(V1::SeqEnum, V2::SeqEnum) -> BoolElt
intrinsic DiscriminantFromShiodaInvariants(JI::SeqEnum) -> .
intrinsic ShiodaAlgebraicInvariants(PrimaryInvariants::SeqEnum :
    ratsolve := true) -> SeqEnum

intrinsic MaedaInvariants(f::RngUPolElt) -> SeqEnum
intrinsic MaedaInvariants(C::CrvHyp) -> SeqEnum
```

### Reconstruct curve from invariants

#### Genus 2

```
intrinsic HyperellipticCurveFromIgusaInvariants(II::SeqEnum :
    RationalModel := false) -> CrvHyp, GrpPerm
intrinsic HyperellipticCurveFromG2Invariants(GI::SeqEnum :
    RationalModel := false) -> CrvHyp, GrpPerm
```

#### Genus 3

```
intrinsic HyperellipticCurveFromShiodaInvariants(JI::SeqEnum :
    RationalModel := true) -> CrvHyp, GrpPerm
intrinsic HyperellipticPolynomialFromShiodaInvariants(JI::SeqEnum :
    RationalModel := true) -> SeqEnum, GrpPerm
intrinsic HyperellipticPolynomialsFromShiodaInvariants(JI::SeqEnum :
    RationalModel := true) -> SeqEnum, GrpPerm
```

### Compute twists

#### Genus 2

```
intrinsic TwistsFromIgusaInvariants(JI::SeqEnum[FldFinElt]) -> SeqEnum[CrvHyp], GrpPerm
intrinsic TwistsFromG2Invariants(II::SeqEnum[FldFinElt]) -> SeqEnum[CrvHyp], GrpPerm
```

#### Genus 3

```
intrinsic TwistsFromShiodaInvariants(JI::SeqEnum[FldFinElt]) -> SeqEnum[CrvHyp], GrpPerm
intrinsic TwistedHyperellipticPolynomialsFromShiodaInvariants(JI::SeqEnum[FldFinElt] :
    RationalModel := true, Deterministic := false) -> SeqEnum, GrpPerm
```

#### Any genus

```
intrinsic Twists(H::CrvHyp :
    AutomorphismGroup := false) -> SeqEnum[CrvHyp], GrpPerm

intrinsic TwistsOfHyperellipticPolynomials(f::RngUPolElt, d::RngIntElt :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm
intrinsic TwistsOfHyperellipticPolynomials(f::RngUPolElt :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm
intrinsic TwistsOfHyperellipticPolynomials(fh::SeqEnum[RngUPolElt] :
    AutomorphismGroup := false) -> SeqEnum, GrpPerm

intrinsic TwistsOverFiniteField(C::Crv, Aut::.) -> SeqEnum[Crv]
```

### Isomorphisms

#### Genus 1

```
intrinsic GeometricAutomorphismGroup(Ec::CrvEll) -> GrpPerm
```

#### Genus 2

```
intrinsic GeometricAutomorphismGroupFromIgusaInvariants(II::SeqEnum) -> GrpPerm
intrinsic GeometricAutomorphismGroupGenus2Classification(FF::FldFin) -> SeqEnum, SeqEnum
```

#### Genus 3

```
intrinsic GeometricAutomorphismGroupFromShiodaInvariants(JI::SeqEnum) -> GrpPerm
intrinsic GeometricAutomorphismGroupGenus3Classification(FF::FldFin) -> SeqEnum, SeqEnum
```

#### Any genus

```
intrinsic IsIsomorphicHyperelliptic(X1::CrvHyp, X2::CrvHyp :
    geometric := false, covariant := true, commonfield := false) ->  BoolElt, List
intrinsic IsGL2EquivalentExtended(f1::RngUPolElt, f2::RngUPolElt, deg::RngIntElt :
    geometric := false, covariant := true, commonfield := false) -> BoolElt, List

intrinsic AutomorphismGroupHyperellipticPolynomial(f::RngUPolElt, deg::RngIntElt :
    geometric := false, covariant := true, commonfield := false) -> List
intrinsic AutomorphismGroupHyperelliptic(X::CrvHyp :
    geometric := false, covariant := true, commonfield := false) -> List

intrinsic GeometricAutomorphismGroup(H::CrvHyp) -> GrpPerm
```

### Weighted projective space toolbox

```
intrinsic WPSEnumInit(FF::FldFin, Wght::SeqEnum) -> Rec
intrinsic WPSEnumNext(~V, ~WPSCtxt::Rec)
intrinsic WPSEqual(Wght::SeqEnum, V1::SeqEnum, V2::SeqEnum) -> BoolElt
intrinsic WPSNormalize(Wght::SeqEnum, V::SeqEnum) -> SeqEnuminsic;
intrinsic WPSMultiply(Wght::SeqEnum, V::SeqEnum, lambda::.) -> SeqEnum
intrinsic WPSMinimize(Wght::SeqEnum, V::SeqEnum) -> SeqEnum
```
