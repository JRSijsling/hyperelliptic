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

### Reconstruct curve models from invariants

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
