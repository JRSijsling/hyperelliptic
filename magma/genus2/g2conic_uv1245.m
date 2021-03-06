//freeze;

/***
 *  Parameterized equations of the conic and the quartic for the covariants [1, 2, 4-5].
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
 *  Copyright 2020, R. Lercier & C. Ritzenthaler & J. Sijsling
 */

function Genus2ConicAndCubicUV1245(uv, JI : models := true)

    FF := Parent((Universe(JI)!1)/1); /* We want a field */

    J2, J4, J6, _, J10, J15 := Explode(JI);

    u, v := Explode(uv);

    P3 := PolynomialRing(Universe(uv), 3); w1 := P3.1; w2 := P3.2; w3 := P3.3;

    R :=
    (188160*u-75600*v)*J2^4*J10+(-10035200*u+4416000*v)*J2^2*J4*J10+(-225792000*u-36000000*v)*J2*J6*J10-61440000*v*J4^2*J10+(2352*u+135*v)*J2^6*J6+(784*u+45*v)*J2^5*J4^2+(-134848*u-17340*v)*J2^4*J4*J6+(-50176*u-2880*v)*J2^3*J4^3+(592704*u+5220*v)*J2^3*J6^2+(1806336*u+650880*v)*J2^2*J4^2*J6+(802816*u+46080*v)*J2*J4^4+(-13547520*u-1814400*v)*J2*J4*J6^2-6451200*v*J4^3*J6+15552000*v*J6^3;
    if R eq 0 then
        return R, R, R;
    end if;

    if not models then return R, R, R; end if;

    C :=
        153837509765625*(3*J2^3-160*J2*J4-3600*J6)*w1^2+
        9845600625000*(3*J2^4-180*J4*J2^2+3000*J2*J6+3200*J4^2)*w1*w2+
        157529610000*(3*J2^5-200*J4*J2^3-400*J6*J2^2+3200*J4^2*J2-48000*J4*J6-800000*J10)*w2^2+
        186046875*v*J15*w2*w3+
        4*(-3686400000000*v^2*J10^2-240*(784*u+45*v)*(784*u-675*v)*J2^5*J10+12800*(614656*u^2-540960*u*v-43875*v^2)*J2^3*J4*J10+288000*(614656*u^2+196000*u*v+6025*v^2)*J2^2*J6*J10+860160000*v*(15*v+112*u)*J2*J4^2*J10-608256000000*v^2*J10*J6*J4-3*(784*u+45*v)^2*J2^7*J6-(784*u+45*v)^2*J2^6*J4^2+4*(33712*u+6735*v)*(784*u+45*v)*J2^5*J4*J6+64*(784*u+45*v)^2*J2^4*J4^3-36*(12907776*u^2+227360*u*v+130525*v^2)*J2^4*J6^2-9216*(153664*u^2+110740*u*v+7725*v^2)*J2^3*J4^2*J6-1024*(784*u+45*v)^2*J2^2*J4^4+120960*(87808*u^2+23520*u*v+2775*v^2)*J2^2*J4*J6^2+1843200*v*(815*v+5488*u)*J2*J4^3*J6-217728000*v*(-25*v+112*u)*J2*J6^3-19906560000*v^2*J6^2*J4^2)*w3^2;

    M :=
        195385944403125000000000*(J2^5-100*J4*J2^3+600*J6*J2^2+3200*J4^2*J2-288000*J4*J6-3200000*J10)*w1^3+
        6252350220900000000000*(3*J2^6-320*J4*J2^4-1600*J6*J2^3+10400*J4^2*J2^2+432000*J6*J4*J2-128000*J4^3+4000000*J2*J10+2160000*J6^2)*w1^2*w2+
        184605011718750000*J15*J2*(-784*u-25*v)*w1^2*w3+
        200075207068800000000*(3*J2^7-340*J4*J2^5+13000*J4^2*J2^3-180000*J6*J4*J2^2-160000*J4^3*J2-2400000*J2^2*J10-2700000*J6^2*J2+32000000*J10*J4)*w1*w2^2+
        4922800312500000*((784*u+69*v)*J15*J2^2-1920*v*J4*J15)*w1*w2*w3+
        2540160000*(-224*(784*u+45*v)^2*J2^6*J4^3+12*(124775168*u^2-96063520*u*v-5444925*v^2)*J2^6*J6^2+11264*(784*u+45*v)^2*J2^4*J4^4+648000*(614656*u^2-556640*u*v+6665*v^2)*J2^3*J6^3-163840*(784*u+45*v)^2*J2^2*J4^5-1536000000*(614656*u^2+39200*u*v+5025*v^2)*J2^2*J10^2+3*(784*u+45*v)^2*J2^9*J6+(784*u+45*v)^2*J2^8*J4^2+165888000000*v*(75*v+784*u)*J2*J6^2*J10+12*(784*u+45*v)*(46256*u+1055*v)*J2^7*J4*J6-147456000*v*(745*v+8624*u)*J2*J4^4*J6+12441600000*v*(11*v+784*u)*J2*J4*J6^3-69120000*(1843968*u^2+54880*u*v+33875*v^2)*J2^2*J4*J6*J10-9830400000*v*(75*v+784*u)*J2*J4^3*J10-16*(2108884736*u^2+277841760*u*v+7919775*v^2)*J2^5*J4^2*J6+2880*(20283648*u^2+26538400*u*v+2424425*v^2)*J2^4*J4*J6^2+53760*(10624768*u^2+2621920*u*v+119775*v^2)*J2^3*J4^3*J6-5529600*(614656*u^2+235200*u*v+46875*v^2)*J2^2*J4^2*J6^2+960*(784*u+45*v)*(10192*u+405*v)*J2^7*J10-48000*(10449152*u^2+1224608*u*v+32985*v^2)*J2^5*J4*J10+96000*(16595712*u^2+525280*u*v+229875*v^2)*J2^4*J6*J10+4608000*(1843968*u^2+305760*u*v+11875*v^2)*J2^3*J4^2*J10+3450470400000*v^2*J4^3*J6^2+294912000000000*v^2*J4*J10^2-4478976000000*v^2*J6^4+66355200000000*v^2*J4^2*J6*J10)*w1*w3^2+
        6402406626201600000*(J2^8-120*J2^6*J4+2200*J2^5*J6+5800*J2^4*J4^2-100000*J2^3*J4*J6-136000*J2^2*J4^3+400000*J2^3*J10+1140000*J2^2*J6^2+1760000*J2*J4^2*J6+1280000*J4^4-8000000*J2*J4*J10-7200000*J4*J6^2-80000000*J6*J10)*w2^3+
        15752961000000*(-3*(784*u+45*v)*J15*J2^3+280*(15*v+112*u)*J15*J2*J4+72000*v*J6*J15)*w2^2*w3+
        81285120*(3*(784*u+45*v)^2*J2^10*J6+(784*u+45*v)^2*J2^9*J4^2-64*(784*u+45*v)^2*J2^7*J4^3-4*(301796096*u^2-104593440*u*v-8437725*v^2)*J2^7*J6^2+1024*(784*u+45*v)^2*J2^5*J4^4-86400*(4302592*u^2-1230880*u*v-69625*v^2)*J2^4*J6^3+640000000*(614656*u^2+108192*u*v+1305*v^2)*J2^3*J10^2-5529600000*v*(145*v+784*u)*J2*J4^3*J6^2+120960000*(614656*u^2+70560*u*v+3625*v^2)*J2^3*J4*J6*J10-2457600000000*v*(15*v+784*u)*J2*J4*J10^2+143360000*(784*u+825*v)*(15*v+112*u)*J2^2*J4^3*J10-4*(159152*u+13935*v)*(784*u+45*v)*J2^8*J4*J6-18662400000*v*(-275*v+784*u)*J2*J6^4-240*(18032*u+1755*v)*(784*u+45*v)*J2^8*J10+1600*(143214848*u^2+30364320*u*v+1357425*v^2)*J2^6*J4*J10-16000*(69456128*u^2+1952160*u*v+516825*v^2)*J2^5*J6*J10-4480000*(965888*u^2+333984*u*v+20925*v^2)*J2^4*J4^2*J10+172800000*(614656*u^2-117600*u*v-32775*v^2)*J2^2*J6^2*J10+32*(852527872*u^2+181127520*u*v+8127675*v^2)*J2^6*J4^2*J6-9600*(614656*u^2+2654624*u*v+363945*v^2)*J2^5*J4*J6^2-640*(705010432*u^2+273137760*u*v+17674875*v^2)*J2^4*J4^3*J6+115200*(12907776*u^2+4257120*u*v+934625*v^2)*J2^3*J4^2*J6^2+1536000*(614656*u^2+1168160*u*v+143025*v^2)*J2^2*J4^4*J6+5184000*(614656*u^2-148960*u*v-87375*v^2)*J2^2*J4*J6^3-1474560000000*v^2*J4^5*J6-4976640000000*v^2*J4^2*J6^3-11796480000000*v^2*J4^4*J10-2211840000000000*v^2*J6*J10^2-110592000000*v*(65*v+2352*u)*J2*J4^2*J6*J10-298598400000000*v^2*J4*J6^2*J10)*w2*w3^2+
        ((784*u+45*v)^3*J15*J2^8-32*(1568*u+165*v)*(784*u+45*v)^2*J15*J2^6*J4+384*(784*u+45*v)*(307328*u^2+84280*u*v-11175*v^2)*J15*J2^5*J6+1024*(195*v+784*u)*(784*u+45*v)^2*J15*J2^4*J4^2-5760*(481890304*u^3+193616640*u^2*v-1352400*u*v^2-3024375*v^3)*J15*J2^3*J4*J6-128000*(-15*v+784*u)*(614656*u^2+70560*u*v+16425*v^2)*J15*J2^3*J10-2457600*v*(784*u+45*v)^2*J15*J2^2*J4^3+10368000*v*(784*u+5*v)*(784*u-155*v)*J15*J2^2*J6^2+1105920000*v^2*(-365*v+2352*u)*J15*J2*J4^2*J6+73728000000*v^2*(-15*v+784*u)*J15*J2*J4*J10+1990656000000*v^3*J4*J6^2*J15+44236800000000*v^3*J6*J10*J15)*w3^3;

    return R, C, M;

end function;
