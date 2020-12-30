//freeze;

/***
 *  Parameterized equations of the conic and the quartic for the covariants [1, 2, 4-6].
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

function Genus2ConicAndCubicUV1246(uv, JI : models := true)

    FF := Parent((Universe(JI)!1)/1); /* We want a field */

    J2, J4, J6, _, J10, J15 := Explode(JI);

    u, v := Explode(uv);

    P3 := PolynomialRing(Universe(uv), 3); w1 := P3.1; w2 := P3.2; w3 := P3.3;

    /* Conic [1, 2, 4*J4*u + 6*v ] */
    R :=
        -1536000000*J10^2*v+960*(13*v+u)*J2^5*J10-54400*(2*u+15*v)*J4*J2^3*J10-
        288000*(-9*v+4*u)*J6*J2^2*J10+1536000*(9*v+2*u)*J4^2*J2*J10+69120000*(u-3*v)*J4*J6*J10+
        3*(4*u+v)*J6*J2^7+(4*u+v)*J4^2*J2^6-4*(352*u-177*v)*J6*J4*J2^5-16*(31*u+14*v)*J4^3*J2^4+
        84*(36*u+29*v)*J6^2*J2^4+16*(3156*u-3431*v)*J6*J4^2*J2^3+1024*(19*u+11*v)*J4^4*J2^2-
        8640*(-11*v+29*u)*J6^2*J4*J2^2-7680*(72*u-121*v)*J6*J4^3*J2+648000*J6^3*J2*v-
        81920*(3*u+2*v)*J4^5+1382400*(3*u-4*v)*J6^2*J4^2
        ;
    if R eq 0 then
        return R, R, R;
    end if;

    if not models then return R, R, R; end if;

    C :=
        158939263916015625*(3*J2^3-160*J2*J4-3600*J6)*w1^2+
        10172112890625000*(3*J2^4-180*J4*J2^2+3000*J2*J6+3200*J4^2)*w1*w2+
        299003906250*v*J15*w1*w3+
        162753806250000*(3*J2^5-200*J4*J2^3-400*J6*J2^2+3200*J4^2*J2-48000*J4*J6-800000*J10)*w2^2-
        3986718750*v*J2*J15*w2*w3+
        1024*(-43200*v*(611*v+120*u)*J6^3*J2^3+81920*(3*u+2*v)*(23*u+12*v)*J4^5*J2^2-122880000000*v*(v+6*u)*J4*J10^2+32*(23*u+12*v)*(4*u+v)*J4^3*J2^6-240*(4*u+103*v)*(4*u+v)*J2^7*J10+256000000*v*(11*v+48*u)*J2^2*J10^2-11520000*(48*u^2-252*u*v-47*v^2)*J2^2*J4*J6*J10+4*(4*u+v)*(532*u-397*v)*J6*J4*J2^7+62208000*v*(13*v+5*u)*J6^3*J4*J2-4*(3024*u^2+4872*u*v+29929*v^2)*J2^6*J6^2-256*(769*u^2+772*u*v+184*v^2)*J2^4*J4^4+110592000*(9*u^2-24*u*v+5*v^2)*J4^3*J6^2-1866240000*v^2*J6^4-(4*u+v)^2*J4^2*J2^8-3*(4*u+v)^2*J6*J2^9+13824000000*J2*J6^2*J10*v^2-6553600*(3*u+2*v)^2*J4^6+3200*(208*u^2+3912*u*v+863*v^2)*J2^5*J4*J10+32000*(144*u^2-648*u*v-137*v^2)*J2^4*J6*J10-1536000*(25*u^2+327*u*v+68*v^2)*J2^3*J4^2*J10+81920000*(9*u^2+81*u*v+17*v^2)*J2*J4^3*J10+5529600000*(3*u^2-18*u*v-v^2)*J10*J4^2*J6-16*(33744*u^2-48688*u*v-7731*v^2)*J2^5*J4^2*J6+2880*(600*u^2+142*u*v+2693*v^2)*J2^4*J4*J6^2+2560*(5598*u^2-13197*u*v-1096*v^2)*J2^3*J4^3*J6-691200*(111*u^2-130*u*v+207*v^2)*J2^2*J4^2*J6^2-3686400*(36*u^2-121*u*v-5*v^2)*J2*J4^4*J6)*w3^2;

    M :=
        50691691485214233398437500*(J2^5-100*J4*J2^3+600*J6*J2^2+3200*J4^2*J2-288000*J4*J6-3200000*J10)*w1^3+
        1622134127526855468750000*(3*J2^6-320*J4*J2^4-1600*J6*J2^3+10400*J4^2*J2^2+432000*J6*J4*J2-128000*J4^3+4000000*J2*J10+2160000*J6^2)*w1^2*w2+
        3973481597900390625*(-(24*u+11*v)*J15*J2^2+480*(3*u+v)*J15*J4)*w1^2*w3+
        51908292080859375000000*(3*J2^7-340*J4*J2^5+13000*J4^2*J2^3-180000*J6*J4*J2^2-160000*J4^3*J2-2400000*J2^2*J10-2700000*J6^2*J2+32000000*J10*J4)*w1*w2^2+
        508605644531250000*((5*u+3*v)*J15*J2^3-20*(15*u+8*v)*J15*J2*J4-3600*v*J6*J15)*w1*w2*w3+
        163296000000*(3*(4*u+v)^2*J2^11*J6+(4*u+v)^2*J4^2*J2^10-1048576000*(3*u+2*v)^2*J4^7+3276800000*(243*u^2+204*u*v+85*v^2)*J2*J4^4*J10-147456000000*(81*u^2+18*u*v+25*v^2)*J4^3*J6*J10-64*(32264*u^2+27372*u*v+3189*v^2)*J2^7*J4^2*J6-640*(4932*u^2+186735*u*v+703*v^2)*J2^6*J4*J6^2+640*(251640*u^2+127844*u*v+15521*v^2)*J2^5*J4^3*J6-230400*(567*u^2-21126*u*v+4885*v^2)*J2^4*J4^2*J6^2-102400*(48303*u^2+14592*u*v+3550*v^2)*J2^3*J4^4*J6-1728000*(720*u^2+11076*u*v-2755*v^2)*J2^3*J4*J6^3+55296000*(291*u^2-1322*u*v+520*v^2)*J2^2*J4^3*J6^2+147456000*(363*u^2+60*u*v+37*v^2)*J2*J4^5*J6+2488320000*(15*u^2+156*u*v-79*v^2)*J2*J4^2*J6^3-12800*(2892*u^2+2653*u*v+457*v^2)*J2^7*J4*J10+64000*(648*u^2+548*u*v+7449*v^2)*J2^6*J6*J10+128000*(19584*u^2+16884*u*v+3977*v^2)*J2^5*J4^2*J10-20480000*(3591*u^2+2992*u*v+974*v^2)*J2^3*J4^3*J10-115200000*v*(-539*v+960*u)*J2^3*J6^2*J10+61440000000*(2*u+v)*(24*u+7*v)*J2^2*J4*J10^2+4*(348*u+617*v)*(4*u+v)*J2^9*J4*J6+480*(104*u+77*v)*(4*u+v)*J2^9*J10-32*(43*u+17*v)*(4*u+v)*J4^3*J2^8+6553600*(3*u+2*v)*(49*u+26*v)*J4^6*J2^2+933120000*v*(-3*v+16*u)*J2^2*J6^4-298598400000*v*(-2*v+3*u)*J4*J6^4+110592000000*v*(-53*v+60*u)*J2*J4*J6^2*J10+9953280000000*v^2*J6^3*J10-17694720000*(18*u^2-15*u*v+7*v^2)*J4^4*J6^2-512000000*(48*u^2+44*u*v+9*v^2)*J2^4*J10^2-29491200000000*(3*u^2+2*u*v+v^2)*J4^2*J10^2+4*(9744*u^2+239432*u*v+88149*v^2)*J2^8*J6^2+256*(2609*u^2+2192*u*v+424*v^2)*J4^4*J2^6+172800*(60*u^2+1222*u*v+267*v^2)*J2^5*J6^3-40960*(907*u^2+936*u*v+232*v^2)*J4^5*J2^4-55296000000000*v^2*J2*J6*J10^2-2560000*(3240*u^2+2514*u*v+12809*v^2)*J2^4*J4*J6*J10+921600000*(594*u^2+330*u*v+709*v^2)*J2^2*J4^2*J6*J10)*w1*w3^2+
        1661065346587500000000*(J2^8-120*J2^6*J4+2200*J2^5*J6+5800*J2^4*J4^2-100000*J2^3*J4*J6-136000*J2^2*J4^3+400000*J2^3*J10+1140000*J2^2*J6^2+1760000*J2*J4^2*J6+1280000*J4^4-8000000*J2*J4*J10-7200000*J4*J6^2-80000000*J6*J10)*w2^3+
        2034422578125000*(-3*(4*u-7*v)*J15*J2^4+20*(44*u-63*v)*J15*J4*J2^2+21000*v*J2*J6*J15-3200*(3*u-7*v)*J15*J4^2)*w2^2*w3+
        1275750*(209715200*(5717280*u^2-3816646*u*v+161031*v^2)*J6*J4*J2^5*J10-838860800*(95482260*u^2-64933272*u*v-1000583*v^2)*J6*J4^2*J2^3*J10+226492416000*(929280*u^2-682516*u*v+50951*v^2)*J6^2*J4*J2^2*J10+402653184000*(3500235*u^2-2404392*u*v-75013*v^2)*J6*J4^3*J2*J10+26843545600*(3*u+2*v)^2*J4^6*J2^3+12288*(4*u+v)^2*J6*J2^12+4096*(4*u+v)^2*J4^2*J2^11-3669177139200*(51960*u^2-36362*u*v-343*v^2)*J6^5-16106127360000000*(51960*u^2-36362*u*v-1593*v^2)*J10^3+75*(51960*u^2-36362*u*v+2157*v^2)*J15^2-1258291200*(1783560*u^2-1336343*u*v+49628*v^2)*J6^2*J4^3*J2^3+10066329600*(85515*u^2-55264*u*v-3116*v^2)*J6*J4^5*J2^2-169869312000*(105000*u^2-68674*u*v-5141*v^2)*J6^3*J4^2*J2^2+362387865600*(105495*u^2-73999*u*v-986*v^2)*J6^2*J4^4*J2+1223059046400*(103920*u^2-70849*u*v-2561*v^2)*J6^4*J4*J2+6553600*(10352*u^2-6088*u*v-1965*v^2)*J4*J2^8*J10-65536000*(1808*u^2+3144*u*v+14613*v^2)*J6*J2^7*J10-1048576000*(4478*u^2-2553*u*v-484*v^2)*J4^2*J2^6*J10+1677721600*(8285*u^2+6868*u*v+952*v^2)*J4^3*J2^4*J10+5033164800*(158130*u^2-77586*u*v-47779*v^2)*J6^2*J2^4*J10+26843545600*(243885*u^2-181672*u*v-14308*v^2)*J4^4*J2^2*J10+13589544960000*(155880*u^2-109086*u*v-3529*v^2)*J6^3*J2*J10-36238786560000*(310635*u^2-212922*u*v-6683*v^2)*J6^2*J4^2*J10+65536*(152624*u^2-38848*u*v-12701*v^2)*J6*J4^2*J2^8+3932160*(3768*u^2+55678*u*v-21601*v^2)*J6^2*J4*J2^7-20971520*(28073*u^2-13824*u*v-1834*v^2)*J6*J4^3*J2^6+78643200*(412764*u^2-400948*u*v+63397*v^2)*J6^2*J4^2*J2^5+2831155200*(417480*u^2-277216*u*v-13229*v^2)*J6^3*J4*J2^4-167772160000*(2835840*u^2-1980548*u*v-91647*v^2)*J4*J2^3*J10^2+15099494400000*(103920*u^2-70724*u*v+4439*v^2)*J6*J2^2*J10^2+80530636800000*(105795*u^2-73724*u*v-3061*v^2)*J4^2*J2*J10^2-3623878656000000*(51960*u^2-35862*u*v-1343*v^2)*J6*J4*J10^2-16384*(4*u+v)*(1172*u-237*v)*J6*J4*J2^10-983040*(4*u+v)*(92*u-79*v)*J2^10*J10-131072*(23*u+12*v)*(4*u+v)*J4^3*J2^9-335544320*(3*u+2*v)*(23*u+12*v)*J4^5*J2^5-8493465600*(415680*u^2-286396*u*v-6119*v^2)*J6^4*J2^3-322122547200*(50835*u^2-34112*u*v-1468*v^2)*J6*J4^6-1087163596800*(102795*u^2-70474*u*v-1811*v^2)*J6^3*J4^3+3355443200*(2116880*u^2-1480161*u*v-67329*v^2)*J2^5*J10^2-6012954214400*(21465*u^2-15048*u*v-772*v^2)*J4^5*J10-16384*(7856*u^2+109368*u*v-2049*v^2)*J6^2*J2^9+1048576*(769*u^2+772*u*v+184*v^2)*J4^4*J2^7-19660800*(833376*u^2-559256*u*v-13129*v^2)*J6^3*J2^6)*w2*w3^2+
        32*((4*u+v)^3*J2^11*J15-4*(244*u+51*v)*(4*u+v)^2*J2^9*J4*J15+8*(4*u+v)*(384*u^2-308*u*v-14551*v^2)*J2^8*J6*J15+128*(4*u+v)*(2918*u^2+1229*u*v-172*v^2)*J2^7*J4^2*J15-640*(4032*u^3-768*u^2*v-93680*u*v^2+20041*v^3)*J15*J2^6*J4*J6-128000*(64*u^3+88*u^2*v-36*u*v^2+2443*v^3)*J15*J10*J2^6-5120*(13644*u^3+8695*u^2*v-1892*u*v^2-1272*v^3)*J15*J2^5*J4^3-115200*v*(60*u^2+944*u*v-129*v^2)*J2^5*J6^2*J15+38400*(5184*u^3+948*u^2*v-66812*u*v^2+32925*v^3)*J15*J2^4*J4^2*J6+2560000*(576*u^3+720*u^2*v-276*u*v^2+8191*v^3)*J15*J10*J2^4*J4+1638400*(3*u+2*v)*(324*u^2-19*u*v-100*v^2)*J2^3*J4^4*J15+3456000*v*(240*u^2+2864*u*v-1875*v^2)*J2^3*J4*J6^2*J15+25600000*v^2*(-2129*v+2160*u)*J2^3*J6*J10*J15-18432000*(360*u^3+207*u^2*v-2327*u*v^2+1380*v^3)*J15*J2^2*J4^3*J6-409600000*(216*u^3+243*u^2*v-135*u*v^2+916*v^3)*J15*J10*J2^2*J4^2-2488320000*v^2*(-3*v+4*u)*J2^2*J6^3*J15-393216000*(-3*v+4*u)*(3*u+2*v)^2*J2*J4^5*J15-1658880000*v*(u+9*v)*(15*u-13*v)*J2*J4^2*J6^2*J15-73728000000*v^2*(-61*v+45*u)*J2*J4*J6*J10*J15+36864000000000*v^3*J2*J10^2*J15+8847360000*(-v+u)*(-v+3*u)*(7*v+3*u)*J4^4*J6*J15+65536000000*(3*u+5*v)*(-v+3*u)^2*J4^3*J10*J15+597196800000*v^2*(-v+u)*J4*J6^3*J15-8847360000000*v^3*J6^2*J10*J15)*w3^3;

    return R, C, M;

end function;