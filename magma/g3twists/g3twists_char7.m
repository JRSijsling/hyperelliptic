//freeze;

/***
 *  SL2-invariants of binary octics in characteristic 7
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
 *  Copyright 2013, R. Basson & R. Lercier & C. Ritzenthaler & J. Sijsling
 */

import "../toolbox/sl2invtools.m"    : PowerRepresentativesInFiniteFields, ShiodaInvariantsAppend;
import "g3twists_charp.m" : ShiodaInvariantsCharp;

function ShiodaInvariantsChar7(f :
    PrimaryOnly := false, degmax := Infinity(), degmin := 1)

    JI := []; Wght := [];

    /* Degree 2 */
    if degmax lt 2 then	return JI, Wght; end if;

    if Characteristic(CoefficientRing(f)) eq 0 then

        J2s, J3s, J4s, J5s, J6s, J7s, J8s, J9s, J10s :=
            Explode(ShiodaInvariantsCharp(f : scaled := false));

        /* Removing the content of these invariants */
        j2  := (1/288)                 * J2s;	/* [ <2, 5>,  <3, 2> ]          */
        j3  := (1/103680)              * J3s;	/* [ <2, 8>,  <3, 4>, <5, 1> ]  */
        j4  := (1/2073600)             * J4s;	/* [ <2, 10>, <3, 4>, <5, 2> ]  */
        j5  := (1/298598400)           * J5s;	/* [ <2, 14>, <3, 6>, <5, 2> ] */
        j6  := (1/17915904000)         * J6s;	/* [ <2, 16>, <3, 7>, <5, 3> ] */
        j7  := (1/2579890176000)       * J7s;	/* [ <2, 20>, <3, 9>, <5, 3> ] */
        j8  := (1/154793410560000)     * J8s;	/* [ <2, 22>, <3, 10>, <5, 4> ] */
        j9  := (1/22290251120640000)   * J9s;	/* [ <2, 26>, <3, 12>, <5, 4> ] */
        j10 := (1/1337415067238400000) * J10s;	/* [ <2, 28>, <3, 13>, <5, 5> ] */

        // J2
        if not PrimaryOnly and degmin le 2 then
            Kx := j2;
            Append(~JI, Kx); Append(~Wght, 2);
        end if;
        if degmax le 2 then return JI, Wght; end if;

        // J3
        if degmin le 3 then
            Kx := j3;
            Append(~JI, Kx); Append(~Wght, 3);
        end if;
        if degmax le 3 then return JI, Wght; end if;

        // J4
        if degmin le 4 then
            Kx := -2/7*j2^2 + 1/7*j4;
            Append(~JI, Kx); Append(~Wght, 4);
        end if;
        if degmax le 4 then return JI, Wght; end if;

        // J5
        if degmin le 5 then
            Kx := -5/7*j2*j3 + 1/7*j5;
            Append(~JI, Kx); Append(~Wght, 5);
        end if;
        if degmax le 5 then return JI, Wght; end if;

        // J6
        if degmin le 6 then
            Kx := -18/49*j2^3 - 18/49*j3^2 - 1/49*j2*j4 + 1/49*j6;
            Append(~JI, Kx); Append(~Wght, 6);
        end if;
        if degmax le 6 then return JI, Wght; end if;

        // J7
        if not PrimaryOnly and degmin le 7 then
            Kx := -3/49*j2^2*j3 - 1/49*j3*j4 - 3/49*j2*j5 + 1/49*j7;
            Append(~JI, Kx); Append(~Wght, 7);
        end if;
        if degmax le 7 then return JI, Wght; end if;

        // J8
        if not PrimaryOnly and degmin le 8 then
            Kx := -1/343*j2^4 - 53/343*j2*j3^2 - 11/343*j2^2*j4 - 43/2401*j4^2 - 32/343*j3*j5
                - 6/343*j2*j6 + 1/2401*j8;
            Append(~JI, Kx); Append(~Wght, 8);
        end if;
        if degmax le 8 then return JI, Wght; end if;

        // J9
        if not PrimaryOnly and degmin le 9 then
            Kx := -127/343*j2^3*j3 - 166/343*j3^3 - 15/343*j2*j3*j4 + 9/343*j2^2*j5 -
                3/343*j3*j6 - 6/343*j2*j7 + 1/343*j9;
            Append(~JI, Kx); Append(~Wght, 9);
        end if;
        if degmax le 9 then return JI, Wght; end if;

        // J10
        if degmin le 10 then
            Kx := -654/2401*j2^5 + 915/2401*j2^2*j3^2 + 127/2401*j2^3*j4 - 5/2401*j3^2*j4 +
                258/16807*j2*j4^2 + 45/2401*j2*j3*j5 + 11/2401*j2^2*j6 - 43/16807*j4*j6
                - 15/2401*j3*j7 - 6/16807*j2*j8 + 1/16807*j10;
            Append(~JI, Kx); Append(~Wght, 10);
        end if;
        if degmax le 10 then return JI, Wght; end if;

        // J11
        if not PrimaryOnly and degmin le 11 then
            Kx := 2460/2401*j2^4*j3 + 221/2401*j2*j3^3 + 6/2401*j2^2*j3*j4 + 193/16807*j3*j4^2
                - 169/2401*j2^3*j5 + 169/2401*j3^2*j5 + 1/343*j2*j4*j5 +
                22/2401*j2*j3*j6 - 5/2401*j5*j6 + 8/2401*j2^2*j7 - 3/2401*j4*j7 -
                4/16807*j3*j8 - 5/2401*j2*j9;
            Append(~JI, Kx); Append(~Wght, 11);
        end if;
        if degmax le 11 then return JI, Wght; end if;

        // J13
        if not PrimaryOnly and degmin le 13 then
            Kx := -3305/16807*j2^5*j3 + 3002/16807*j2^2*j3^3 + 1210/16807*j2^3*j3*j4 +
                151/16807*j3^3*j4 - 926/117649*j2*j3*j4^2 - 1483/16807*j2^4*j5 -
                774/16807*j2*j3^2*j5 - 23/2401*j2^2*j4*j5 + 75/16807*j2^2*j3*j6 -
                67/117649*j3*j4*j6 + 2414/16807*j2*j5*j6 + 246/16807*j2^3*j7 +
                19/16807*j3^2*j7 - 2/16807*j2*j4*j7 - 4/16807*j6*j7 + 23/117649*j2*j3*j8
                - 11/16807*j2^2*j9 - 1/16807*j4*j9 - 3/117649*j3*j10;
            Append(~JI, Kx); Append(~Wght, 13);
        end if;
        if degmax le 13 then return JI, Wght; end if;

        // J14
        if degmin le 14 then
            Kx := 6905/823543*j2^7 - 139103/117649*j2^4*j3^2 + 21970/117649*j2*j3^4 -
                890/117649*j2^5*j4 + 7337/117649*j2^2*j3^2*j4 - 12564/823543*j2^3*j4^2 -
                58717/5764801*j3^2*j4^2 - 86/823543*j2*j4^3 - 955/117649*j2^3*j3*j5 -
                60010/823543*j3^3*j5 - 40/16807*j2*j3*j4*j5 - 9/2401*j2^2*j5^2 -
                1535/117649*j2^4*j6 + 3081/117649*j2*j3^2*j6 + 1872/823543*j2^2*j4*j6 +
                2582/5764801*j4^2*j6 + 321/117649*j3*j5*j6 + 6/117649*j2*j6^2 -
                466/117649*j2^2*j3*j7 + 257/117649*j3*j4*j7 + 6/2401*j2*j5*j7 -
                228/823543*j7^2 + 356/823543*j2^3*j8 - 104/5764801*j3^2*j8 -
                12/823543*j2*j4*j8 - 47/5764801*j6*j8 + 27/117649*j2*j3*j9 -
                115/823543*j5*j9 - 31/823543*j2^2*j10 - 29/5764801*j4*j10;
            Append(~JI, Kx); Append(~Wght, 14);
        end if;
        if degmax le 14 then return JI, Wght; end if;

        // J15
        if not PrimaryOnly and degmin le 15 then
            Kx := -94324/117649*j2^6*j3 - 32142/117649*j2^3*j3^3 + 5066/117649*j2^4*j3*j4 +
                1264/117649*j2*j3^3*j4 + 9442/823543*j2^2*j3*j4^2 + 723/823543*j3*j4^3 +
                12645/117649*j2^5*j5 - 6575/117649*j2^2*j3^2*j5 + 62/16807*j2^3*j4*j5 +
                507/117649*j3^2*j4*j5 - 157/117649*j2*j4^2*j5 + 51/16807*j2*j3*j5^2 -
                836/117649*j2^3*j3*j6 + 108/823543*j2*j3*j4*j6 - 2482/117649*j2^2*j5*j6
                - 81/823543*j4*j5*j6 - 173/117649*j2^4*j7 + 698/117649*j2*j3^2*j7 +
                104/117649*j2^2*j4*j7 - 34/823543*j4^2*j7 - 17/16807*j3*j5*j7 -
                4/117649*j2*j6*j7 - 20/823543*j2^2*j3*j8 - 23/823543*j3*j4*j8 +
                16810/117649*j2*j5*j8 - 1/823543*j7*j8 + 218/117649*j2^3*j9 -
                15/117649*j2*j4*j9 + 18/823543*j2*j3*j10 - 3/823543*j5*j10;
            Append(~JI, Kx); Append(~Wght, 15);
        end if;

        return JI, Wght;

    end if;

    a0:= Coefficient(f, 0); a1:= Coefficient(f, 1);
    a2:= Coefficient(f, 2); a3:= Coefficient(f, 3);
    a4:= Coefficient(f, 4); a5:= Coefficient(f, 5);
    a6:= Coefficient(f, 6); a7:= Coefficient(f, 7);
    a8:= Coefficient(f, 8);

    i1:= a4;
    i2:= a0*a8;
    j2:= a1*a7;
    k2:= a2*a6;
    l2:= a3*a5;
    i3:= a0*a5*a7 + a1*a3*a8;
    j3:= a0*a6^2  + a2^2*a8;
    k3:= a1*a5*a6 + a2*a3*a7;
    l3:= a2*a5^2  + a3^2*a6;
    i4:= a0*a6*a5^2  + a8*a2*a3^2;
    j4:= a0*a7*a6*a3 + a8*a1*a2*a5;
    k4:= a0*a7^2*a2  + a8*a1^2*a6;
    l4:= a1*a5^3     + a7*a3^3;
    m4:= a1*a6^2*a3  + a7*a2^2*a5;
    i5:= a0^2*a7^2*a6 + a8^2*a1^2*a2;
    j5:= a0*a5^4      + a8*a3^4;
    k5:= a0*a7^2*a3^2 + a8*a1^2*a5^2;
    l5:= a1^2*a6^3    + a7^2*a2^3;
    i6:= a0^2*a7^3*a3 + a8^2*a1^3*a5;
    i7:= a0^3*a7^4 + a8^3*a1^4;

    /* Degree 2 */
    if not PrimaryOnly and degmin le 2 then
        J2:= 2*i1^2 + 3*k2 + 2*l2;
        Append(~JI, J2); Append(~Wght, 2);
    end if;
    if degmax le 2 then	return JI, Wght; end if;

    /* Degree 3 */
    if degmin le 3 then
        J3:= 5*i1^3 + 2*i1*k2 + 4*i1*l2 + 5*l3;
        Append(~JI, J3); Append(~Wght, 3);
    end if;
    if degmax le 3 then	return JI, Wght; end if;

    /* Degree 4 */
    if degmin le 4 then
        J4:= 4*i1^4 + 3*i1^2*i2 + 4*i1^2*j2 + 6*i1^2*k2 + i2*k2 + 6*j2*k2 + k2^2 + i1^2*l2 + 3*i2*l2 + 4*j2*l2 +
            2*k2*l2 + 6*l2^2 + 6*i1*j3 + 5*i1*k3 + 4*i1*l3 + i4 + 4*l4 + 5*m4;
        Append(~JI, J4); Append(~Wght, 4);
    end if;
    if degmax le 4 then	return JI, Wght; end if;

    /* Degree 5 */
    if degmin le 5 then
        J5:= i1^5 + 3*i1^3*i2 + 4*i1^3*j2 + 3*i1^3*k2 + 4*i1*i2*k2 + 3*i1*j2*k2 + 6*i1^3*l2 + i1*i2*l2 + 6*i1*j2*l2 +
            4*i1*k2*l2 + 3*i1*l2^2 + 3*i1^2*j3 + 3*k2*j3 + 2*l2*j3 + 5*i1^2*k3 + 3*k2*k3 + 5*l2*k3 + 4*i1^2*l3 + 3*i2*l3 +
            4*j2*l3 + 3*k2*l3 + 2*l2*l3 + 5*i1*i4 + 4*i1*l4 + 4*i1*m4 + j5;
        Append(~JI, J5); Append(~Wght, 5);
    end if;
    if degmax le 5 then	return JI, Wght; end if;

    /* Degree 6 */
    if degmin le 6 then
        J6:= 4*i1^4*i2 + 2*i1^2*i2^2 + 3*i1^2*i2*j2 + 2*i1^2*j2^2 + i1^4*k2 + 2*i1^2*i2*k2 + 3*i2^2*k2 + i2*j2*k2 +
            3*j2^2*k2 + 2*i2*k2^2 + 2*j2*k2^2 + 3*i1^2*i2*l2 + 2*i2^2*l2 + 2*i1^2*j2*l2 + 3*i2*j2*l2 + 2*j2^2*l2 +
            4*i1^2*k2*l2 + 6*i2*k2*l2 + 4*j2*k2*l2 + k2^2*l2 + 2*i1^2*l2^2 + 6*i2*l2^2 + j2*l2^2 + 5*k2*l2^2 + 5*l2^3 +
            6*i1^3*i3 + 6*i1*k2*i3 + 2*i1*l2*i3 + 5*i1^3*j3 + i1*i2*j3 + 6*i1*j2*j3 + 4*i1*l2*j3 + 6*j3^2 + 2*i1*i2*k3 +
            5*i1*j2*k3 + i1*k2*k3 + 5*j3*k3 + k3^2 + 4*i1^3*l3 + 6*i1*i2*l3 + 4*i1*j2*l3 + 3*i1*k2*l3 + 6*i1*l2*l3 +
            6*i3*l3 + 4*j3*l3 + 2*l3^2 + 6*i2*i4 + j2*i4 + k2*i4 + 5*i1^2*j4 + 4*k2*j4 + 3*i2*l4 + 4*j2*l4 + 4*i1^2*m4 +
            2*i2*m4 + 5*j2*m4 + 4*k2*m4 + 2*l2*m4 + 5*i1*l5;
        Append(~JI, J6); Append(~Wght, 6);
    end if;
    if degmax le 6 then	return JI, Wght; end if;

    /* Degree 7 */
    if not PrimaryOnly and degmin le 7 then
        J7:= 5*i1^7 + i1^5*i2 + 6*i1^3*i2^2 + 2*i1^5*j2 + 2*i1^3*i2*j2 + 6*i1^3*j2^2 + i1^5*k2 + 6*i1^3*i2*k2 +
            i1*i2^2*k2 + i1^3*j2*k2 + 5*i1*i2*j2*k2 + i1*j2^2*k2 + 2*i1^3*k2^2 + i1*i2*k2^2 + 5*i1*j2*k2^2 + 4*i1^3*i2*l2
            + 2*i1*i2^2*l2 + 3*i1^3*j2*l2 + 3*i1*i2*j2*l2 + 2*i1*j2^2*l2 + 5*i1^3*k2*l2 + 4*i1*i2*k2*l2 + 3*i1*j2*k2*l2 +
            2*i1^3*l2^2 + i1*j2*l2^2 + 5*i1*k2*l2^2 + i1*l2^3 + i1^4*i3 + 4*i1^2*k2*i3 + 5*k2^2*i3 + 2*i1^2*l2*i3 +
            3*k2*l2*i3 + 6*l2^2*i3 + 3*i1^4*j3 + 5*i1^2*i2*j3 + 2*i1^2*j2*j3 + 5*i1^2*k2*j3 + 5*i2*k2*j3 + 2*j2*k2*j3 +
            i2*l2*j3 + 6*j2*l2*j3 + 5*k2*l2*j3 + 5*l2^2*j3 + 3*i1*j3^2 + 6*i1^4*k3 + 6*i1^2*i2*k3 + i1^2*j2*k3 +
            i1^2*k2*k3 + 5*i2*k2*k3 + 2*j2*k2*k3 + 4*k2^2*k3 + 6*i1^2*l2*k3 + 6*i2*l2*k3 + j2*l2*k3 + k2*l2*k3 + l2^2*k3 +
            2*i1*j3*k3 + 3*i1*k3^2 + 4*i1^4*l3 + 5*i1^2*i2*l3 + 6*i2^2*l3 + 5*i1^2*j2*l3 + 2*i2*j2*l3 + 6*j2^2*l3 +
            2*i1^2*k2*l3 + 3*i2*k2*l3 + j2*k2*l3 + 4*k2^2*l3 + 5*i2*l2*l3 + j2*l2*l3 + 3*l2^2*l3 + 3*i1*i3*l3 + 6*i1*j3*l3
            + 4*i1*l3^2 + 2*i1^3*i4 + 6*i1*i2*i4 + i1*j2*i4 + 3*i1*k2*i4 + 6*i1*l2*i4 + 4*j3*i4 + k3*i4 + 2*i1^3*j4 +
            5*i1*k2*j4 + i1*l2*j4 + 2*i1^3*l4 + 2*i1*i2*l4 + 5*i1*j2*l4 + 5*i1*k2*l4 + 4*i1*l2*l4 + 5*j3*l4 + 4*l3*l4 +
            6*i1^3*m4 + 2*i1*i2*m4 + 5*i1*j2*m4 + 5*i1*k2*m4 + 3*i1*l2*m4 + 6*j3*m4 + 2*l3*m4 + 3*i1^2*j5 + 4*i2*j5 +
            3*j2*j5 + 2*k2*j5 + 3*l2*j5 + i1^2*l5 + k2*l5 + 2*l2*l5;
        Append(~JI, J7); Append(~Wght, 7);
    end if;
    if degmax le 7 then	return JI, Wght; end if;

    /* Degree 8 */
    if not PrimaryOnly and degmin le 8 then
        J8:= i1^8 + i1^6*i2 + 3*i1^4*i2^2 + 6*i1^2*i2^3 + 4*i1^6*j2 + 2*i1^4*i2*j2 + 3*i1^2*i2^2*j2 + 2*i1^4*j2^2 +
            4*i1^2*i2*j2^2 + i1^2*j2^3 + i1^4*i2*k2 + 4*i1^2*i2^2*k2 + 2*i2^3*k2 + 2*i1^4*j2*k2 + 3*i1^2*i2*j2*k2 +
            i2^2*j2*k2 + 6*i2*j2^2*k2 + 5*j2^3*k2 + 5*i1^4*k2^2 + 3*i1^2*i2*k2^2 + 4*i2^2*k2^2 + 6*i1^2*j2*k2^2 +
            3*j2^2*k2^2 + 2*j2*k2^3 + 4*i1^6*l2 + 2*i1^4*i2*l2 + 3*i1^2*i2^2*l2 + 6*i2^3*l2 + 4*i1^4*j2*l2 +
            4*i1^2*i2*j2*l2 + 3*i2^2*j2*l2 + 4*i2*j2^2*l2 + j2^3*l2 + 3*i1^4*k2*l2 + 5*i1^2*i2*k2*l2 + i2^2*k2*l2 +
            6*i1^2*j2*k2*l2 + 4*i2*j2*k2*l2 + 2*j2^2*k2*l2 + 2*i1^2*k2^2*l2 + 2*j2*k2^2*l2 + k2^3*l2 + i1^4*l2^2 +
            5*i1^2*i2*l2^2 + 2*i2^2*l2^2 + 3*i1^2*j2*l2^2 + 3*i2*j2*l2^2 + 2*j2^2*l2^2 + 6*i1^2*k2*l2^2 + 6*i2*k2*l2^2 +
            5*j2*k2*l2^2 + 4*k2^2*l2^2 + i1^2*l2^3 + 4*i2*l2^3 + 6*j2*l2^3 + 6*k2*l2^3 + 2*l2^4 + 4*i1^5*i3 + 5*i1^3*i2*i3
            + 2*i1^3*j2*i3 + 5*i1^3*k2*i3 + 5*i1*i2*k2*i3 + 2*i1*j2*k2*i3 + 5*i1*k2^2*i3 + 3*i1^3*l2*i3 + 4*i1*i2*l2*i3 +
            3*i1*j2*l2*i3 + 2*i1*k2*l2*i3 + 4*i1*l2^2*i3 + i1^5*j3 + 2*i1^3*i2*j3 + i1*i2^2*j3 + i1^3*j2*j3 +
            5*i1*i2*j2*j3 + i1*j2^2*j3 + i1^3*k2*j3 + 6*i1*i2*k2*j3 + i1*j2*k2*j3 + 2*i1*k2^2*j3 + 4*i1*i2*l2*j3 +
            4*i1*j2*l2*j3 + 6*i1*k2*l2*j3 + i1*l2^2*j3 + i1^2*i3*j3 + 4*k2*i3*j3 + 4*l2*i3*j3 + 2*i1^2*j3^2 + 5*i2*j3^2 +
            2*j2*j3^2 + 5*i1^5*k3 + 2*i1*i2^2*k3 + 2*i1^3*j2*k3 + 3*i1*i2*j2*k3 + 2*i1*j2^2*k3 + i1^3*k2*k3 +
            2*i1*i2*k2*k3 + i1*j2*k2*k3 + 2*i1*k2^2*k3 + i1^3*l2*k3 + 2*i1*i2*l2*k3 + 5*i1*j2*l2*k3 + 2*i1*k2*l2*k3 +
            2*i1*l2^2*k3 + 6*i1^2*j3*k3 + 3*i2*j3*k3 + 4*j2*j3*k3 + 6*k2*j3*k3 + 5*l2*j3*k3 + i1^2*k3^2 + 2*i2*k3^2 +
            5*j2*k3^2 + 3*k2*k3^2 + 4*i1^5*l3 + 5*i1^3*i2*l3 + 2*i1*i2^2*l3 + 4*i1^3*j2*l3 + 2*i1*i2*j2*l3 + 3*i1*j2^2*l3
            + 6*i1^3*k2*l3 + 2*i1*i2*k2*l3 + 4*i1*k2^2*l3 + i1^3*l2*l3 + i1*i2*l2*l3 + 2*i1*j2*l2*l3 + 6*i1*k2*l2*l3 +
            5*i2*i3*l3 + 2*j2*i3*l3 + 6*k2*i3*l3 + 4*l2*i3*l3 + 4*i1^2*j3*l3 + 5*i2*j3*l3 + 6*j2*j3*l3 + 2*i1^2*l3^2 +
            2*i2*l3^2 + 4*j2*l3^2 + k2*l3^2 + 3*l2*l3^2 + 2*i1^4*i4 + 6*i1^2*i2*i4 + 6*i2^2*i4 + 3*i1^2*j2*i4 + 2*i2*j2*i4
            + 6*j2^2*i4 + 4*i1^2*k2*i4 + i2*k2*i4 + 5*j2*k2*i4 + 4*k2^2*i4 + 6*i1^2*l2*i4 + 5*j2*l2*i4 + 5*k2*l2*i4 +
            3*l2^2*i4 + 4*i1*i3*i4 + 3*i1*j3*i4 + 4*i1*k3*i4 + i1^4*j4 + 3*i1^2*i2*j4 + 4*i1^2*j2*j4 + i2*k2*j4 +
            6*j2*k2*j4 + k2^2*j4 + 6*i1^2*l2*j4 + k2*l2*j4 + l2^2*j4 + 2*i1*j3*j4 + 5*i1^4*k4 + 5*i1^2*k2*k4 + 5*k2^2*k4 +
            2*i1^2*l2*k4 + 5*k2*l2*k4 + 3*l2^2*k4 + 4*i1^4*l4 + 6*i1^2*i2*l4 + 3*i2^2*l4 + 6*i1^2*j2*l4 + i2*j2*l4 +
            3*j2^2*l4 + 3*i1^2*k2*l4 + 3*i2*k2*l4 + 6*j2*k2*l4 + 6*k2^2*l4 + 4*i1^2*l2*l4 + 5*j2*l2*l4 + 2*k2*l2*l4 +
            5*l2^2*l4 + 6*i1*j3*l4 + 5*i1*k3*l4 + 5*i1*l3*l4 + 2*l4^2 + i1^4*m4 + 3*i1^2*i2*m4 + 2*i2^2*m4 + 2*i1^2*j2*m4
            + 3*i2*j2*m4 + 2*j2^2*m4 + 6*i1^2*k2*m4 + 3*i2*k2*m4 + 3*j2*k2*m4 + 4*k2^2*m4 + 5*i2*l2*m4 + 2*j2*l2*m4 +
            2*k2*l2*m4 + 6*l2^2*m4 + 2*i1*i3*m4 + 4*i1*j3*m4 + 4*i1*l3*m4 + 2*i1^3*j5 + 2*i1*i2*j5 + 4*i1*j2*j5 +
            2*i1*k2*j5 + 3*i1*l2*j5 + 4*i3*j5 + 2*j3*j5 + k3*j5 + 2*l3*j5 + 6*i1^3*k5 + 6*i1*k2*k5 + 6*i1*l2*k5 + l3*k5 +
            i1^3*l5 + 3*i1*i2*l5 + 4*i1*j2*l5 + 3*i1*k2*l5 + 4*i1*l2*l5 + 5*j3*l5 + 5*k3*l5 + 5*l3*l5;

        Append(~JI, J8); Append(~Wght, 8);
    end if;
    if degmax le 8 then	return JI, Wght; end if;

    /* Degree 9 */
    if not PrimaryOnly and degmin le 9 then
        J9:= 5*i1^5*i2^2 + 2*i1^3*i2^3 + i1^7*j2 + i1^3*i2^2*j2 + 2*i1^5*j2^2 + 6*i1^3*i2*j2^2 + 5*i1^3*j2^3 +
            5*i1^7*k2 + 6*i1^5*i2*k2 + 3*i1^3*i2^2*k2 + 5*i1*i2^3*k2 + 2*i1^5*j2*k2 + i1^3*i2*j2*k2 + 6*i1*i2^2*j2*k2 +
            3*i1^3*j2^2*k2 + i1*i2*j2^2*k2 + 2*i1*j2^3*k2 + 5*i1^3*i2*k2^2 + 6*i1*i2^2*k2^2 + i1*i2*j2*k2^2 + 6*i1^3*k2^3
            + 6*i1*i2*k2^3 + 6*i1*j2*k2^3 + 2*i1*k2^4 + 4*i1^5*i2*l2 + 3*i1*i2^3*l2 + 4*i1^5*j2*l2 + 5*i1*i2^2*j2*l2 +
            2*i1*i2*j2^2*l2 + 4*i1*j2^3*l2 + 5*i1^5*k2*l2 + 4*i1^3*i2*k2*l2 + 2*i1^3*j2*k2*l2 + 3*i1*i2*k2^2*l2 +
            4*i1*j2*k2^2*l2 + 3*i1^5*l2^2 + 6*i1^3*i2*l2^2 + 3*i1*i2^2*l2^2 + 5*i1^3*j2*l2^2 + 2*i1*i2*j2*l2^2 +
            2*i1*j2^2*l2^2 + 6*i1^3*k2*l2^2 + 6*i1*i2*k2*l2^2 + 6*i1*j2*k2*l2^2 + 5*i1*k2^2*l2^2 + 6*i1^3*l2^3 +
            5*i1*i2*l2^3 + 3*i1*j2*l2^3 + 5*i1*k2*l2^3 + i1*l2^4 + 5*i1^6*i3 + i1^4*i2*i3 + 6*i1^4*j2*i3 + 4*i1^2*i2*k2*i3
            + 3*i1^2*j2*k2*i3 + 6*i1^2*k2^2*i3 + 5*i2*k2^2*i3 + 2*j2*k2^2*i3 + 2*k2^3*i3 + i1^4*l2*i3 + 2*i1^2*i2*l2*i3 +
            5*i1^2*j2*l2*i3 + 3*i2*k2*l2*i3 + 4*j2*k2*l2*i3 + 6*k2^2*l2*i3 + 4*i1^2*l2^2*i3 + 6*i2*l2^2*i3 + j2*l2^2*i3 +
            2*k2*l2^2*i3 + 4*l2^3*i3 + 4*i1^6*j3 + i1^4*i2*j3 + 6*i1^2*i2^2*j3 + i1^4*j2*j3 + 2*i1^2*i2*j2*j3 +
            6*i1^2*j2^2*j3 + 5*i1^4*k2*j3 + 3*i1^2*i2*k2*j3 + 6*i2^2*k2*j3 + 3*i1^2*j2*k2*j3 + 2*i2*j2*k2*j3 +
            6*j2^2*k2*j3 + 6*i1^2*k2^2*j3 + 2*i2*k2^2*j3 + 4*j2*k2^2*j3 + 2*k2^3*j3 + 4*i1^4*l2*j3 + 4*i2^2*l2*j3 +
            6*i1^2*j2*l2*j3 + 6*i2*j2*l2*j3 + 4*j2^2*l2*j3 + 2*i2*k2*l2*j3 + 5*j2*k2*l2*j3 + 2*k2^2*l2*j3 + 3*i1^2*l2^2*j3
            + i2*l2^2*j3 + 3*j2*l2^2*j3 + 6*k2*l2^2*j3 + 2*l2^3*j3 + 3*i1^3*i3*j3 + 6*i1*l2*i3*j3 + 3*i1*i2*j3^2 +
            4*i1*j2*j3^2 + 2*i1*k2*j3^2 + 5*i1*l2*j3^2 + 2*j3^3 + 3*i1^6*k3 + 5*i1^4*i2*k3 + 3*i1^2*i2^2*k3 + 4*i1^4*j2*k3
            + i1^2*i2*j2*k3 + 3*i1^2*j2^2*k3 + 2*i1^4*k2*k3 + i1^2*i2*k2*k3 + 6*i2^2*k2*k3 + 2*i2*j2*k2*k3 + 6*j2^2*k2*k3
            + 3*i1^2*k2^2*k3 + 6*i2*k2^2*k3 + 5*k2^3*k3 + 5*i1^4*l2*k3 + 3*i2^2*l2*k3 + 3*i1^2*j2*l2*k3 + i2*j2*l2*k3 +
            3*j2^2*l2*k3 + 6*i1^2*k2*l2*k3 + 6*i2*k2*l2*k3 + 6*j2*k2*l2*k3 + k2^2*l2*k3 + 4*i1^2*l2^2*k3 + 2*i2*l2^2*k3 +
            5*j2*l2^2*k3 + 6*k2*l2^2*k3 + 2*l2^3*k3 + 5*i1^3*j3*k3 + 2*i1*i2*j3*k3 + 5*i1*j2*j3*k3 + 6*i1*k2*j3*k3 +
            2*i1*l2*j3*k3 + 6*j3^2*k3 + 3*i1*i2*k3^2 + 4*i1*j2*k3^2 + 2*i1*k2*k3^2 + 6*i1*l2*k3^2 + j3*k3^2 + 6*i1^6*l3 +
            4*i1^4*i2*l3 + 3*i1^2*i2^2*l3 + 2*i2^3*l3 + 2*i1^4*j2*l3 + 4*i1^2*i2*j2*l3 + i2^2*j2*l3 + 6*i2*j2^2*l3 +
            5*j2^3*l3 + 2*i1^4*k2*l3 + 2*i1^2*i2*k2*l3 + 6*i2^2*k2*l3 + 3*i1^2*j2*k2*l3 + 6*i2*j2*k2*l3 + 2*j2^2*k2*l3 +
            6*i1^2*k2^2*l3 + i2*k2^2*l3 + j2*k2^2*l3 + 3*k2^3*l3 + 6*i1^2*j2*l2*l3 + 6*i2*j2*l2*l3 + j2^2*l2*l3 +
            5*j2*k2*l2*l3 + 2*k2^2*l2*l3 + 2*i1^2*l2^2*l3 + 2*i2*l2^2*l3 + 3*j2*l2^2*l3 + 4*k2*l2^2*l3 + 3*l2^3*l3 +
            6*i1^3*i3*l3 + 3*i1*i2*i3*l3 + 4*i1*j2*i3*l3 + 6*i1*k2*i3*l3 + i1*l2*i3*l3 + 3*i1^3*j3*l3 + 6*i1*i2*j3*l3 +
            2*i1*l2*j3*l3 + 6*j3^2*l3 + 4*i1^3*l3^2 + 6*i1*j2*l3^2 + 5*i1*k2*l3^2 + 2*i3*l3^2 + j3*l3^2 + 2*l3^3 +
            5*i1^5*i4 + i1^3*i2*i4 + 3*i1*i2^2*i4 + 6*i1^3*j2*i4 + i1*i2*j2*i4 + 3*i1*j2^2*i4 + 3*i1^3*k2*i4 + i1*i2*k2*i4
            + 3*i1*j2*k2*i4 + 5*i1*k2^2*i4 + i1^3*l2*i4 + 3*i1*j2*l2*i4 + 5*i1*k2*l2*i4 + 3*i1*l2^2*i4 + 2*i1^2*i3*i4 +
            2*k2*i3*i4 + 2*l2*i3*i4 + 3*i1^2*j3*i4 + 4*i2*j3*i4 + 3*j2*j3*i4 + 5*k2*j3*i4 + 3*l2*j3*i4 + 5*i1^2*k3*i4 +
            i2*k3*i4 + 6*j2*k3*i4 + 4*k2*k3*i4 + 4*l2*k3*i4 + 3*i1^5*j4 + 2*i1^3*i2*j4 + 5*i1^3*j2*j4 + 5*i1^3*k2*j4 +
            5*i1*i2*k2*j4 + 2*i1*j2*k2*j4 + 5*i1^3*l2*j4 + i1*i2*l2*j4 + 6*i1*j2*l2*j4 + i1*k2*l2*j4 + 6*i1*l2^2*j4 +
            6*i1^2*j3*j4 + 6*k2*j3*j4 + l2*j3*j4 + 4*i1^5*k4 + i1^3*k2*k4 + 2*i1*k2^2*k4 + 4*i1*k2*l2*k4 + 6*i1*l2^2*k4 +
            i1^5*l4 + 6*i1^3*i2*l4 + i1*i2^2*l4 + 5*i1*i2*j2*l4 + i1*j2^2*l4 + 2*i1^3*k2*l4 + i1*i2*k2*l4 + i1*j2*k2*l4 +
            4*i1*k2^2*l4 + i1^3*l2*l4 + 2*i1*i2*l2*l4 + 3*i1*j2*l2*l4 + 4*i1*k2*l2*l4 + 6*i1^2*j3*l4 + 5*i2*j3*l4 +
            2*j2*j3*l4 + 3*k2*j3*l4 + 5*k2*k3*l4 + 5*i1^2*l3*l4 + 6*i2*l3*l4 + 6*j2*l3*l4 + 4*k2*l3*l4 + 5*l2*l3*l4 +
            5*i1^5*m4 + i1^3*i2*m4 + i1*i2^2*m4 + 6*i1^3*j2*m4 + 5*i1*i2*j2*m4 + i1*j2^2*m4 + 5*i1^3*k2*m4 + 2*i1*i2*k2*m4
            + 2*i1*j2*k2*m4 + 4*i1*k2^2*m4 + i1^3*l2*m4 + i1*i2*l2*m4 + i1*k2*l2*m4 + 5*i1*l2^2*m4 + i1^2*i3*m4 +
            3*k2*i3*m4 + 6*l2*i3*m4 + 3*i1^2*j3*m4 + 6*i2*j3*m4 + j2*j3*m4 + 6*k2*j3*m4 + 2*l2*j3*m4 + 5*i1^2*l3*m4 +
            2*i2*l3*m4 + 2*j2*l3*m4 + 5*k2*l3*m4 + 2*l2*l3*m4 + 6*i1^4*j5 + 6*i1^2*i2*j5 + 2*i2^2*j5 + 2*i1^2*j2*j5 +
            3*i2*j2*j5 + 2*j2^2*j5 + 4*j2*k2*j5 + 5*i1^2*l2*j5 + 6*i2*l2*j5 + 2*j2*l2*j5 + 6*l2^2*j5 + 2*i1^4*k5 +
            5*i1^2*k2*k5 + 3*k2^2*k5 + 6*i1^2*l2*k5 + l2^2*k5 + 6*i1*l3*k5 + 3*i1^4*l5 + i1^2*i2*l5 + 6*i1^2*j2*l5 +
            4*i1^2*k2*l5 + i2*k2*l5 + 6*j2*k2*l5 + 6*k2^2*l5 + 3*i1^2*l2*l5 + 2*i2*l2*l5 + 5*j2*l2*l5 + 2*k2*l2*l5 +
            2*l2^2*l5 + 5*i1*j3*l5 + 4*i1*k3*l5 + 2*i1*l3*l5 + 5*m4*l5;
        Append(~JI, J9); Append(~Wght, 9);
    end if;
    if degmax le 9 then	return JI, Wght; end if;

    /* Degree 10 */
    if degmin le 10 then
        J10:= 5*i1^10 + 2*i1^8*i2 + 3*i1^6*i2^2 + 2*i1^4*i2^3 + 5*i1^2*i2^4 + 5*i1^4*i2^2*j2 + i1^2*i2^3*j2 +
            2*i1^6*j2^2 + 5*i1^4*i2*j2^2 + 2*i1^2*i2^2*j2^2 + 2*i1^4*j2^3 + i1^2*i2*j2^3 + 5*i1^2*j2^4 + 2*i1^8*k2 +
            i1^6*i2*k2 + 2*i1^4*i2^2*k2 + 6*i1^2*i2^3*k2 + 4*i2^4*k2 + 4*i1^6*j2*k2 + 4*i1^4*i2*j2*k2 + 5*i1^2*i2^2*j2*k2
            + 5*i2^3*j2*k2 + 3*i2^2*j2^2*k2 + 3*i1^2*j2^3*k2 + 5*i2*j2^3*k2 + 4*j2^4*k2 + 2*i1^6*k2^2 + 5*i2^3*k2^2 +
            5*i1^4*j2*k2^2 + 5*i1^2*i2*j2*k2^2 + 3*i2^2*j2*k2^2 + 5*i1^2*j2^2*k2^2 + 6*j2^3*k2^2 + 3*i1^4*k2^3 +
            6*i2^2*k2^3 + 6*i1^2*j2*k2^3 + 3*i2*j2*k2^3 + 5*i1^2*k2^4 + 5*i2*k2^4 + 2*j2*k2^4 + 2*k2^5 + 4*i1^8*l2 +
            2*i1^6*i2*l2 + 5*i1^4*i2^2*l2 + 6*i1^2*i2^3*l2 + 5*i2^4*l2 + i1^6*j2*l2 + i1^2*i2^2*j2*l2 + i2^3*j2*l2 +
            i1^4*j2^2*l2 + i1^2*i2*j2^2*l2 + 2*i2^2*j2^2*l2 + 6*i1^2*j2^3*l2 + i2*j2^3*l2 + 5*j2^4*l2 + 6*i1^4*i2*k2*l2 +
            i1^2*i2^2*k2*l2 + 5*i2^3*k2*l2 + 5*i1^4*j2*k2*l2 + 3*i1^2*i2*j2*k2*l2 + 2*i2^2*j2*k2*l2 + 6*i1^2*j2^2*k2*l2 +
            2*i2*j2^2*k2*l2 + 5*j2^3*k2*l2 + 3*i1^4*k2^2*l2 + 5*i1^2*i2*k2^2*l2 + 2*i2^2*k2^2*l2 + 6*i1^2*j2*k2^2*l2 +
            5*i2*j2*k2^2*l2 + 6*i2*k2^3*l2 + j2*k2^3*l2 + 3*k2^4*l2 + 3*i1^4*i2*l2^2 + 3*i1^2*i2^2*l2^2 + 3*i2^3*l2^2 +
            4*i1^2*i2*j2*l2^2 + 5*i2^2*j2*l2^2 + 2*i2*j2^2*l2^2 + 4*j2^3*l2^2 + i1^4*k2*l2^2 + 3*i1^2*i2*k2*l2^2 +
            2*i1^2*j2*k2*l2^2 + 4*i2*j2*k2*l2^2 + 6*j2^2*k2*l2^2 + 2*i1^2*k2^2*l2^2 + 3*j2*k2^2*l2^2 + 3*k2^3*l2^2 +
            4*i1^2*i2*l2^3 + 3*i1^2*j2*l2^3 + 3*i2*j2*l2^3 + 4*j2^2*l2^3 + 3*i1^2*k2*l2^3 + i2*k2*l2^3 + j2*k2*l2^3 +
            2*k2^2*l2^3 + 5*i1^2*l2^4 + 5*i2*l2^4 + 2*j2*l2^4 + 3*i1^7*i3 + 2*i1^5*i2*i3 + 6*i1^3*i2^2*i3 + 6*i1^5*j2*i3 +
            2*i1^3*i2*j2*i3 + 6*i1^3*j2^2*i3 + 2*i1^5*k2*i3 + 3*i1^3*i2*k2*i3 + 6*i1*i2^2*k2*i3 + 2*i1^3*j2*k2*i3 +
            2*i1*i2*j2*k2*i3 + 6*i1*j2^2*k2*i3 + i1*i2*k2^2*i3 + 2*i1*j2*k2^2*i3 + i1*k2^3*i3 + 6*i1^3*i2*l2*i3 +
            2*i1*i2^2*l2*i3 + 2*i1^3*j2*l2*i3 + 3*i1*i2*j2*l2*i3 + 2*i1*j2^2*l2*i3 + 4*i1^3*k2*l2*i3 + 6*i1*i2*k2*l2*i3 +
            6*i1*j2*k2*l2*i3 + 5*i1*k2^2*l2*i3 + 6*i1^3*l2^2*i3 + 6*i1*i2*l2^2*i3 + 2*i1*j2*l2^2*i3 + 4*i1*k2*l2^2*i3 +
            i1*l2^3*i3 + 5*i1^4*i3^2 + 3*i1^2*k2*i3^2 + 2*k2^2*i3^2 + 3*i1^2*l2*i3^2 + k2*l2*i3^2 + 3*l2^2*i3^2 +
            4*i1^7*j3 + 6*i1^5*i2*j3 + 5*i1^3*i2^2*j3 + 5*i1*i2^3*j3 + 3*i1^5*j2*j3 + 6*i1*i2^2*j2*j3 + 2*i1^3*j2^2*j3 +
            i1*i2*j2^2*j3 + 2*i1*j2^3*j3 + i1^5*k2*j3 + 4*i1^3*i2*k2*j3 + 4*i1*i2^2*k2*j3 + 5*i1^3*j2*k2*j3 +
            6*i1*i2*j2*k2*j3 + 4*i1*j2^2*k2*j3 + 3*i1^3*k2^2*j3 + 3*i1*i2*k2^2*j3 + 2*i1*j2*k2^2*j3 + 6*i1^5*l2*j3 +
            6*i1^3*i2*l2*j3 + i1*i2^2*l2*j3 + 5*i1^3*j2*l2*j3 + 6*i1*i2*j2*l2*j3 + 4*i1^3*k2*l2*j3 + 5*i1*i2*k2*l2*j3 +
            6*i1*k2^2*l2*j3 + 4*i1^3*l2^2*j3 + 5*i1*j2*l2^2*j3 + 4*i1*k2*l2^2*j3 + i1*l2^3*j3 + 5*i1^4*i3*j3 +
            i1^2*i2*i3*j3 + 6*i1^2*j2*i3*j3 + 3*i1^2*k2*i3*j3 + 4*i2*k2*i3*j3 + 3*j2*k2*i3*j3 + 3*k2^2*i3*j3 +
            3*i1^2*l2*i3*j3 + 4*i2*l2*i3*j3 + 3*j2*l2*i3*j3 + k2*l2*i3*j3 + 2*l2^2*i3*j3 + 4*i1^4*j3^2 + 6*i1^2*i2*j3^2 +
            6*i2^2*j3^2 + 2*i1^2*j2*j3^2 + 2*i2*j2*j3^2 + 6*j2^2*j3^2 + 6*i1^2*k2*j3^2 + 6*i2*k2*j3^2 + 4*j2*k2*j3^2 +
            5*k2^2*j3^2 + 3*i1^2*l2*j3^2 + 5*k2*l2*j3^2 + 5*l2^2*j3^2 + 5*i1*i3*j3^2 + 2*i1*j3^3 + 2*i1^3*i2^2*k3 +
            3*i1*i2^3*k3 + 5*i1^5*j2*k3 + 5*i1^3*i2*j2*k3 + 5*i1*i2^2*j2*k3 + 2*i1*i2*j2^2*k3 + 4*i1*j2^3*k3 +
            3*i1^5*k2*k3 + 5*i1^3*i2*k2*k3 + 4*i1*i2^2*k2*k3 + i1^3*j2*k2*k3 + 2*i1*i2*j2*k2*k3 + i1*j2^2*k2*k3 +
            5*i1^3*k2^2*k3 + i1*i2*k2^2*k3 + i1*j2*k2^2*k3 + 2*i1*k2^3*k3 + 4*i1^5*l2*k3 + 5*i1^3*i2*l2*k3 +
            4*i1*i2^2*l2*k3 + 4*i1^3*j2*l2*k3 + 6*i1*i2*j2*l2*k3 + 4*i1*j2^2*l2*k3 + 6*i1^3*k2*l2*k3 + 2*i1*i2*k2*l2*k3 +
            5*i1*j2*k2*l2*k3 + 4*i1*k2^2*l2*k3 + 3*i1^3*l2^2*k3 + 3*i1*i2*l2^2*k3 + 5*i1*j2*l2^2*k3 + 6*i1*k2*l2^2*k3 +
            5*i1*l2^3*k3 + i1^4*j3*k3 + 2*i1^2*i2*j3*k3 + 5*i2^2*j3*k3 + 5*i1^2*j2*j3*k3 + 4*i2*j2*j3*k3 + 5*j2^2*j3*k3 +
            2*i1^2*k2*j3*k3 + 4*i2*k2*j3*k3 + 5*j2*k2*j3*k3 + 4*k2^2*j3*k3 + 2*i1^2*l2*j3*k3 + 2*i2*l2*j3*k3 +
            5*j2*l2*j3*k3 + 3*k2*l2*j3*k3 + 5*l2^2*j3*k3 + 5*i1*j3^2*k3 + 4*i1^2*i2*k3^2 + i2^2*k3^2 + 5*i2*j2*k3^2 +
            j2^2*k3^2 + 2*i1^2*k2*k3^2 + i2*k2*k3^2 + 4*j2*k2*k3^2 + 2*i2*l2*k3^2 + 2*j2*l2*k3^2 + k2*l2*k3^2 + l2^2*k3^2
            + 5*i1*j3*k3^2 + 2*i1*k3^3 + 5*i1^7*l3 + i1^3*i2^2*l3 + 4*i1*i2^3*l3 + 4*i1^5*j2*l3 + 4*i1^3*i2*j2*l3 +
            5*i1*i2^2*j2*l3 + 2*i1^3*j2^2*l3 + 6*i1*i2*j2^2*l3 + 6*i1*j2^3*l3 + 4*i1^5*k2*l3 + i1^3*i2*k2*l3 +
            5*i1*i2^2*k2*l3 + 2*i1^3*j2*k2*l3 + 4*i1*i2*j2*k2*l3 + 3*i1*j2^2*k2*l3 + 3*i1*i2*k2^2*l3 + 4*i1*j2*k2^2*l3 +
            5*i1*k2^3*l3 + 4*i1^3*i2*l2*l3 + 2*i1*i2^2*l2*l3 + 6*i1*i2*j2*l2*l3 + 2*i1*j2^2*l2*l3 + 2*i1^3*k2*l2*l3 +
            i1*i2*k2*l2*l3 + 4*i1*j2*k2*l2*l3 + 6*i1*k2^2*l2*l3 + 6*i1^3*l2^2*l3 + 4*i1*i2*l2^2*l3 + 2*i1*k2*l2^2*l3 +
            2*i1*l2^3*l3 + 5*i1^4*i3*l3 + i1^2*i2*i3*l3 + 6*i2^2*i3*l3 + i1^2*j2*i3*l3 + 2*i2*j2*i3*l3 + 6*j2^2*i3*l3 +
            5*i1^2*k2*i3*l3 + 4*i2*k2*i3*l3 + 6*i1^2*l2*i3*l3 + 2*i2*l2*i3*l3 + 2*j2*l2*i3*l3 + 6*k2*l2*i3*l3 + l2^2*i3*l3
            + 3*i1*i3^2*l3 + 3*i1^4*j3*l3 + 5*i1^2*i2*j3*l3 + i2^2*j3*l3 + 2*i2*j2*j3*l3 + 4*j2^2*j3*l3 + 6*i2*k2*j3*l3 +
            5*j2*k2*j3*l3 + k2^2*j3*l3 + 4*i1^2*l2*j3*l3 + 5*i2*l2*j3*l3 + 6*j2*l2*j3*l3 + 2*k2*l2*j3*l3 + l2^2*j3*l3 +
            6*i1*j3^2*l3 + i1^4*l3^2 + 6*i1^2*i2*l3^2 + 4*i2^2*l3^2 + 3*i2*j2*l3^2 + 5*j2^2*l3^2 + i1^2*k2*l3^2 +
            6*j2*k2*l3^2 + 2*k2^2*l3^2 + i2*l2*l3^2 + 4*j2*l2*l3^2 + 6*k2*l2*l3^2 + 3*l2^2*l3^2 + 3*i1*i3*l3^2 +
            2*i1*j3*l3^2 + 6*i1*l3^3 + 5*i1^6*i4 + 2*i1^4*i2*i4 + 5*i1^2*i2^2*i4 + 2*i2^3*i4 + 4*i1^4*j2*i4 +
            6*i1^2*i2*j2*i4 + i2^2*j2*i4 + 3*i1^2*j2^2*i4 + 6*i2*j2^2*i4 + 5*j2^3*i4 + i1^4*k2*i4 + 2*i2^2*k2*i4 +
            4*i1^2*j2*k2*i4 + 2*i2*j2*k2*i4 + 3*j2^2*k2*i4 + 5*i1^2*k2^2*i4 + 6*j2*k2^2*i4 + 3*k2^3*i4 + 4*i1^4*l2*i4 +
            i1^2*i2*l2*i4 + 6*i2^2*l2*i4 + 5*i1^2*j2*l2*i4 + j2^2*l2*i4 + 2*i1^2*k2*l2*i4 + 5*i2*k2*l2*i4 + 3*j2*k2*l2*i4
            + 2*k2^2*l2*i4 + 2*i1^2*l2^2*i4 + i2*l2^2*i4 + 6*j2*l2^2*i4 + k2*l2^2*i4 + 4*i1*i2*i3*i4 + 3*i1*j2*i3*i4 +
            5*i1*k2*i3*i4 + 3*i1*l2*i3*i4 + 4*i1^3*j3*i4 + i1*i2*j3*i4 + 4*i1*k2*j3*i4 + 3*i1*l2*j3*i4 + 6*i3*j3*i4 +
            5*j3^2*i4 + i1^3*k3*i4 + 2*i1*j2*k3*i4 + 3*i1*k2*k3*i4 + 5*i1*l2*k3*i4 + 4*j3*k3*i4 + 6*i1^6*j4 + 5*i1^4*i2*j4
            + 5*i1^2*i2^2*j4 + 6*i1^4*j2*j4 + 4*i1^2*i2*j2*j4 + 5*i1^2*j2^2*j4 + 5*i1^4*k2*j4 + i1^2*i2*k2*j4 +
            4*i2^2*k2*j4 + 6*i1^2*j2*k2*j4 + 6*i2*j2*k2*j4 + 4*j2^2*k2*j4 + 6*i1^2*k2^2*j4 + 5*j2*k2^2*j4 + 2*k2^3*j4 +
            5*i1^4*l2*j4 + 5*i1^2*i2*l2*j4 + 6*i1^2*j2*l2*j4 + 2*i1^2*k2*l2*j4 + 3*j2*k2*l2*j4 + 5*i1^2*l2^2*j4 +
            6*i2*l2^2*j4 + 5*j2*l2^2*j4 + 5*k2*l2^2*j4 + 4*l2^3*j4 + 3*i1^3*j3*j4 + 2*i1*i2*j3*j4 + 5*i1*j2*j3*j4 +
            i1*l2*j3*j4 + 3*j3^2*j4 + 5*i1^4*i2*k4 + 2*i1^4*j2*k4 + 5*i1^2*i2*k2*k4 + 2*i1^2*j2*k2*k4 + 6*i1^2*k2^2*k4 +
            5*i2*k2^2*k4 + 2*j2*k2^2*k4 + 4*k2^3*k4 + 4*i1^4*l2*k4 + 2*i1^2*i2*l2*k4 + 5*i1^2*j2*l2*k4 + 4*i1^2*k2*l2*k4 +
            5*i2*k2*l2*k4 + 2*j2*k2*l2*k4 + 6*k2^2*l2*k4 + i1^2*l2^2*k4 + 3*i2*l2^2*k4 + 4*j2*l2^2*k4 + 5*k2*l2^2*k4 +
            4*l2^3*k4 + 6*i1^3*k3*k4 + 2*i1*k2*k3*k4 + 2*i1*l2*k3*k4 + 2*i1^2*i2^2*l4 + i2^3*l4 + 4*i1^4*j2*l4 +
            i1^2*i2*j2*l4 + 4*i2^2*j2*l4 + 4*i1^2*j2^2*l4 + 3*i2*j2^2*l4 + 6*j2^3*l4 + 4*i1^4*k2*l4 + 4*i2^2*k2*l4 +
            i1^2*j2*k2*l4 + i2*j2*k2*l4 + 2*j2^2*k2*l4 + 4*i1^2*k2^2*l4 + 3*i2*k2^2*l4 + 5*j2*k2^2*l4 + 2*k2^3*l4 +
            i1^4*l2*l4 + 3*i1^2*i2*l2*l4 + 3*i2^2*l2*l4 + 6*i2*j2*l2*l4 + 5*j2^2*l2*l4 + 2*i1^2*k2*l2*l4 + i2*k2*l2*l4 +
            2*j2*k2*l2*l4 + 5*k2^2*l2*l4 + 6*i1^2*l2^2*l4 + 4*i2*l2^2*l4 + 4*k2*l2^2*l4 + 4*i1^3*j3*l4 + 2*i1*i2*j3*l4 +
            4*i1*j2*j3*l4 + 2*i1*k2*j3*l4 + 3*i1*l2*j3*l4 + 3*i1^3*k3*l4 + 6*i1*i2*k3*l4 + 2*i1*k2*k3*l4 + 3*i1*l2*k3*l4 +
            3*j3*k3*l4 + 4*k3^2*l4 + 4*i1^3*l3*l4 + 5*i1*i2*l3*l4 + 4*i1*j2*l3*l4 + i1*k2*l3*l4 + 3*i1*l2*l3*l4 +
            6*l3^2*l4 + 4*i1^2*l4^2 + 2*j2*l4^2 + 6*k2*l4^2 + 4*l2*l4^2 + 6*i1^6*m4 + i1^2*i2^2*m4 + 3*i2^3*m4 +
            4*i1^4*j2*m4 + 3*i1^2*i2*j2*m4 + 5*i2^2*j2*m4 + 3*i1^2*j2^2*m4 + 2*i2*j2^2*m4 + 4*j2^3*m4 + 2*i1^4*k2*m4 +
            5*i1^2*i2*k2*m4 + 2*i2^2*k2*m4 + 2*i1^2*j2*k2*m4 + 2*i2*j2*k2*m4 + 3*j2^2*k2*m4 + 6*i1^2*k2^2*m4 +
            4*j2*k2^2*m4 + 4*k2^3*m4 + 2*i1^4*l2*m4 + 4*i1^2*i2*l2*m4 + 5*i2^2*l2*m4 + 2*i1^2*j2*l2*m4 + 4*i2*j2*l2*m4 +
            5*j2^2*l2*m4 + 3*i1^2*k2*l2*m4 + 4*i2*k2*l2*m4 + j2*k2*l2*m4 + k2^2*l2*m4 + 3*i1^2*l2^2*m4 + 6*j2*l2^2*m4 +
            5*k2*l2^2*m4 + 5*l2^3*m4 + 2*i1^3*i3*m4 + 2*i1*i2*i3*m4 + 5*i1*j2*i3*m4 + 2*i1*k2*i3*m4 + 6*i1*l2*i3*m4 +
            i1*i2*j3*m4 + 5*i1*j2*j3*m4 + 6*i1*k2*j3*m4 + 6*i1*l2*j3*m4 + 5*j3^2*m4 + 3*i1*i2*l3*m4 + 6*i1*j2*l3*m4 +
            3*i1*k2*l3*m4 + 6*i3*l3*m4 + 5*j3*l3*m4 + 4*l3^2*m4 + 5*i1^5*i5 + 6*i1^3*k2*i5 + i1*k2^2*i5 + i1^3*l2*i5 +
            4*i1*k2*l2*i5 + 3*i1*l2^2*i5 + 3*i1^5*j5 + 4*i1^3*i2*j5 + 2*i1*i2^2*j5 + 4*i1^3*j2*j5 + 2*i1*i2*j2*j5 +
            3*i1*j2^2*j5 + 4*i1^3*k2*j5 + 2*i1*i2*k2*j5 + 5*i1*j2*k2*j5 + 3*i1*k2^2*j5 + 4*i1^3*l2*j5 + 4*i1*i2*l2*j5 +
            2*i1*j2*l2*j5 + 6*i1*k2*l2*j5 + i1*l2^2*j5 + i1^2*i3*j5 + 4*i2*i3*j5 + 3*j2*i3*j5 + 4*k2*i3*j5 + l2*i3*j5 +
            4*i1^2*j3*j5 + i2*j3*j5 + j2*j3*j5 + 4*k2*j3*j5 + 4*l2*j3*j5 + 2*i1^2*k3*j5 + 5*j2*k3*j5 + 2*k2*k3*j5 +
            2*l2*k3*j5 + 3*i1^2*l3*j5 + 3*j2*l3*j5 + 6*k2*l3*j5 + 3*l2*l3*j5 + 6*i1^3*i2*k5 + i1^3*j2*k5 + 6*i1^3*k2*k5 +
            6*i1*i2*k2*k5 + i1*j2*k2*k5 + 3*i1*k2^2*k5 + 6*i1^3*l2*k5 + 6*i1*i2*l2*k5 + i1*j2*l2*k5 + 4*i1*k2*l2*k5 +
            i1*l2^2*k5 + 4*i1^2*j3*k5 + 2*k2*j3*k5 + 6*l2*j3*k5 + 5*i1^2*k3*k5 + 2*k2*k3*k5 + 3*l2*k3*k5 + i1^2*l3*k5 +
            i2*l3*k5 + 6*j2*l3*k5 + 6*k2*l3*k5 + 6*l2*l3*k5 + i1*l4*k5 + 2*i1^5*l5 + 3*i1^3*i2*l5 + 5*i1*i2^2*l5 +
            4*i1^3*j2*l5 + 4*i1*i2*j2*l5 + 5*i1*j2^2*l5 + 5*i1^3*k2*l5 + 4*i1*i2*k2*l5 + i1*j2*k2*l5 + 3*i1^3*l2*l5 +
            5*i1*i2*l2*l5 + 5*i1*k2*l2*l5 + 4*i1*l2^2*l5 + 5*i1^2*i3*l5 + k2*i3*l5 + 6*l2*i3*l5 + 2*i1^2*j3*l5 +
            5*i2*j3*l5 + 2*j2*j3*l5 + 2*k2*j3*l5 + 2*l2*j3*l5 + 5*i1^2*k3*l5 + 5*i2*k3*l5 + 2*j2*k3*l5 + k2*k3*l5 +
            2*i1^2*l3*l5 + 4*i2*l3*l5 + 6*j2*l3*l5 + 5*k2*l3*l5 + 6*l2*l3*l5 + 6*i1*m4*l5 + 6*l5^2;
        Append(~JI, J10); Append(~Wght, 10);
    end if;
    if degmax le 10 then return JI, Wght; end if;

    /* Degree 11 */
    if not PrimaryOnly and degmin le 11 then
        J11:= 2*i1^11 + 2*i1^9*i2 + 4*i1^7*i2^2 + 3*i1^5*i2^3 + 4*i1^3*i2^4 + 3*i1^9*j2 + 4*i1^7*i2*j2 +
            3*i1^5*i2^2*j2 + 5*i1^3*i2^3*j2 + 4*i1^7*j2^2 + 6*i1^5*i2*j2^2 + 3*i1^3*i2^2*j2^2 + 2*i1^5*j2^3 +
            5*i1^3*i2*j2^3 + 4*i1^3*j2^4 + 2*i1^9*k2 + 6*i1^7*i2*k2 + i1^5*i2^2*k2 + 5*i1^3*i2^3*k2 + 3*i1*i2^4*k2 +
            i1^5*i2*j2*k2 + 6*i1^3*i2^2*j2*k2 + 2*i1*i2^3*j2*k2 + i1^5*j2^2*k2 + i1^3*i2*j2^2*k2 + 4*i1*i2^2*j2^2*k2 +
            2*i1^3*j2^3*k2 + 2*i1*i2*j2^3*k2 + 3*i1*j2^4*k2 + 2*i1^7*k2^2 + 4*i1^5*i2*k2^2 + i1^3*i2^2*k2^2 + i1^5*j2*k2^2
            + 2*i1^3*i2*j2*k2^2 + 3*i1*i2^2*j2*k2^2 + 2*i1^3*j2^2*k2^2 + i1*i2*j2^2*k2^2 + 3*i1*j2^3*k2^2 + 2*i1^5*k2^3 +
            6*i1^3*i2*k2^3 + 4*i1*i2^2*k2^3 + i1*i2*j2*k2^3 + i1^3*k2^4 + 3*i1*i2*k2^4 + i1*k2^5 + 4*i1^9*l2 + i1^7*i2*l2
            + 6*i1^5*i2^2*l2 + 3*i1^3*i2^3*l2 + 6*i1*i2^4*l2 + 2*i1^7*j2*l2 + 6*i1^5*i2*j2*l2 + 5*i1^3*i2^2*j2*l2 +
            4*i1*i2^3*j2*l2 + 6*i1^5*j2^2*l2 + 2*i1^3*i2*j2^2*l2 + i1*i2^2*j2^2*l2 + 4*i1^3*j2^3*l2 + 4*i1*i2*j2^3*l2 +
            6*i1*j2^4*l2 + 4*i1^5*i2*k2*l2 + 5*i1^3*i2^2*k2*l2 + 6*i1^5*j2*k2*l2 + 5*i1^3*i2*j2*k2*l2 + 3*i1^3*j2^2*k2*l2
            + 5*i1^5*k2^2*l2 + i1^3*j2*k2^2*l2 + 5*i1*i2*j2*k2^2*l2 + i1*j2^2*k2^2*l2 + 5*i1^3*k2^3*l2 + 6*i1*i2*k2^3*l2 +
            i1*j2*k2^3*l2 + 6*i1*k2^4*l2 + i1^7*l2^2 + 3*i1^5*i2*l2^2 + 6*i1^3*i2^2*l2^2 + 2*i1^5*j2*l2^2 +
            i1^3*i2*j2*l2^2 + 4*i1*i2^2*j2*l2^2 + 6*i1^3*j2^2*l2^2 + 6*i1*i2*j2^2*l2^2 + 4*i1*j2^3*l2^2 + 2*i1^5*k2*l2^2 +
            i1^3*i2*k2*l2^2 + 2*i1*i2^2*k2*l2^2 + 5*i1^3*j2*k2*l2^2 + 4*i1*i2*j2*k2*l2^2 + 5*i1*j2^2*k2*l2^2 +
            4*i1*i2*k2^2*l2^2 + 2*i1*j2*k2^2*l2^2 + 2*i1*k2^3*l2^2 + 5*i1^5*l2^3 + 5*i1^3*i2*l2^3 + 3*i1*i2^2*l2^3 +
            4*i1*j2^2*l2^3 + 2*i1^3*k2*l2^3 + 4*i1*i2*k2*l2^3 + 3*i1*j2*k2*l2^3 + 6*i1*k2^2*l2^3 + 6*i1^3*l2^4 +
            3*i1*i2*l2^4 + 2*i1*j2*l2^4 + 5*i1*k2*l2^4 + 4*i1^8*i3 + 4*i1^6*i2*i3 + 4*i1^4*i2^2*i3 + 4*i1^6*j2*i3 +
            6*i1^4*i2*j2*i3 + 4*i1^4*j2^2*i3 + 6*i1^6*k2*i3 + 2*i1^2*i2^2*k2*i3 + i1^4*j2*k2*i3 + 3*i1^2*i2*j2*k2*i3 +
            2*i1^2*j2^2*k2*i3 + 3*i1^4*k2^2*i3 + 4*i1^2*i2*k2^2*i3 + 6*i2^2*k2^2*i3 + 2*i1^2*j2*k2^2*i3 + 2*i2*j2*k2^2*i3
            + 6*j2^2*k2^2*i3 + 5*i1^2*k2^3*i3 + 3*i2*k2^3*i3 + j2*k2^3*i3 + 4*k2^4*i3 + 2*i1^6*l2*i3 + 4*i1^4*i2*l2*i3 +
            i1^2*i2^2*l2*i3 + 4*i1^4*j2*l2*i3 + 5*i1^2*i2*j2*l2*i3 + i1^2*j2^2*l2*i3 + 3*i1^4*k2*l2*i3 +
            5*i1^2*i2*k2*l2*i3 + 5*i2^2*k2*l2*i3 + 6*i1^2*j2*k2*l2*i3 + 4*i2*j2*k2*l2*i3 + 5*j2^2*k2*l2*i3 +
            i1^2*k2^2*l2*i3 + i2*k2^2*l2*i3 + 2*j2*k2^2*l2*i3 + 5*k2^3*l2*i3 + i1^4*l2^2*i3 + i1^2*i2*l2^2*i3 +
            3*i2^2*l2^2*i3 + i2*j2*l2^2*i3 + 3*j2^2*l2^2*i3 + 6*i1^2*k2*l2^2*i3 + 5*i2*k2*l2^2*i3 + 2*j2*k2*l2^2*i3 +
            4*i1^2*l2^3*i3 + 6*i2*l2^3*i3 + j2*l2^3*i3 + 6*l2^4*i3 + 2*i1^5*i3^2 + 3*i1*k2^2*i3^2 + 5*i1^3*l2*i3^2 +
            3*i1*k2*l2*i3^2 + 3*i1*l2^2*i3^2 + 3*i1^8*j3 + 6*i1^6*i2*j3 + 3*i1^4*i2^2*j3 + 2*i1^2*i2^3*j3 + 5*i1^6*j2*j3 +
            3*i1^4*i2*j2*j3 + i1^2*i2^2*j2*j3 + i1^4*j2^2*j3 + 6*i1^2*i2*j2^2*j3 + 5*i1^2*j2^3*j3 + 5*i1^6*k2*j3 +
            5*i1^4*i2*k2*j3 + 4*i1^2*i2^2*k2*j3 + 2*i2^3*k2*j3 + 5*i1^2*i2*j2*k2*j3 + i2^2*j2*k2*j3 + 5*i1^2*j2^2*k2*j3 +
            6*i2*j2^2*k2*j3 + 5*j2^3*k2*j3 + i1^4*k2^2*j3 + 3*i1^2*i2*k2^2*j3 + 4*i2^2*k2^2*j3 + 3*i1^2*j2*k2^2*j3 +
            5*i2*j2*k2^2*j3 + 5*j2^2*k2^2*j3 + 3*i1^2*k2^3*j3 + 6*i2*k2^3*j3 + 4*k2^4*j3 + 5*i1^6*l2*j3 + i1^2*i2^2*l2*j3
            + 6*i2^3*l2*j3 + 6*i1^4*j2*l2*j3 + 4*i1^2*i2*j2*l2*j3 + 3*i2^2*j2*l2*j3 + 2*i1^2*j2^2*l2*j3 + 4*i2*j2^2*l2*j3
            + j2^3*l2*j3 + i1^4*k2*l2*j3 + 4*i1^2*i2*k2*l2*j3 + i2^2*k2*l2*j3 + 2*i1^2*j2*k2*l2*j3 + 5*i2*j2*k2*l2*j3 +
            j2^2*k2*l2*j3 + 4*i2*k2^2*l2*j3 + 4*j2*k2^2*l2*j3 + 5*k2^3*l2*j3 + 2*i1^2*i2*l2^2*j3 + 5*i2^2*l2^2*j3 +
            2*i1^2*j2*l2^2*j3 + i2*j2*l2^2*j3 + j2^2*l2^2*j3 + 2*i1^2*k2*l2^2*j3 + 6*i2*k2*l2^2*j3 + 3*j2*k2*l2^2*j3 +
            6*k2^2*l2^2*j3 + 5*i1^2*l2^3*j3 + 6*i2*l2^3*j3 + j2*l2^3*j3 + 2*k2*l2^3*j3 + 3*l2^4*j3 + 3*i1^3*i2*i3*j3 +
            4*i1^3*j2*i3*j3 + 6*i1^3*k2*i3*j3 + 5*i1^3*l2*i3*j3 + 6*i1*i2*l2*i3*j3 + i1*j2*l2*i3*j3 + i1*k2*l2*i3*j3 +
            3*i1^5*j3^2 + 4*i1^3*i2*j3^2 + 5*i1*i2^2*j3^2 + 6*i1^3*j2*j3^2 + 4*i1*i2*j2*j3^2 + 5*i1*j2^2*j3^2 +
            3*i1^3*k2*j3^2 + 4*i1*j2*k2*j3^2 + 4*i1*k2^2*j3^2 + 3*i1^3*l2*j3^2 + 5*i1*i2*l2*j3^2 + 5*i1*j2*l2*j3^2 +
            4*i1*k2*l2*j3^2 + 4*i1*l2^2*j3^2 + i1^2*i3*j3^2 + 3*k2*i3*j3^2 + 6*l2*i3*j3^2 + i1^2*j3^3 + 2*i2*j3^3 +
            5*j2*j3^3 + 2*k2*j3^3 + 6*l2*j3^3 + 2*i1^8*k3 + 6*i1^6*i2*k3 + i1^2*i2^3*k3 + 2*i1^4*i2*j2*k3 +
            4*i1^2*i2^2*j2*k3 + 5*i1^4*j2^2*k3 + 3*i1^2*i2*j2^2*k3 + 6*i1^2*j2^3*k3 + 3*i1^6*k2*k3 + 6*i1^4*i2*k2*k3 +
            2*i1^2*i2^2*k2*k3 + 2*i2^3*k2*k3 + 6*i1^4*j2*k2*k3 + 4*i1^2*i2*j2*k2*k3 + i2^2*j2*k2*k3 + i1^2*j2^2*k2*k3 +
            6*i2*j2^2*k2*k3 + 5*j2^3*k2*k3 + 2*i1^4*k2^2*k3 + i1^2*i2*k2^2*k3 + 6*i2^2*k2^2*k3 + 3*i1^2*j2*k2^2*k3 +
            i2*j2*k2^2*k3 + 6*i1^2*k2^3*k3 + 2*j2*k2^3*k3 + 3*k2^4*k3 + 5*i1^6*l2*k3 + 4*i1^4*i2*l2*k3 + 6*i1^2*i2^2*l2*k3
            + i2^3*l2*k3 + 3*i1^4*j2*l2*k3 + 5*i1^2*i2*j2*l2*k3 + 4*i2^2*j2*l2*k3 + 3*i1^2*j2^2*l2*k3 + 3*i2*j2^2*l2*k3 +
            6*j2^3*l2*k3 + 5*i1^4*k2*l2*k3 + 6*i1^2*i2*k2*l2*k3 + 2*i2^2*k2*l2*k3 + 4*i1^2*j2*k2*l2*k3 + i2*j2*k2*l2*k3 +
            4*j2^2*k2*l2*k3 + 6*i1^2*k2^2*l2*k3 + 6*i2*k2^2*l2*k3 + 5*j2*k2^2*l2*k3 + 5*k2^3*l2*k3 + 2*i1^2*i2*l2^2*k3 +
            i2^2*l2^2*k3 + 5*i1^2*j2*l2^2*k3 + 5*i2*j2*l2^2*k3 + j2^2*l2^2*k3 + i1^2*k2*l2^2*k3 + 5*i2*k2*l2^2*k3 +
            4*k2^2*l2^2*k3 + i1^2*l2^3*k3 + 5*k2*l2^3*k3 + 4*l2^4*k3 + i1^5*j3*k3 + i1*i2^2*j3*k3 + 3*i1^3*j2*j3*k3 +
            5*i1*i2*j2*j3*k3 + i1*j2^2*j3*k3 + 5*i1^3*k2*j3*k3 + 6*i1*i2*k2*j3*k3 + 6*i1*j2*k2*j3*k3 + 5*i1*k2^2*j3*k3 +
            2*i1^3*l2*j3*k3 + 3*i1*j2*l2*j3*k3 + 5*i1*k2*l2*j3*k3 + 5*i1^2*j3^2*k3 + 6*i2*j3^2*k3 + j2*j3^2*k3 +
            5*l2*j3^2*k3 + 6*i1^3*i2*k3^2 + 5*i1*i2^2*k3^2 + 3*i1^3*j2*k3^2 + 4*i1*i2*j2*k3^2 + 5*i1*j2^2*k3^2 +
            i1^3*k2*k3^2 + 6*i1*i2*k2*k3^2 + 5*i1*j2*k2*k3^2 + 2*i1*k2^2*k3^2 + i1^3*l2*k3^2 + 3*i1*i2*l2*k3^2 +
            5*i1*l2^2*k3^2 + i2*j3*k3^2 + 6*j2*j3*k3^2 + 4*k2*j3*k3^2 + 2*i1^2*k3^3 + 2*k2*k3^3 + 5*i1^8*l3 + 3*i1^6*i2*l3
            + 5*i1^4*i2^2*l3 + i1^2*i2^3*l3 + 4*i2^4*l3 + i1^6*j2*l3 + i1^4*i2*j2*l3 + 2*i1^2*i2^2*j2*l3 + 5*i2^3*j2*l3 +
            5*i1^4*j2^2*l3 + 3*i2^2*j2^2*l3 + 4*i1^2*j2^3*l3 + 5*i2*j2^3*l3 + 4*j2^4*l3 + 6*i1^6*k2*l3 + 5*i1^4*i2*k2*l3 +
            i1^2*i2^2*k2*l3 + 5*i2^3*k2*l3 + 6*i1^4*j2*k2*l3 + i1^2*i2*j2*k2*l3 + i2^2*j2*k2*l3 + 6*i1^2*j2^2*k2*l3 +
            4*i2*j2^2*k2*l3 + 4*j2^3*k2*l3 + 2*i1^4*k2^2*l3 + i1^2*i2*k2^2*l3 + 6*i1^2*j2*k2^2*l3 + i1^2*k2^3*l3 +
            6*i2*k2^3*l3 + 6*j2*k2^3*l3 + 2*k2^4*l3 + i1^6*l2*l3 + 6*i1^4*i2*l2*l3 + 2*i1^2*i2^2*l2*l3 + 6*i2^3*l2*l3 +
            4*i1^4*j2*l2*l3 + 2*i1^2*i2*j2*l2*l3 + 6*i2^2*j2*l2*l3 + 3*i1^2*j2^2*l2*l3 + 5*i2*j2^2*l2*l3 + 4*j2^3*l2*l3 +
            4*i1^4*k2*l2*l3 + i1^2*i2*k2*l2*l3 + 6*i2^2*k2*l2*l3 + 2*i1^2*j2*k2*l2*l3 + 2*j2^2*k2*l2*l3 +
            4*i1^2*k2^2*l2*l3 + i2*k2^2*l2*l3 + 2*j2*k2^2*l2*l3 + 2*k2^3*l2*l3 + i1^4*l2^2*l3 + 2*i1^2*i2*l2^2*l3 +
            2*i1^2*j2*l2^2*l3 + 5*i2*j2*l2^2*l3 + i1^2*k2*l2^2*l3 + i2*k2*l2^2*l3 + 5*j2*k2*l2^2*l3 + 5*k2^2*l2^2*l3 +
            5*i1^2*l2^3*l3 + 4*j2*l2^3*l3 + 3*k2*l2^3*l3 + 3*i1^5*i3*l3 + 5*i1*i2^2*i3*l3 + i1^3*j2*i3*l3 +
            4*i1*i2*j2*i3*l3 + 5*i1*j2^2*i3*l3 + 6*i1^3*k2*i3*l3 + 6*i1*i2*k2*i3*l3 + 3*i1*j2*k2*i3*l3 + 4*i1*k2^2*i3*l3 +
            2*i1^3*l2*i3*l3 + 6*i1*i2*l2*i3*l3 + 5*i1*j2*l2*i3*l3 + 6*i1*k2*l2*i3*l3 + 5*i1*l2^2*i3*l3 + 2*i1^2*i3^2*l3 +
            k2*i3^2*l3 + 2*l2*i3^2*l3 + 4*i1^5*j3*l3 + 6*i1^3*i2*j3*l3 + 4*i1*i2^2*j3*l3 + 3*i1^3*j2*j3*l3 +
            5*i1*i2*j2*j3*l3 + 5*i1*j2^2*j3*l3 + i1^3*k2*j3*l3 + 6*i1*j2*k2*j3*l3 + i1*k2^2*j3*l3 + 6*i1^3*l2*j3*l3 +
            3*i1*i2*l2*j3*l3 + 3*i1*j2*l2*j3*l3 + 4*i1*l2^2*j3*l3 + i1^2*j3^2*l3 + 3*i2*j3^2*l3 + j2*j3^2*l3 + 2*i1^5*l3^2
            + 6*i1*i2^2*l3^2 + i1^3*j2*l3^2 + 4*i1*i2*j2*l3^2 + 3*i1*j2^2*l3^2 + 5*i1^3*k2*l3^2 + i1*i2*k2*l3^2 +
            6*i1*k2^2*l3^2 + 4*i1^3*l2*l3^2 + 4*i1*i2*l2*l3^2 + 5*i1*j2*l2*l3^2 + 5*i1*k2*l2*l3^2 + 2*i1*l2^2*l3^2 +
            6*i1^2*i3*l3^2 + 3*i2*i3*l3^2 + j2*i3*l3^2 + 4*l2*i3*l3^2 + 5*i1^2*j3*l3^2 + 2*i2*j3*l3^2 + 6*j2*j3*l3^2 +
            i1^2*l3^3 + 6*i2*l3^3 + 5*j2*l3^3 + 5*k2*l3^3 + l2*l3^3 + 6*i1^5*i2*i4 + 2*i1^3*i2^2*i4 + i1*i2^3*i4 +
            4*i1^5*j2*i4 + 3*i1^3*i2*j2*i4 + 4*i1*i2^2*j2*i4 + 2*i1^3*j2^2*i4 + 3*i1*i2*j2^2*i4 + 6*i1*j2^3*i4 +
            4*i1^5*k2*i4 + 6*i1^3*i2*k2*i4 + 5*i1*i2^2*k2*i4 + 5*i1^3*j2*k2*i4 + i1*i2*j2*k2*i4 + i1*j2^2*k2*i4 +
            i1*i2*k2^2*i4 + 6*i1*j2*k2^2*i4 + i1*k2^3*i4 + 2*i1^5*l2*i4 + 5*i1^3*i2*l2*i4 + 2*i1*i2^2*l2*i4 +
            3*i1^3*j2*l2*i4 + 6*i1*i2*j2*l2*i4 + 6*i1*j2^2*l2*i4 + i1^3*k2*l2*i4 + 4*i1*j2*k2*l2*i4 + 2*i1^3*l2^2*i4 +
            3*i1*i2*l2^2*i4 + 3*i1*j2*l2^2*i4 + 4*i1*k2*l2^2*i4 + 2*i1^4*i3*i4 + 2*i1^2*i2*i3*i4 + 5*i1^2*j2*i3*i4 +
            3*i1^2*k2*i3*i4 + 2*i2*k2*i3*i4 + 5*j2*k2*i3*i4 + 6*k2^2*i3*i4 + 2*i2*l2*i3*i4 + 5*j2*l2*i3*i4 + 2*k2*l2*i3*i4
            + 3*l2^2*i3*i4 + 5*i1^4*j3*i4 + 2*i1^2*i2*j3*i4 + 2*i2^2*j3*i4 + 4*i1^2*j2*j3*i4 + 3*i2*j2*j3*i4 +
            2*j2^2*j3*i4 + 3*i1^2*k2*j3*i4 + i2*k2*j3*i4 + 5*j2*k2*j3*i4 + 3*k2^2*j3*i4 + 4*i1^2*l2*j3*i4 + 5*i2*l2*j3*i4
            + 4*j2*l2*j3*i4 + 2*k2*l2*j3*i4 + 3*l2^2*j3*i4 + 2*i1*i3*j3*i4 + 3*i1*j3^2*i4 + 3*i1^4*k3*i4 + 4*i1^2*i2*k3*i4
            + 4*i2^2*k3*i4 + 2*i1^2*j2*k3*i4 + 6*i2*j2*k3*i4 + 4*j2^2*k3*i4 + 2*i1^2*k2*k3*i4 + 4*i2*k2*k3*i4 +
            4*j2*k2*k3*i4 + 5*k2^2*k3*i4 + 2*i1^2*l2*k3*i4 + 2*j2*l2*k3*i4 + 3*k2*l2*k3*i4 + 2*l2^2*k3*i4 + 3*i1*j3*k3*i4
            + i1^7*j4 + 5*i1^5*i2*j4 + i1^3*i2^2*j4 + 5*i1^5*j2*j4 + 5*i1^3*i2*j2*j4 + i1^3*j2^2*j4 + 5*i1^5*k2*j4 +
            3*i1^3*i2*k2*j4 + 6*i1*i2^2*k2*j4 + 3*i1^3*j2*k2*j4 + 2*i1*i2*j2*k2*j4 + 6*i1*j2^2*k2*j4 + 2*i1^3*k2^2*j4 +
            i1*i2*k2^2*j4 + 2*i1*j2*k2^2*j4 + 4*i1^5*l2*j4 + 6*i1^3*i2*l2*j4 + 4*i1*i2^2*l2*j4 + 6*i1^3*j2*l2*j4 +
            6*i1*i2*j2*l2*j4 + 4*i1*j2^2*l2*j4 + 2*i1^3*k2*l2*j4 + 6*i1*i2*k2*l2*j4 + 6*i1*j2*k2*l2*j4 + 3*i1*k2^2*l2*j4 +
            2*i1^3*l2^2*j4 + i1*i2*l2^2*j4 + 2*i1*j2*l2^2*j4 + 4*i1*k2*l2^2*j4 + 4*i1*l2^3*j4 + 6*i1^2*i2*j3*j4 +
            i1^2*j2*j3*j4 + 4*i1^2*k2*j3*j4 + 6*i2*k2*j3*j4 + j2*k2*j3*j4 + 6*k2^2*j3*j4 + 6*i1^2*l2*j3*j4 + i2*l2*j3*j4 +
            6*j2*l2*j3*j4 + 5*k2*l2*j3*j4 + l2^2*j3*j4 + 2*i1*j3^2*j4 + 5*i1^7*k4 + 4*i1^5*i2*k4 + 3*i1^5*j2*k4 +
            4*i1^5*k2*k4 + i1^3*i2*k2*k4 + 6*i1^3*j2*k2*k4 + i1^3*k2^2*k4 + 2*i1*i2*k2^2*k4 + 5*i1*j2*k2^2*k4 +
            4*i1*k2^3*k4 + 5*i1^5*l2*k4 + 4*i1*i2*k2*l2*k4 + 3*i1*j2*k2*l2*k4 + i1*k2^2*l2*k4 + i1^3*l2^2*k4 +
            6*i1*i2*l2^2*k4 + i1*j2*l2^2*k4 + 4*i1*k2*l2^2*k4 + 2*i1*l2^3*k4 + 2*i1^4*k3*k4 + 4*i1^2*k2*k3*k4 +
            6*k2^2*k3*k4 + 2*i1^2*l2*k3*k4 + 3*k2*l2*k3*k4 + 2*l2^2*k3*k4 + 3*i1^7*l4 + i1^5*i2*l4 + 2*i1^3*i2^2*l4 +
            5*i1*i2^3*l4 + 3*i1^5*j2*l4 + 2*i1^3*i2*j2*l4 + 6*i1*i2^2*j2*l4 + 3*i1^3*j2^2*l4 + i1*i2*j2^2*l4 +
            2*i1*j2^3*l4 + 4*i1^5*k2*l4 + i1^3*i2*k2*l4 + 6*i1*i2^2*k2*l4 + 5*i1^3*j2*k2*l4 + 4*i1*i2*j2*k2*l4 +
            4*i1*j2^2*k2*l4 + 3*i1^3*k2^2*l4 + 3*i1*k2^3*l4 + 3*i1^5*l2*l4 + 6*i1^3*i2*l2*l4 + 2*i1*i2^2*l2*l4 +
            2*i1^3*j2*l2*l4 + i1*i2*j2*l2*l4 + 4*i1*j2^2*l2*l4 + 5*i1^3*k2*l2*l4 + 3*i1*i2*k2*l2*l4 + 5*i1*j2*k2*l2*l4 +
            i1*k2^2*l2*l4 + 4*i1^3*l2^2*l4 + 5*i1*i2*l2^2*l4 + i1*j2*l2^2*l4 + i1*k2*l2^2*l4 + 4*i1*l2^3*l4 + 2*i1^4*j3*l4
            + 2*i1^2*i2*j3*l4 + 6*i2^2*j3*l4 + 4*i1^2*j2*j3*l4 + 2*i2*j2*j3*l4 + 6*j2^2*j3*l4 + 5*i1^2*k2*j3*l4 +
            5*i2*k2*j3*l4 + 3*j2*k2*j3*l4 + 3*k2^2*j3*l4 + 5*i1^2*l2*j3*l4 + 5*i2*l2*j3*l4 + 3*j2*l2*j3*l4 + 2*k2*l2*j3*l4
            + 2*l2^2*j3*l4 + 4*i1^4*k3*l4 + 4*i1^2*i2*k3*l4 + 4*i1^2*j2*k3*l4 + i1^2*k2*k3*l4 + 2*i2*k2*k3*l4 +
            2*j2*k2*k3*l4 + 6*k2^2*k3*l4 + 6*i1^2*l2*k3*l4 + 4*i2*l2*k3*l4 + 5*j2*l2*k3*l4 + 6*k2*l2*k3*l4 + 2*l2^2*k3*l4
            + 5*i1*j3*k3*l4 + 5*i1*k3^2*l4 + 2*i1^4*l3*l4 + 2*i1^2*i2*l3*l4 + 2*i1^2*j2*l3*l4 + 5*i2*j2*l3*l4 +
            2*j2^2*l3*l4 + i1^2*k2*l3*l4 + i2*k2*l3*l4 + 4*j2*k2*l3*l4 + 3*i1^2*l2*l3*l4 + 4*i2*l2*l3*l4 + 6*j2*l2*l3*l4 +
            3*k2*l2*l3*l4 + l2^2*l3*l4 + 4*i1*j3*l3*l4 + 3*i1^3*l4^2 + 5*i1*i2*l4^2 + 6*i1*j2*l4^2 + i1*k2*l4^2 +
            3*i1*l2*l4^2 + j3*l4^2 + 5*i1^7*m4 + 6*i1^5*i2*m4 + i1^3*i2^2*m4 + 5*i1*i2^3*m4 + 4*i1^5*j2*m4 +
            5*i1^3*i2*j2*m4 + 6*i1*i2^2*j2*m4 + i1^3*j2^2*m4 + i1*i2*j2^2*m4 + 2*i1*j2^3*m4 + 3*i1^5*k2*m4 +
            3*i1^3*i2*k2*m4 + 6*i1*i2^2*k2*m4 + 6*i1*i2*j2*k2*m4 + 2*i1*j2^2*k2*m4 + 3*i1^3*k2^2*m4 + 3*i1*i2*k2^2*m4 +
            5*i1*j2*k2^2*m4 + 5*i1*k2^3*m4 + 3*i1^5*l2*m4 + 4*i1^3*i2*l2*m4 + 2*i1*i2^2*l2*m4 + 4*i1^3*j2*l2*m4 +
            4*i1*i2*j2*l2*m4 + i1*j2^2*l2*m4 + 3*i1^3*k2*l2*m4 + 3*i1*i2*k2*l2*m4 + 3*i1*j2*k2*l2*m4 + 3*i1*k2^2*l2*m4 +
            4*i1^3*l2^2*m4 + 6*i1*j2*l2^2*m4 + 2*i1*k2*l2^2*m4 + 2*i1*l2^3*m4 + 4*i1^4*i3*m4 + i1^2*i2*i3*m4 +
            6*i1^2*j2*i3*m4 + 3*i1^2*k2*i3*m4 + 3*i2*k2*i3*m4 + 4*j2*k2*i3*m4 + 5*k2^2*i3*m4 + 4*i1^2*l2*i3*m4 +
            6*i2*l2*i3*m4 + j2*l2*i3*m4 + 6*k2*l2*i3*m4 + l2^2*i3*m4 + 4*i1^4*j3*m4 + 5*i1^2*i2*j3*m4 + 3*i2^2*j3*m4 +
            3*i1^2*j2*j3*m4 + i2*j2*j3*m4 + 3*j2^2*j3*m4 + 6*j2*k2*j3*m4 + 3*k2^2*j3*m4 + i1^2*l2*j3*m4 + 6*i2*l2*j3*m4 +
            2*k2*l2*j3*m4 + 6*l2^2*j3*m4 + 4*i1*j3^2*m4 + 5*i1^4*l3*m4 + 6*i1^2*i2*l3*m4 + 3*i2^2*l3*m4 + i1^2*j2*l3*m4 +
            5*i2*j2*l3*m4 + 6*j2^2*l3*m4 + 2*i1^2*k2*l3*m4 + 5*i2*k2*l3*m4 + 6*j2*k2*l3*m4 + 3*k2^2*l3*m4 +
            6*i1^2*l2*l3*m4 + i2*l2*l3*m4 + 3*j2*l2*l3*m4 + 3*k2*l2*l3*m4 + 3*l2^2*l3*m4 + i1*i3*l3*m4 + 4*i1*j3*l3*m4 +
            5*i1*l3^2*m4 + 5*l3*m4^2 + 6*i1^6*i5 + 5*i1^4*k2*i5 + 3*i1^2*k2^2*i5 + 5*i1^4*l2*i5 + i1^2*k2*l2*i5 +
            k2^2*l2*i5 + 4*i1^2*l2^2*i5 + 5*k2*l2^2*i5 + 2*l2^3*i5 + 2*i1^6*j5 + i1^4*i2*j5 + 3*i1^2*i2^2*j5 + 3*i2^3*j5 +
            i1^4*j2*j5 + 2*i1^2*i2*j2*j5 + 5*i2^2*j2*j5 + 2*i1^2*j2^2*j5 + 2*i2*j2^2*j5 + 4*j2^3*j5 + i1^2*i2*k2*j5 +
            i2^2*k2*j5 + 5*i1^2*j2*k2*j5 + 2*i2*j2*k2*j5 + 4*j2^2*k2*j5 + 4*i1^2*k2^2*j5 + 2*j2*k2^2*j5 + 3*k2^3*j5 +
            6*i1^4*l2*j5 + 3*i2^2*l2*j5 + 5*i1^2*j2*l2*j5 + 2*i2*j2*l2*j5 + 2*j2^2*l2*j5 + 3*i1^2*k2*l2*j5 + 5*i2*k2*l2*j5
            + j2*k2*l2*j5 + 6*k2^2*l2*j5 + 2*i1^2*l2^2*j5 + 5*j2*l2^2*j5 + 4*k2*l2^2*j5 + 5*l2^3*j5 + i1^3*i3*j5 +
            4*i1*k2*i3*j5 + 5*i1*l2*i3*j5 + i1^3*j3*j5 + 5*i1*j2*j3*j5 + 3*i1*k2*j3*j5 + 2*i1*l2*j3*j5 + 2*i3*j3*j5 +
            2*j3^2*j5 + 5*i1^3*k3*j5 + 6*i1*i2*k3*j5 + 4*i1*j2*k3*j5 + i1*k2*k3*j5 + 5*i1*l2*k3*j5 + 2*k3^2*j5 +
            6*i1^3*l3*j5 + 3*i1*i2*l3*j5 + 4*i1*j2*l3*j5 + 2*i1*k2*l3*j5 + 6*i1*l2*l3*j5 + i3*l3*j5 + 5*i1^2*i4*j5 +
            5*i2*i4*j5 + 4*j2*i4*j5 + 6*k2*i4*j5 + 5*l2*i4*j5 + 6*i1^2*l4*j5 + 6*i2*l4*j5 + 5*j2*l4*j5 + 2*k2*l4*j5 +
            6*l2*l4*j5 + 6*i1^6*k5 + 2*i1^4*i2*k5 + 5*i1^4*j2*k5 + 6*i1^4*k2*k5 + 5*i1^2*i2*k2*k5 + 2*i1^2*j2*k2*k5 +
            4*i1^2*k2^2*k5 + 3*i2*k2^2*k5 + 4*j2*k2^2*k5 + k2^3*k5 + 3*i1^4*l2*k5 + 6*i1^2*i2*l2*k5 + i1^2*j2*l2*k5 +
            4*i1^2*k2*l2*k5 + k2^2*l2*k5 + 6*i1^2*l2^2*k5 + i2*l2^2*k5 + 6*j2*l2^2*k5 + 5*l2^3*k5 + 4*i1^3*j3*k5 +
            2*i1*k2*j3*k5 + i1*l2*j3*k5 + 4*i1^3*k3*k5 + 5*i1*k2*k3*k5 + 5*i1*l2*k3*k5 + 3*i1^3*l3*k5 + 6*i1*i2*l3*k5 +
            i1*j2*l3*k5 + 2*i1*k2*l3*k5 + 3*i1*l2*l3*k5 + j3*l3*k5 + 5*i1^2*l4*k5 + 6*k2*l4*k5 + 4*l2*l4*k5 + 4*i1^4*i2*l5
            + 4*i1^2*i2^2*l5 + 5*i1^4*j2*l5 + 6*i1^2*i2*j2*l5 + 4*i1^2*j2^2*l5 + 6*i1^4*k2*l5 + 3*i1^2*i2*k2*l5 +
            4*i2^2*k2*l5 + 2*i1^2*j2*k2*l5 + 6*i2*j2*k2*l5 + 4*j2^2*k2*l5 + i1^2*k2^2*l5 + 2*i2*k2^2*l5 + 4*j2*k2^2*l5 +
            5*k2^3*l5 + 5*i1^4*l2*l5 + i2^2*l2*l5 + 2*i1^2*j2*l2*l5 + 5*i2*j2*l2*l5 + j2^2*l2*l5 + 3*i1^2*k2*l2*l5 +
            5*i2*k2*l2*l5 + j2*k2*l2*l5 + 5*k2^2*l2*l5 + 5*i1^2*l2^2*l5 + 6*i2*l2^2*l5 + 5*k2*l2^2*l5 + 6*l2^3*l5 +
            5*i1*k2*i3*l5 + 4*i1*l2*i3*l5 + 5*i1^3*j3*l5 + 5*i1*i2*j3*l5 + 2*i1*j2*j3*l5 + 5*i1*k2*j3*l5 + 5*i1*l2*j3*l5 +
            2*j3^2*l5 + i1^3*k3*l5 + 4*i1*i2*k3*l5 + 3*i1*j2*k3*l5 + 4*i1*k2*k3*l5 + 3*i1*l2*k3*l5 + 4*j3*k3*l5 +
            4*i1^3*l3*l5 + 4*i1*i2*l3*l5 + 3*i1*k2*l3*l5 + 2*i1*l2*l3*l5 + 4*i3*l3*l5 + 3*j3*l3*l5 + 6*i1^2*m4*l5 +
            5*i2*m4*l5 + 2*j2*m4*l5 + k2*m4*l5 + 6*l2*m4*l5 + 2*i1*l5^2;
        Append(~JI, J11); Append(~Wght, 11);
    end if;
    if degmax le 11 then return JI, Wght; end if;

    /* Degree 13 */
    if not PrimaryOnly and degmin le 13 then
        J13:= 3*i1^13 + 3*i1^11*i2 + 6*i1^9*i2^2 + 4*i1^7*i2^3 + 2*i1^5*i2^4 + 2*i1^3*i2^5 + 3*i1^11*j2 + i1^9*i2*j2
            + 5*i1^7*i2^2*j2 + 2*i1^5*i2^3*j2 + 4*i1^3*i2^4*j2 + 3*i1^9*j2^2 + i1^7*i2*j2^2 + 3*i1^5*i2^2*j2^2 +
            6*i1^3*i2^3*j2^2 + 4*i1^7*j2^3 + i1^5*i2*j2^3 + i1^3*i2^2*j2^3 + 6*i1^5*j2^4 + 3*i1^3*i2*j2^4 + 5*i1^3*j2^5 +
            6*i1^11*k2 + 3*i1^9*i2*k2 + i1^7*i2^2*k2 + i1^5*i2^3*k2 + 6*i1^3*i2^4*k2 + 5*i1*i2^5*k2 + 4*i1^9*j2*k2 +
            i1^7*i2*j2*k2 + i1^5*i2^2*j2*k2 + 4*i1^3*i2^3*j2*k2 + 3*i1*i2^4*j2*k2 + 3*i1^7*j2^2*k2 + 6*i1^5*i2*j2^2*k2 +
            i1^3*i2^2*j2^2*k2 + i1*i2^3*j2^2*k2 + 6*i1^5*j2^3*k2 + 4*i1^3*i2*j2^3*k2 + 6*i1*i2^2*j2^3*k2 + 6*i1^3*j2^4*k2
            + 4*i1*i2*j2^4*k2 + 2*i1*j2^5*k2 + 2*i1^9*k2^2 + 5*i1^7*i2*k2^2 + 4*i1^5*i2^2*k2^2 + 4*i1^3*i2^3*k2^2 +
            2*i1*i2^4*k2^2 + 6*i1^7*j2*k2^2 + 4*i1^5*i2*j2*k2^2 + i1^3*i2^2*j2*k2^2 + 5*i1*i2^3*j2*k2^2 + 2*i1^5*j2^2*k2^2
            + 2*i1^3*i2*j2^2*k2^2 + i1*i2^2*j2^2*k2^2 + 3*i1*i2*j2^3*k2^2 + 3*i1*j2^4*k2^2 + 4*i1^5*i2*k2^3 +
            5*i1*i2^3*k2^3 + 3*i1^5*j2*k2^3 + 2*i1^3*i2*j2*k2^3 + 2*i1*i2^2*j2*k2^3 + 6*i1^3*j2^2*k2^3 + 4*i1*i2*j2^2*k2^3
            + 3*i1*j2^3*k2^3 + i1^5*k2^4 + i1^3*i2*k2^4 + i1*i2^2*k2^4 + 5*i1*i2*j2*k2^4 + 2*i1*j2^2*k2^4 + 3*i1*k2^6 +
            2*i1^11*l2 + 2*i1^9*i2*l2 + 2*i1^7*i2^2*l2 + 5*i1^5*i2^3*l2 + 3*i1^3*i2^4*l2 + 3*i1*i2^5*l2 + 2*i1^9*j2*l2 +
            3*i1^7*i2*j2*l2 + 3*i1^5*i2^2*j2*l2 + 2*i1^3*i2^3*j2*l2 + 6*i1*i2^4*j2*l2 + 3*i1^7*j2^2*l2 + 3*i1^5*i2*j2^2*l2
            + 4*i1^3*i2^2*j2^2*l2 + 2*i1*i2^3*j2^2*l2 + 3*i1^5*j2^3*l2 + 2*i1^3*i2*j2^3*l2 + 5*i1*i2^2*j2^3*l2 +
            3*i1^3*j2^4*l2 + i1*i2*j2^4*l2 + 4*i1*j2^5*l2 + 3*i1^9*k2*l2 + 6*i1^7*i2*k2*l2 + 2*i1^5*i2^2*k2*l2 +
            4*i1*i2^4*k2*l2 + 4*i1^7*j2*k2*l2 + 5*i1^3*i2^2*j2*k2*l2 + 5*i1*i2^3*j2*k2*l2 + 5*i1^3*i2*j2^2*k2*l2 +
            3*i1*i2^2*j2^2*k2*l2 + 4*i1^3*j2^3*k2*l2 + 5*i1*i2*j2^3*k2*l2 + 4*i1*j2^4*k2*l2 + i1^7*k2^2*l2 +
            2*i1^5*i2*k2^2*l2 + 5*i1^3*i2^2*k2^2*l2 + 6*i1*i2^3*k2^2*l2 + 3*i1^5*j2*k2^2*l2 + 6*i1^3*i2*j2*k2^2*l2 +
            4*i1*i2^2*j2*k2^2*l2 + 3*i1^3*j2^2*k2^2*l2 + 3*i1*i2*j2^2*k2^2*l2 + i1*j2^3*k2^2*l2 + i1^5*k2^3*l2 +
            3*i1^3*i2*k2^3*l2 + 4*i1*i2^2*k2^3*l2 + 2*i1^3*j2*k2^3*l2 + 5*i1*j2^2*k2^3*l2 + 6*i1^3*k2^4*l2 +
            5*i1*i2*k2^4*l2 + i1*j2*k2^4*l2 + 3*i1*k2^5*l2 + i1^9*l2^2 + 3*i1^7*i2*l2^2 + 3*i1^5*i2^2*l2^2 +
            6*i1^3*i2^3*l2^2 + 4*i1*i2^4*l2^2 + 5*i1^7*j2*l2^2 + i1^5*i2*j2*l2^2 + 6*i1^3*i2^2*j2*l2^2 + 6*i1*i2^3*j2*l2^2
            + 6*i1^3*i2*j2^2*l2^2 + 3*i1^3*j2^3*l2^2 + i1*i2*j2^3*l2^2 + 3*i1*j2^4*l2^2 + 2*i1^7*k2*l2^2 +
            5*i1^5*i2*k2*l2^2 + 2*i1^3*i2^2*k2*l2^2 + 6*i1*i2^3*k2*l2^2 + 2*i1^5*j2*k2*l2^2 + i1^3*i2*j2*k2*l2^2 +
            4*i1*i2^2*j2*k2*l2^2 + 6*i1^3*j2^2*k2*l2^2 + 5*i1*i2*j2^2*k2*l2^2 + 6*i1*j2^3*k2*l2^2 + i1^5*k2^2*l2^2 +
            4*i1^3*i2*k2^2*l2^2 + 2*i1*i2^2*k2^2*l2^2 + 4*i1^3*j2*k2^2*l2^2 + i1*i2*j2*k2^2*l2^2 + 5*i1*j2^2*k2^2*l2^2 +
            4*i1^3*k2^3*l2^2 + 3*i1*i2*k2^3*l2^2 + 5*i1*j2*k2^3*l2^2 + 5*i1^5*i2*l2^3 + 2*i1^3*i2^2*l2^3 + 4*i1*i2^3*l2^3
            + 2*i1^5*j2*l2^3 + i1^3*i2*j2*l2^3 + 2*i1*i2^2*j2*l2^3 + 3*i1^3*j2^2*l2^3 + 5*i1*i2*j2^2*l2^3 + 3*i1*j2^3*l2^3
            + i1^5*k2*l2^3 + 5*i1*i2^2*k2*l2^3 + 3*i1^3*j2*k2*l2^3 + 5*i1*i2*j2*k2*l2^3 + 5*i1*j2^2*k2*l2^3 +
            3*i1^3*k2^2*l2^3 + i1*i2*k2^2*l2^3 + 2*i1*k2^3*l2^3 + 3*i1^5*l2^4 + 6*i1^3*i2*l2^4 + 4*i1*i2^2*l2^4 +
            6*i1^3*j2*l2^4 + 3*i1*j2^2*l2^4 + 4*i1^3*k2*l2^4 + i1*i2*k2*l2^4 + 3*i1*j2*k2*l2^4 + 2*i1*k2^2*l2^4 +
            3*i1^3*l2^5 + 3*i1*i2*l2^5 + 2*i1*j2*l2^5 + i1*k2*l2^5 + 4*i1*l2^6 + 2*i1^10*i3 + 2*i1^8*i2*i3 + i1^6*i2^2*i3
            + i1^4*i2^3*i3 + 4*i1^6*i2*j2*i3 + 4*i1^4*i2^2*j2*i3 + 2*i1^6*j2^2*i3 + 3*i1^4*i2*j2^2*i3 + 6*i1^4*j2^3*i3 +
            5*i1^8*k2*i3 + 5*i1^6*i2*k2*i3 + 4*i1^4*i2^2*k2*i3 + 4*i1^2*i2^3*k2*i3 + 6*i1^6*j2*k2*i3 + 5*i1^4*i2*j2*k2*i3
            + 2*i1^2*i2^2*j2*k2*i3 + 5*i1^4*j2^2*k2*i3 + 5*i1^2*i2*j2^2*k2*i3 + 3*i1^2*j2^3*k2*i3 + i1^6*k2^2*i3 +
            2*i1^4*i2*k2^2*i3 + 3*i1^2*i2^2*k2^2*i3 + 5*i2^3*k2^2*i3 + 2*i1^2*i2*j2*k2^2*i3 + 6*i2^2*j2*k2^2*i3 +
            2*i1^2*j2^2*k2^2*i3 + i2*j2^2*k2^2*i3 + 2*j2^3*k2^2*i3 + 5*i1^4*k2^3*i3 + 6*i1^2*i2*k2^3*i3 + 6*i2^2*k2^3*i3 +
            3*i1^2*j2*k2^3*i3 + 5*i2*j2*k2^3*i3 + 3*j2^2*k2^3*i3 + 5*i1^2*k2^4*i3 + 6*i2*k2^4*i3 + 5*j2*k2^4*i3 +
            6*k2^5*i3 + 3*i1^8*l2*i3 + 4*i1^4*i2^2*l2*i3 + 2*i1^2*i2^3*l2*i3 + 6*i1^6*j2*l2*i3 + 5*i1^4*i2*j2*l2*i3 +
            i1^2*i2^2*j2*l2*i3 + 5*i1^4*j2^2*l2*i3 + 6*i1^2*i2*j2^2*l2*i3 + 5*i1^2*j2^3*l2*i3 + i1^4*i2*k2*l2*i3 +
            3*i1^2*i2^2*k2*l2*i3 + 3*i2^3*k2*l2*i3 + 5*i1^4*j2*k2*l2*i3 + 4*i1^2*i2*j2*k2*l2*i3 + 5*i2^2*j2*k2*l2*i3 +
            2*i2*j2^2*k2*l2*i3 + 4*j2^3*k2*l2*i3 + 2*i1^4*k2^2*l2*i3 + 5*i1^2*i2*k2^2*l2*i3 + 3*i2^2*k2^2*l2*i3 +
            5*i1^2*j2*k2^2*l2*i3 + 5*i2*j2*k2^2*l2*i3 + 6*j2^2*k2^2*l2*i3 + 4*i1^2*k2^3*l2*i3 + 3*i2*k2^3*l2*i3 +
            j2*k2^3*l2*i3 + 6*k2^4*l2*i3 + i1^6*l2^2*i3 + i1^4*i2*l2^2*i3 + 5*i1^2*i2^2*l2^2*i3 + 6*i2^3*l2^2*i3 +
            3*i1^2*i2*j2*l2^2*i3 + 3*i2^2*j2*l2^2*i3 + 6*i1^2*j2^2*l2^2*i3 + 4*i2*j2^2*l2^2*i3 + j2^3*l2^2*i3 +
            4*i1^4*k2*l2^2*i3 + 3*i1^2*i2*k2*l2^2*i3 + i2^2*k2*l2^2*i3 + i1^2*j2*k2*l2^2*i3 + 5*i2*j2*k2*l2^2*i3 +
            j2^2*k2*l2^2*i3 + i1^2*k2^2*l2^2*i3 + 6*i2*k2^2*l2^2*i3 + 4*j2*k2^2*l2^2*i3 + 5*k2^3*l2^2*i3 + 5*i1^4*l2^3*i3
            + 4*i1^2*i2*l2^3*i3 + i2^2*l2^3*i3 + 6*i1^2*j2*l2^3*i3 + 5*i2*j2*l2^3*i3 + j2^2*l2^3*i3 + 2*i1^2*k2*l2^3*i3 +
            6*i2*k2*l2^3*i3 + j2*k2*l2^3*i3 + 5*k2^2*l2^3*i3 + 3*i1^2*l2^4*i3 + 2*i2*l2^4*i3 + 3*j2*l2^4*i3 + 5*k2*l2^4*i3
            + l2^5*i3 + 2*i1^7*i3^2 + 5*i1^5*i2*i3^2 + 2*i1^5*j2*i3^2 + i1^5*k2*i3^2 + 4*i1*i2*k2^2*i3^2 +
            3*i1*j2*k2^2*i3^2 + 5*i1*k2^3*i3^2 + 2*i1^3*i2*l2*i3^2 + 5*i1^3*j2*l2*i3^2 + 3*i1^3*k2*l2*i3^2 +
            4*i1*i2*k2*l2*i3^2 + 3*i1*j2*k2*l2*i3^2 + 5*i1*k2^2*l2*i3^2 + 5*i1^3*l2^2*i3^2 + 4*i1*i2*l2^2*i3^2 +
            3*i1*j2*l2^2*i3^2 + 6*i1*k2*l2^2*i3^2 + 5*i1*l2^3*i3^2 + 3*i1^8*i2*j3 + i1^6*i2^2*j3 + 3*i1^4*i2^3*j3 +
            3*i1^2*i2^4*j3 + 2*i1^8*j2*j3 + 5*i1^6*i2*j2*j3 + 4*i1^4*i2^2*j2*j3 + 2*i1^2*i2^3*j2*j3 + 5*i1^6*j2^2*j3 +
            4*i1^4*i2*j2^2*j3 + 4*i1^2*i2^2*j2^2*j3 + 3*i1^4*j2^3*j3 + 2*i1^2*i2*j2^3*j3 + 3*i1^2*j2^4*j3 + 2*i1^8*k2*j3 +
            4*i1^4*i2^2*k2*j3 + 3*i2^4*k2*j3 + i1^6*j2*k2*j3 + i1^4*i2*j2*k2*j3 + 4*i1^2*i2^2*j2*k2*j3 + 2*i2^3*j2*k2*j3 +
            5*i1^4*j2^2*k2*j3 + 6*i1^2*i2*j2^2*k2*j3 + 4*i2^2*j2^2*k2*j3 + 4*i1^2*j2^3*k2*j3 + 2*i2*j2^3*k2*j3 +
            3*j2^4*k2*j3 + 5*i1^4*i2*k2^2*j3 + 4*i1^4*j2*k2^2*j3 + i1^2*i2*j2*k2^2*j3 + 4*i2^2*j2*k2^2*j3 +
            4*i1^2*j2^2*k2^2*j3 + 6*i2*j2^2*k2^2*j3 + 4*j2^3*k2^2*j3 + 2*i1^4*k2^3*j3 + 2*i1^2*i2*k2^3*j3 + 4*i2^2*k2^3*j3
            + 3*i2*j2*k2^3*j3 + j2^2*k2^3*j3 + 2*i1^2*k2^4*j3 + i2*k2^4*j3 + 6*j2*k2^4*j3 + 2*k2^5*j3 + 4*i1^8*l2*j3 +
            3*i1^6*i2*l2*j3 + 3*i1^4*i2^2*l2*j3 + 6*i1^2*i2^3*l2*j3 + 2*i2^4*l2*j3 + 3*i1^4*i2*j2*l2*j3 + 6*i2^3*j2*l2*j3
            + 2*i1^4*j2^2*l2*j3 + 3*i1^2*i2*j2^2*l2*j3 + 5*i2^2*j2^2*l2*j3 + 5*i1^2*j2^3*l2*j3 + 6*i2*j2^3*l2*j3 +
            2*j2^4*l2*j3 + 5*i1^6*k2*l2*j3 + 5*i1^4*i2*k2*l2*j3 + i1^2*i2^2*k2*l2*j3 + 2*i2^3*k2*l2*j3 + i2^2*j2*k2*l2*j3
            + 6*i2*j2^2*k2*l2*j3 + 5*j2^3*k2*l2*j3 + 6*i1^4*k2^2*l2*j3 + 2*i1^2*i2*k2^2*l2*j3 + 2*i2^2*k2^2*l2*j3 +
            3*i1^2*j2*k2^2*l2*j3 + 3*i2*j2*k2^2*l2*j3 + 5*j2^2*k2^2*l2*j3 + 4*i1^2*k2^3*l2*j3 + i2*k2^3*l2*j3 +
            4*j2*k2^3*l2*j3 + i1^6*l2^2*j3 + i1^2*i2^2*l2^2*j3 + 6*i2^3*l2^2*j3 + 3*i1^4*j2*l2^2*j3 + 4*i1^2*i2*j2*l2^2*j3
            + i2^2*j2*l2^2*j3 + 4*i1^2*j2^2*l2^2*j3 + i2*j2^2*l2^2*j3 + 6*j2^3*l2^2*j3 + i1^4*k2*l2^2*j3 +
            i1^2*i2*k2*l2^2*j3 + i2^2*k2*l2^2*j3 + 3*i1^2*j2*k2*l2^2*j3 + 2*i2*j2*k2*l2^2*j3 + 3*j2^2*k2*l2^2*j3 +
            4*i1^2*k2^2*l2^2*j3 + 2*i2*k2^2*l2^2*j3 + 6*j2*k2^2*l2^2*j3 + 5*k2^3*l2^2*j3 + 4*i1^4*l2^3*j3 +
            i1^2*i2*l2^3*j3 + 2*i2^2*l2^3*j3 + 5*i1^2*j2*l2^3*j3 + 2*i2*j2*l2^3*j3 + 3*j2^2*l2^3*j3 + 5*i1^2*k2*l2^3*j3 +
            3*i2*k2*l2^3*j3 + 5*j2*k2*l2^3*j3 + k2^2*l2^3*j3 + 2*i1^2*l2^4*j3 + 4*i2*l2^4*j3 + 5*j2*l2^4*j3 + l2^5*j3 +
            3*i1^7*i3*j3 + 6*i1^5*i2*i3*j3 + 2*i1^3*i2^2*i3*j3 + 3*i1^5*j2*i3*j3 + 3*i1^3*i2*j2*i3*j3 + 2*i1^3*j2^2*i3*j3
            + 4*i1^5*k2*i3*j3 + 4*i1^3*i2*k2*i3*j3 + 3*i1^3*j2*k2*i3*j3 + 6*i1^3*k2^2*i3*j3 + i1*i2*k2^2*i3*j3 +
            6*i1*j2*k2^2*i3*j3 + 4*i1*k2^3*i3*j3 + 5*i1^3*i2*l2*i3*j3 + 4*i1*i2^2*l2*i3*j3 + 3*i1^3*j2*l2*i3*j3 +
            6*i1*i2*j2*l2*i3*j3 + 4*i1*j2^2*l2*i3*j3 + 2*i1^3*k2*l2*i3*j3 + 6*i1*i2*k2*l2*i3*j3 + 3*i1*j2*k2*l2*i3*j3 +
            6*i1*k2^2*l2*i3*j3 + 5*i1^3*l2^2*i3*j3 + 5*i1*i2*l2^2*i3*j3 + 4*i1*k2*l2^2*i3*j3 + 6*i1*l2^3*i3*j3 +
            5*i1^4*i3^2*j3 + 3*i1^2*k2*i3^2*j3 + 4*k2^2*i3^2*j3 + i1^2*l2*i3^2*j3 + 6*l2^2*i3^2*j3 + 5*i1^7*j3^2 +
            2*i1^5*i2*j3^2 + 3*i1*i2^3*j3^2 + 4*i1^3*i2*j2*j3^2 + 5*i1*i2^2*j2*j3^2 + 3*i1^3*j2^2*j3^2 + 2*i1*i2*j2^2*j3^2
            + 4*i1*j2^3*j3^2 + 3*i1^5*k2*j3^2 + 3*i1^3*i2*k2*j3^2 + 5*i1^3*j2*k2*j3^2 + 3*i1*i2*j2*k2*j3^2 +
            4*i1*j2^2*k2*j3^2 + i1*i2*k2^2*j3^2 + 5*i1*j2*k2^2*j3^2 + 5*i1^3*i2*l2*j3^2 + 2*i1*i2^2*l2*j3^2 +
            6*i1^3*j2*l2*j3^2 + 5*i1*j2^2*l2*j3^2 + 2*i1^3*k2*l2*j3^2 + i1*i2*k2*l2*j3^2 + 4*i1*j2*k2*l2*j3^2 +
            3*i1*k2^2*l2*j3^2 + 2*i1^3*l2^2*j3^2 + 6*i1*i2*l2^2*j3^2 + 6*i1*j2*l2^2*j3^2 + i1*k2*l2^2*j3^2 + i1*l2^3*j3^2
            + 3*i1^4*i3*j3^2 + 6*i1^2*i2*i3*j3^2 + i1^2*j2*i3*j3^2 + i1^2*k2*i3*j3^2 + 4*i2*k2*i3*j3^2 + 3*j2*k2*i3*j3^2 +
            k2^2*i3*j3^2 + 4*i1^2*l2*i3*j3^2 + i2*l2*i3*j3^2 + 6*j2*l2*i3*j3^2 + 2*k2*l2*i3*j3^2 + 5*l2^2*i3*j3^2 +
            5*i1^4*j3^3 + 3*i1^2*i2*j3^3 + 6*i2^2*j3^3 + i1^2*j2*j3^3 + 2*i2*j2*j3^3 + 6*j2^2*j3^3 + i2*k2*j3^3 +
            4*j2*k2*j3^3 + 4*k2^2*j3^3 + 2*i1^2*l2*j3^3 + i2*l2*j3^3 + 4*j2*l2*j3^3 + 3*k2*l2*j3^3 + 6*l2^2*j3^3 +
            6*i1*i3*j3^3 + 5*i1*j3^4 + 2*i1^10*k3 + 5*i1^8*i2*k3 + 3*i1^6*i2^2*k3 + 5*i1^2*i2^4*k3 + 5*i1^8*j2*k3 +
            6*i1^4*i2^2*j2*k3 + i1^2*i2^3*j2*k3 + i1^6*j2^2*k3 + 2*i1^4*i2*j2^2*k3 + 2*i1^2*i2^2*j2^2*k3 + 6*i1^4*j2^3*k3
            + i1^2*i2*j2^3*k3 + 5*i1^2*j2^4*k3 + 6*i1^8*k2*k3 + 2*i1^6*i2*k2*k3 + 5*i1^4*i2^2*k2*k3 + 2*i1^2*i2^3*k2*k3 +
            3*i2^4*k2*k3 + i1^4*i2*j2*k2*k3 + 4*i1^2*i2^2*j2*k2*k3 + 2*i2^3*j2*k2*k3 + 3*i1^4*j2^2*k2*k3 +
            4*i2^2*j2^2*k2*k3 + i1^2*j2^3*k2*k3 + 2*i2*j2^3*k2*k3 + 3*j2^4*k2*k3 + 5*i1^6*k2^2*k3 + 4*i1^4*i2*k2^2*k3 +
            5*i1^2*i2^2*k2^2*k3 + 4*i2^3*k2^2*k3 + 3*i1^4*j2*k2^2*k3 + i1^2*i2*j2*k2^2*k3 + 6*i2^2*j2*k2^2*k3 +
            2*i1^2*j2^2*k2^2*k3 + 4*i2*j2^2*k2^2*k3 + 6*i1^4*k2^3*k3 + 6*i1^2*i2*k2^3*k3 + 6*i2^2*k2^3*k3 +
            2*i1^2*j2*k2^3*k3 + 4*i2*j2*k2^3*k3 + 6*j2^2*k2^3*k3 + 3*i2*k2^4*k3 + 4*j2*k2^4*k3 + 4*k2^5*k3 + 5*i1^8*l2*k3
            + 5*i1^6*i2*l2*k3 + 5*i1^4*i2^2*l2*k3 + 5*i2^4*l2*k3 + i1^6*j2*l2*k3 + i1^4*i2*j2*l2*k3 + 2*i1^2*i2^2*j2*l2*k3
            + i2^3*j2*l2*k3 + 6*i1^4*j2^2*l2*k3 + 3*i1^2*i2*j2^2*l2*k3 + 2*i2^2*j2^2*l2*k3 + 2*i1^2*j2^3*l2*k3 +
            i2*j2^3*l2*k3 + 5*j2^4*l2*k3 + 3*i1^6*k2*l2*k3 + i1^4*i2*k2*l2*k3 + 2*i2^3*k2*l2*k3 + 4*i1^4*j2*k2*l2*k3 +
            2*i2^2*j2*k2*l2*k3 + 3*i1^2*j2^2*k2*l2*k3 + 4*i2*j2^2*k2*l2*k3 + 6*j2^3*k2*l2*k3 + 2*i1^4*k2^2*l2*k3 +
            i1^2*i2*k2^2*l2*k3 + 2*j2^2*k2^2*l2*k3 + 3*i2*k2^3*l2*k3 + 6*j2*k2^3*l2*k3 + 3*k2^4*l2*k3 + 6*i1^6*l2^2*k3 +
            5*i1^4*i2*l2^2*k3 + 4*i1^2*i2^2*l2^2*k3 + 2*i2^3*l2^2*k3 + 3*i1^4*j2*l2^2*k3 + 4*i1^2*i2*j2*l2^2*k3 +
            i2^2*j2*l2^2*k3 + 6*i1^2*j2^2*l2^2*k3 + 6*i2*j2^2*l2^2*k3 + 5*j2^3*l2^2*k3 + 6*i1^4*k2*l2^2*k3 +
            3*i1^2*i2*k2*l2^2*k3 + 6*i1^2*j2*k2*l2^2*k3 + 6*i2*j2*k2*l2^2*k3 + 2*j2^2*k2*l2^2*k3 + 6*i1^2*k2^2*l2^2*k3 +
            2*i2*k2^2*l2^2*k3 + 6*k2^3*l2^2*k3 + 4*i1^4*l2^3*k3 + 2*i1^2*i2*l2^3*k3 + i2^2*l2^3*k3 + 6*i1^2*j2*l2^3*k3 +
            5*i2*j2*l2^3*k3 + j2^2*l2^3*k3 + 3*i1^2*k2*l2^3*k3 + 4*i2*k2*l2^3*k3 + 3*j2*k2*l2^3*k3 + 4*k2^2*l2^3*k3 +
            5*i1^2*l2^4*k3 + 5*i2*l2^4*k3 + 4*j2*l2^4*k3 + 3*k2*l2^4*k3 + 3*l2^5*k3 + 4*i1^7*j3*k3 + i1^5*i2*j3*k3 +
            3*i1^3*i2^2*j3*k3 + 2*i1*i2^3*j3*k3 + 6*i1^5*j2*j3*k3 + 5*i1^3*i2*j2*j3*k3 + i1*i2^2*j2*j3*k3 +
            6*i1^3*j2^2*j3*k3 + 6*i1*i2*j2^2*j3*k3 + 5*i1*j2^3*j3*k3 + 3*i1^5*k2*j3*k3 + 5*i1*i2^2*k2*j3*k3 +
            2*i1^3*j2*k2*j3*k3 + 6*i1*i2*j2*k2*j3*k3 + 3*i1*j2^2*k2*j3*k3 + 5*i1^3*k2^2*j3*k3 + 3*i1*j2*k2^2*j3*k3 +
            6*i1^5*l2*j3*k3 + 3*i1^3*i2*l2*j3*k3 + i1*i2^2*l2*j3*k3 + 5*i1^3*j2*l2*j3*k3 + 2*i1*i2*j2*l2*j3*k3 +
            4*i1*j2^2*l2*j3*k3 + 2*i1^3*k2*l2*j3*k3 + 3*i1*i2*k2*l2*j3*k3 + 5*i1*j2*k2*l2*j3*k3 + 2*i1*k2^2*l2*j3*k3 +
            6*i1^3*l2^2*j3*k3 + 5*i1*i2*l2^2*j3*k3 + 2*i1*l2^3*j3*k3 + 4*i1^4*j3^2*k3 + 5*i1^2*i2*j3^2*k3 + 4*i2^2*j3^2*k3
            + 2*i1^2*j2*j3^2*k3 + 6*i2*j2*j3^2*k3 + 4*j2^2*j3^2*k3 + 2*i1^2*k2*j3^2*k3 + 2*i2*k2*j3^2*k3 + 3*j2*k2*j3^2*k3
            + 3*k2^2*j3^2*k3 + 4*i1^2*l2*j3^2*k3 + 3*i2*l2*j3^2*k3 + 2*j2*l2*j3^2*k3 + 2*k2*l2*j3^2*k3 + 2*l2^2*j3^2*k3 +
            6*i1*j3^3*k3 + 3*i1^7*k3^2 + 4*i1^5*i2*k3^2 + 4*i1^3*i2^2*k3^2 + 3*i1*i2^3*k3^2 + 5*i1^5*j2*k3^2 +
            4*i1^3*i2*j2*k3^2 + 5*i1*i2^2*j2*k3^2 + 6*i1^3*j2^2*k3^2 + 2*i1*i2*j2^2*k3^2 + 4*i1*j2^3*k3^2 +
            i1^3*i2*k2*k3^2 + 5*i1*i2^2*k2*k3^2 + 5*i1^3*j2*k2*k3^2 + 2*i1*j2^2*k2*k3^2 + 6*i1^3*k2^2*k3^2 +
            i1*i2*k2^2*k3^2 + 5*i1*j2*k2^2*k3^2 + 5*i1*k2^3*k3^2 + 3*i1*i2^2*l2*k3^2 + 6*i1^3*j2*l2*k3^2 +
            5*i1*i2*j2*l2*k3^2 + 6*i1*j2^2*l2*k3^2 + 2*i1^3*k2*l2*k3^2 + i1*j2*k2*l2*k3^2 + 5*i1*k2^2*l2*k3^2 +
            6*i1^3*l2^2*k3^2 + i1*i2*l2^2*k3^2 + 2*i1*j2*l2^2*k3^2 + i1*k2*l2^2*k3^2 + 4*i1*l2^3*k3^2 + i1^4*j3*k3^2 +
            6*i1^2*i2*j3*k3^2 + 3*i2^2*j3*k3^2 + 6*i1^2*j2*j3*k3^2 + i2*j2*j3*k3^2 + 3*j2^2*j3*k3^2 + 5*i1^2*k2*j3*k3^2 +
            2*i2*k2*j3*k3^2 + j2*k2*j3*k3^2 + i1^2*l2*j3*k3^2 + 4*i2*l2*j3*k3^2 + 4*j2*l2*j3*k3^2 + 4*l2^2*j3*k3^2 +
            4*i1*j3^2*k3^2 + i1^4*k3^3 + 5*i1^2*j2*k3^3 + 2*i1^2*k2*k3^3 + 5*i2*k2*k3^3 + 4*j2*k2*k3^3 + 4*k2^2*k3^3 +
            5*i1^2*l2*k3^3 + i2*l2*k3^3 + j2*l2*k3^3 + 5*k2*l2*k3^3 + i1*j3*k3^3 + 6*i1*k3^4 + 5*i1^10*l3 + 3*i1^8*i2*l3 +
            6*i1^6*i2^2*l3 + 4*i1^4*i2^3*l3 + 4*i1^2*i2^4*l3 + 2*i2^5*l3 + 5*i1^6*i2*j2*l3 + i1^2*i2^3*j2*l3 +
            4*i2^4*j2*l3 + 5*i1^6*j2^2*l3 + 5*i1^4*i2*j2^2*l3 + i1^2*i2^2*j2^2*l3 + 6*i2^3*j2^2*l3 + 5*i1^4*j2^3*l3 +
            i2^2*j2^3*l3 + i1^2*j2^4*l3 + 3*i2*j2^4*l3 + 5*j2^5*l3 + 4*i1^8*k2*l3 + 3*i1^6*i2*k2*l3 + 3*i1^4*i2^2*k2*l3 +
            4*i1^6*j2*k2*l3 + 3*i1^4*i2*j2*k2*l3 + 2*i1^2*i2^2*j2*k2*l3 + 4*i2^3*j2*k2*l3 + 2*i1^2*i2*j2^2*k2*l3 +
            2*i2^2*j2^2*k2*l3 + 3*i1^2*j2^3*k2*l3 + 5*i2*j2^3*k2*l3 + 3*j2^4*k2*l3 + 3*i1^4*i2*k2^2*l3 +
            2*i1^2*i2^2*k2^2*l3 + 2*i2^3*k2^2*l3 + 4*i1^4*j2*k2^2*l3 + 5*i1^2*i2*j2*k2^2*l3 + 2*i1^2*j2^2*k2^2*l3 +
            i2*j2^2*k2^2*l3 + 4*j2^3*k2^2*l3 + 6*i1^4*k2^3*l3 + 5*i1^2*i2*k2^3*l3 + i2^2*k2^3*l3 + 6*i1^2*j2*k2^3*l3 +
            j2^2*k2^3*l3 + 6*i1^2*k2^4*l3 + 6*i2*k2^4*l3 + 5*j2*k2^4*l3 + k2^5*l3 + 6*i1^8*l2*l3 + 3*i1^6*i2*l2*l3 +
            4*i1^4*i2^2*l2*l3 + 6*i1^2*i2^3*l2*l3 + 3*i2^4*l2*l3 + i1^6*j2*l2*l3 + 4*i1^4*i2*j2*l2*l3 + i1^2*i2^2*j2*l2*l3
            + i2^3*j2*l2*l3 + 2*i1^4*j2^2*l2*l3 + i1^2*i2*j2^2*l2*l3 + 6*i1^2*j2^3*l2*l3 + 6*i2*j2^3*l2*l3 + 4*j2^4*l2*l3
            + 6*i1^4*i2*k2*l2*l3 + 2*i1^2*i2^2*k2*l2*l3 + 2*i2^3*k2*l2*l3 + i1^4*j2*k2*l2*l3 + 3*i1^2*i2*j2*k2*l2*l3 +
            6*i2^2*j2*k2*l2*l3 + 5*i1^2*j2^2*k2*l2*l3 + 2*i2*j2^2*k2*l2*l3 + 4*j2^3*k2*l2*l3 + i1^4*k2^2*l2*l3 +
            3*i1^2*i2*k2^2*l2*l3 + 4*i2^2*k2^2*l2*l3 + 5*i1^2*j2*k2^2*l2*l3 + 6*i2*j2*k2^2*l2*l3 + 3*j2^2*k2^2*l2*l3 +
            5*i1^2*k2^3*l2*l3 + 2*i2*k2^3*l2*l3 + j2*k2^3*l2*l3 + 5*k2^4*l2*l3 + 5*i1^6*l2^2*l3 + 4*i1^4*i2*l2^2*l3 +
            i1^4*j2*l2^2*l3 + i1^2*i2*j2*l2^2*l3 + 6*i2^2*j2*l2^2*l3 + 2*i1^2*j2^2*l2^2*l3 + 4*i2*j2^2*l2^2*l3 +
            4*j2^3*l2^2*l3 + 3*i1^4*k2*l2^2*l3 + 4*i1^2*i2*k2*l2^2*l3 + 6*i1^2*j2*k2*l2^2*l3 + 2*i2*j2*k2*l2^2*l3 +
            6*j2^2*k2*l2^2*l3 + 3*i1^2*k2^2*l2^2*l3 + 3*i2*k2^2*l2^2*l3 + 6*j2*k2^2*l2^2*l3 + 3*k2^3*l2^2*l3 +
            4*i1^4*l2^3*l3 + 5*i1^2*i2*l2^3*l3 + 2*i2^2*l2^3*l3 + 5*i2*j2*l2^3*l3 + j2^2*l2^3*l3 + 4*i1^2*k2*l2^3*l3 +
            4*i2*k2*l2^3*l3 + j2*k2*l2^3*l3 + 3*k2^2*l2^3*l3 + 2*i1^2*l2^4*l3 + 4*i2*l2^4*l3 + j2*l2^4*l3 + 2*k2*l2^4*l3 +
            5*l2^5*l3 + 4*i1^7*i3*l3 + 2*i1^5*i2*i3*l3 + 3*i1^3*i2^2*i3*l3 + 3*i1*i2^3*i3*l3 + 3*i1^5*j2*i3*l3 +
            5*i1*i2^2*j2*i3*l3 + 4*i1^3*j2^2*i3*l3 + 2*i1*i2*j2^2*i3*l3 + 4*i1*j2^3*i3*l3 + 3*i1^5*k2*i3*l3 +
            6*i1^3*i2*k2*i3*l3 + i1*i2^2*k2*i3*l3 + 5*i1^3*j2*k2*i3*l3 + 3*i1*i2*j2*k2*i3*l3 + 3*i1*j2^2*k2*i3*l3 +
            4*i1^3*k2^2*i3*l3 + 3*i1*i2*k2^2*i3*l3 + 6*i1*j2*k2^2*i3*l3 + 3*i1*k2^3*i3*l3 + 2*i1^5*l2*i3*l3 +
            4*i1^3*i2*l2*i3*l3 + 3*i1*i2^2*l2*i3*l3 + 5*i1^3*j2*l2*i3*l3 + 4*i1*i2*j2*l2*i3*l3 + 5*i1^3*k2*l2*i3*l3 +
            3*i1*j2*k2*l2*i3*l3 + 6*i1*k2^2*l2*i3*l3 + i1^3*l2^2*i3*l3 + 6*i1*i2*l2^2*i3*l3 + 6*i1*j2*l2^2*i3*l3 +
            6*i1*k2*l2^2*i3*l3 + 4*i1*l2^3*i3*l3 + 3*i1^4*i3^2*l3 + 5*i1^2*i2*i3^2*l3 + 2*i1^2*j2*i3^2*l3 +
            6*i1^2*k2*i3^2*l3 + 6*i2*k2*i3^2*l3 + j2*k2*i3^2*l3 + 3*k2^2*i3^2*l3 + 4*i1^2*l2*i3^2*l3 + 5*i2*l2*i3^2*l3 +
            2*j2*l2*i3^2*l3 + k2*l2*i3^2*l3 + l2^2*i3^2*l3 + 5*i1^7*j3*l3 + 4*i1^5*i2*j3*l3 + 2*i1*i2^3*j3*l3 +
            6*i1^5*j2*j3*l3 + i1^3*i2*j2*j3*l3 + 5*i1*i2^2*j2*j3*l3 + 6*i1^3*j2^2*j3*l3 + 5*i1*i2*j2^2*j3*l3 +
            2*i1*j2^3*j3*l3 + i1^5*k2*j3*l3 + 4*i1^3*i2*k2*j3*l3 + 3*i1*i2^2*k2*j3*l3 + 5*i1^3*j2*k2*j3*l3 +
            6*i1*i2*j2*k2*j3*l3 + 3*i1*j2^2*k2*j3*l3 + 5*i1^3*k2^2*j3*l3 + 2*i1*i2*k2^2*j3*l3 + 4*i1*k2^3*j3*l3 +
            3*i1^5*l2*j3*l3 + 3*i1*i2^2*l2*j3*l3 + i1^3*j2*l2*j3*l3 + 3*i1*i2*j2*l2*j3*l3 + 3*i1*j2^2*l2*j3*l3 +
            4*i1^3*k2*l2*j3*l3 + 6*i1*j2*k2*l2*j3*l3 + 5*i1*k2^2*l2*j3*l3 + 2*i1^3*l2^2*j3*l3 + i1*i2*l2^2*j3*l3 +
            2*i1*j2*l2^2*j3*l3 + 4*i1*k2*l2^2*j3*l3 + 5*i1*l2^3*j3*l3 + 5*i1^4*j3^2*l3 + 4*i1^2*i2*j3^2*l3 +
            4*i2^2*j3^2*l3 + i1^2*j2*j3^2*l3 + 2*i2*j2*j3^2*l3 + j2^2*j3^2*l3 + 4*i1^2*k2*j3^2*l3 + 5*i2*k2*j3^2*l3 +
            3*j2*k2*j3^2*l3 + 4*k2^2*j3^2*l3 + 2*i1^2*l2*j3^2*l3 + i2*l2*j3^2*l3 + j2*l2*j3^2*l3 + 5*l2^2*j3^2*l3 +
            2*i1*j3^3*l3 + 4*i1^7*l3^2 + 3*i1^5*i2*l3^2 + 6*i1^3*i2^2*l3^2 + 2*i1*i2^3*l3^2 + i1^5*j2*l3^2 +
            3*i1*i2^2*j2*l3^2 + 2*i1^3*j2^2*l3^2 + 3*i1*i2*j2^2*l3^2 + 6*i1*j2^3*l3^2 + 4*i1^5*k2*l3^2 + 3*i1^3*j2*k2*l3^2
            + 4*i1*i2*j2*k2*l3^2 + 2*i1*j2^2*k2*l3^2 + i1^3*k2^2*l3^2 + 5*i1*i2*k2^2*l3^2 + i1*j2*k2^2*l3^2 + i1*k2^3*l3^2
            + 4*i1^5*l2*l3^2 + 5*i1^3*i2*l2*l3^2 + 4*i1*i2^2*l2*l3^2 + 6*i1^3*j2*l2*l3^2 + 2*i1*i2*j2*l2*l3^2 +
            i1*j2^2*l2*l3^2 + i1*j2*k2*l2*l3^2 + 2*i1*k2^2*l2*l3^2 + 6*i1^3*l2^2*l3^2 + 2*i1*i2*l2^2*l3^2 +
            5*i1*j2*l2^2*l3^2 + i1*k2*l2^2*l3^2 + 6*i1*l2^3*l3^2 + i1^4*i3*l3^2 + 5*i2^2*i3*l3^2 + i1^2*j2*i3*l3^2 +
            2*j2^2*i3*l3^2 + 6*i1^2*k2*i3*l3^2 + 3*i2*k2*i3*l3^2 + 5*k2^2*i3*l3^2 + 5*i1^2*l2*i3*l3^2 + 2*j2*l2*i3*l3^2 +
            2*k2*l2*i3*l3^2 + 2*l2^2*i3*l3^2 + 6*i1*i3^2*l3^2 + i1^2*i2*j3*l3^2 + 5*i1^2*j2*j3*l3^2 + 3*i1^2*k2*j3*l3^2 +
            6*i2*k2*j3*l3^2 + 5*j2*k2*j3*l3^2 + k2^2*j3*l3^2 + i1^2*l2*j3*l3^2 + 3*i2*l2*j3*l3^2 + 4*j2*l2*j3*l3^2 +
            2*k2*l2*j3*l3^2 + l2^2*j3*l3^2 + 2*i1*j3^2*l3^2 + 2*i1^4*l3^3 + 6*i1^2*i2*l3^3 + 6*i2^2*l3^3 + 2*i1^2*j2*l3^3
            + j2^2*l3^3 + 5*i1^2*k2*l3^3 + i2*k2*l3^3 + 6*j2*k2*l3^3 + 3*k2^2*l3^3 + 2*i1^2*l2*l3^3 + 6*i2*l2*l3^3 +
            5*l2^2*l3^3 + 4*i1*i3*l3^3 + 2*i1*j3*l3^3 + 2*i1*l3^4 + i1^9*i4 + i1^7*i2*i4 + 3*i1^5*i2^2*i4 + i1^3*i2^3*i4 +
            5*i1*i2^4*i4 + i1^7*j2*i4 + 2*i1^5*i2*j2*i4 + 4*i1^3*i2^2*j2*i4 + i1*i2^3*j2*i4 + 3*i1^5*j2^2*i4 +
            3*i1^3*i2*j2^2*i4 + 2*i1*i2^2*j2^2*i4 + 6*i1^3*j2^3*i4 + i1*i2*j2^3*i4 + 5*i1*j2^4*i4 + 2*i1^5*i2*k2*i4 +
            5*i1^3*i2^2*k2*i4 + 4*i1*i2^3*k2*i4 + 2*i1^5*j2*k2*i4 + i1^3*i2*j2*k2*i4 + 3*i1^3*j2^2*k2*i4 +
            2*i1*i2*j2^2*k2*i4 + i1*j2^3*k2*i4 + 4*i1^3*i2*k2^2*i4 + 3*i1*i2^2*k2^2*i4 + 4*i1^3*j2*k2^2*i4 +
            6*i1*i2*j2*k2^2*i4 + 4*i1*j2^2*k2^2*i4 + 5*i1*i2*k2^3*i4 + 5*i1*j2*k2^3*i4 + 6*i1*k2^4*i4 + 5*i1^7*l2*i4 +
            i1^5*i2*l2*i4 + 6*i1^3*i2^2*l2*i4 + 6*i1*i2^3*l2*i4 + 5*i1^5*j2*l2*i4 + 2*i1^3*i2*j2*l2*i4 +
            5*i1*i2^2*j2*l2*i4 + 2*i1^3*j2^2*l2*i4 + 3*i1*j2^3*l2*i4 + 2*i1^5*k2*l2*i4 + 4*i1^3*i2*k2*l2*i4 +
            i1*i2^2*k2*l2*i4 + 3*i1^3*j2*k2*l2*i4 + i1*i2*j2*k2*l2*i4 + 3*i1^3*k2^2*l2*i4 + 5*i1*i2*k2^2*l2*i4 +
            5*i1*j2*k2^2*l2*i4 + 6*i1*k2^3*l2*i4 + 5*i1^5*l2^2*i4 + 6*i1^3*i2*l2^2*i4 + 3*i1*i2^2*l2^2*i4 +
            6*i1^3*j2*l2^2*i4 + 6*i1*i2*j2*l2^2*i4 + 6*i1*j2^2*l2^2*i4 + 6*i1^3*k2*l2^2*i4 + i1*i2*k2*l2^2*i4 +
            4*i1*j2*k2*l2^2*i4 + 5*i1*k2^2*l2^2*i4 + 4*i1^3*l2^3*i4 + 2*i1*i2*l2^3*i4 + 2*i1*j2*l2^3*i4 + 2*i1*k2*l2^3*i4
            + 6*i1*l2^4*i4 + 4*i1^4*i2*i3*i4 + 6*i1^2*i2^2*i3*i4 + 5*i1^4*j2*i3*i4 + 2*i1^2*i2*j2*i3*i4 +
            6*i1^2*j2^2*i3*i4 + 6*i1^4*k2*i3*i4 + 4*i1^2*i2*k2*i3*i4 + 6*i2^2*k2*i3*i4 + 2*i1^2*j2*k2*i3*i4 +
            2*i2*j2*k2*i3*i4 + 6*j2^2*k2*i3*i4 + i1^2*k2^2*i3*i4 + 5*i2*k2^2*i3*i4 + 4*j2*k2^2*i3*i4 + 2*k2^3*i3*i4 +
            4*i1^4*l2*i3*i4 + 6*i1^2*i2*l2*i3*i4 + 6*i2^2*l2*i3*i4 + 6*i1^2*j2*l2*i3*i4 + 2*i2*j2*l2*i3*i4 +
            6*j2^2*l2*i3*i4 + 4*i1^2*k2*l2*i3*i4 + 4*i2*k2*l2*i3*i4 + 2*k2^2*l2*i3*i4 + i2*l2^2*i3*i4 + 5*j2*l2^2*i3*i4 +
            4*k2*l2^2*i3*i4 + 6*l2^3*i3*i4 + 5*i1^3*i3^2*i4 + 5*i1*k2*i3^2*i4 + 2*i1*l2*i3^2*i4 + 5*i1^4*i2*j3*i4 +
            4*i2^3*j3*i4 + 6*i1^4*j2*j3*i4 + i1^2*i2*j2*j3*i4 + 2*i2^2*j2*j3*i4 + 6*i1^2*j2^2*j3*i4 + 5*i2*j2^2*j3*i4 +
            3*j2^3*j3*i4 + 2*i1^4*k2*j3*i4 + 5*i1^2*i2*k2*j3*i4 + 2*i2^2*k2*j3*i4 + 3*i1^2*j2*k2*j3*i4 + 4*i2*j2*k2*j3*i4
            + j2^2*k2*j3*i4 + 4*i1^2*k2^2*j3*i4 + 3*j2*k2^2*j3*i4 + 6*k2^3*j3*i4 + 3*i1^4*l2*j3*i4 + 5*i2^2*l2*j3*i4 +
            5*i1^2*j2*l2*j3*i4 + 2*i2*j2*l2*j3*i4 + 5*i1^2*k2*l2*j3*i4 + 4*i2*k2*l2*j3*i4 + 6*j2*k2*l2*j3*i4 +
            5*k2^2*l2*j3*i4 + 2*i2*l2^2*j3*i4 + 2*j2*l2^2*j3*i4 + 6*k2*l2^2*j3*i4 + 2*l2^3*j3*i4 + 2*i1^3*i3*j3*i4 +
            5*i1*i2*i3*j3*i4 + 2*i1*j2*i3*j3*i4 + 4*i1*k2*i3*j3*i4 + 2*i1*l2*i3*j3*i4 + 5*i1^3*j3^2*i4 + 6*i1*i2*j3^2*i4 +
            4*i1*j2*j3^2*i4 + 3*i1*k2*j3^2*i4 + 2*i1*l2*j3^2*i4 + 3*i3*j3^2*i4 + 2*j3^3*i4 + 5*i1^6*k3*i4 + i2^3*k3*i4 +
            2*i1^4*j2*k3*i4 + i1^2*i2*j2*k3*i4 + 4*i2^2*j2*k3*i4 + 6*i1^2*j2^2*k3*i4 + 3*i2*j2^2*k3*i4 + 6*j2^3*k3*i4 +
            i1^4*k2*k3*i4 + 4*i1^2*i2*k2*k3*i4 + 2*i2^2*k2*k3*i4 + 3*i1^2*j2*k2*k3*i4 + 2*i2*j2*k2*k3*i4 + 3*j2^2*k2*k3*i4
            + 2*i1^2*k2^2*k3*i4 + 2*i2*k2^2*k3*i4 + 6*j2*k2^2*k3*i4 + k2^3*k3*i4 + 2*i1^4*l2*k3*i4 + 6*i2^2*l2*k3*i4 +
            3*i1^2*j2*l2*k3*i4 + j2^2*l2*k3*i4 + 3*i1^2*k2*l2*k3*i4 + 2*i2*k2*l2*k3*i4 + 4*j2*k2*l2*k3*i4 +
            2*k2^2*l2*k3*i4 + 4*i1^2*l2^2*k3*i4 + 2*i2*l2^2*k3*i4 + 3*j2*l2^2*k3*i4 + 3*k2*l2^2*k3*i4 + 6*i1^3*j3*k3*i4 +
            i1*i2*j3*k3*i4 + 4*i1*j2*j3*k3*i4 + 5*i1*k2*j3*k3*i4 + 5*i1*l2*j3*k3*i4 + 4*i1^9*j4 + 5*i1^5*i2^2*j4 +
            2*i1^3*i2^3*j4 + 3*i1^7*j2*j4 + i1^5*i2*j2*j4 + i1^3*i2^2*j2*j4 + i1^5*j2^2*j4 + 6*i1^3*i2*j2^2*j4 +
            5*i1^3*j2^3*j4 + 5*i1^7*k2*j4 + 2*i1^3*i2^2*k2*j4 + 5*i1*i2^3*k2*j4 + 5*i1^5*j2*k2*j4 + 4*i1^3*i2*j2*k2*j4 +
            6*i1*i2^2*j2*k2*j4 + i1^3*j2^2*k2*j4 + i1*i2*j2^2*k2*j4 + 2*i1*j2^3*k2*j4 + 6*i1^5*k2^2*j4 + 3*i1^3*i2*k2^2*j4
            + 2*i1*i2^2*k2^2*j4 + 5*i1^3*j2*k2^2*j4 + 5*i1*j2^2*k2^2*j4 + 4*i1*i2*k2^3*j4 + 2*i1*j2*k2^3*j4 + 3*i1^7*l2*j4
            + i1^5*i2*l2*j4 + 4*i1^3*i2^2*l2*j4 + i1*i2^3*l2*j4 + 2*i1^5*j2*l2*j4 + i1^3*i2*j2*l2*j4 + 4*i1*i2^2*j2*l2*j4
            + 2*i1^3*j2^2*l2*j4 + 3*i1*i2*j2^2*l2*j4 + 6*i1*j2^3*l2*j4 + 4*i1^5*k2*l2*j4 + 5*i1^3*i2*k2*l2*j4 +
            4*i1*i2^2*k2*l2*j4 + i1^3*j2*k2*l2*j4 + i1*i2*j2*k2*l2*j4 + 2*i1*j2^2*k2*l2*j4 + 5*i1^3*k2^2*l2*j4 +
            4*i1*j2*k2^2*l2*j4 + 5*i1*k2^3*l2*j4 + 4*i1^5*l2^2*j4 + 3*i1^3*i2*l2^2*j4 + 2*i1*i2^2*l2^2*j4 +
            2*i1^3*j2*l2^2*j4 + 5*i1*j2^2*l2^2*j4 + 2*i1^3*k2*l2^2*j4 + i1*i2*k2*l2^2*j4 + 3*i1*j2*k2*l2^2*j4 +
            4*i1*k2^2*l2^2*j4 + i1^3*l2^3*j4 + 5*i1*i2*l2^3*j4 + 6*i1*j2*l2^3*j4 + 6*i1*k2*l2^3*j4 + 6*i1^6*j3*j4 +
            5*i1^4*i2*j3*j4 + 4*i1^2*i2^2*j3*j4 + 6*i1^4*j2*j3*j4 + 6*i1^2*i2*j2*j3*j4 + 4*i1^2*j2^2*j3*j4 +
            4*i1^4*k2*j3*j4 + 4*i2^2*k2*j3*j4 + 3*i1^2*j2*k2*j3*j4 + 6*i2*j2*k2*j3*j4 + 4*j2^2*k2*j3*j4 +
            2*i1^2*k2^2*j3*j4 + 4*i2*k2^2*j3*j4 + 5*k2^3*j3*j4 + 4*i1^4*l2*j3*j4 + 5*i1^2*i2*l2*j3*j4 + 3*i2^2*l2*j3*j4 +
            i1^2*j2*l2*j3*j4 + i2*j2*l2*j3*j4 + 3*j2^2*l2*j3*j4 + 6*i1^2*k2*l2*j3*j4 + 3*i2*k2*l2*j3*j4 + 3*j2*k2*l2*j3*j4
            + 4*k2^2*l2*j3*j4 + 6*i1^2*l2^2*j3*j4 + 2*i2*l2^2*j3*j4 + j2*l2^2*j3*j4 + 2*k2*l2^2*j3*j4 + 6*i1^3*j3^2*j4 +
            5*i1*i2*j3^2*j4 + 2*i1*j2*j3^2*j4 + 6*i1*k2*j3^2*j4 + 6*i1*l2*j3^2*j4 + 5*j3^3*j4 + 3*i1^9*k4 + 3*i1^7*i2*k4 +
            5*i1^5*i2^2*k4 + 4*i1^5*i2*j2*k4 + 5*i1^5*j2^2*k4 + i1^7*k2*k4 + 2*i1^5*i2*k2*k4 + 3*i1^3*i2^2*k2*k4 +
            i1^3*i2*j2*k2*k4 + 3*i1^3*j2^2*k2*k4 + 3*i1^3*i2*k2^2*k4 + 6*i1*i2^2*k2^2*k4 + 3*i1^3*j2*k2^2*k4 +
            2*i1*i2*j2*k2^2*k4 + 6*i1*j2^2*k2^2*k4 + 2*i1*i2*k2^3*k4 + 2*i1*j2*k2^3*k4 + i1*k2^4*k4 + 3*i1^7*l2*k4 +
            6*i1^5*i2*l2*k4 + i1^5*j2*l2*k4 + 4*i1^5*k2*l2*k4 + 3*i1^3*i2*k2*l2*k4 + 5*i1*i2^2*k2*l2*k4 +
            3*i1^3*j2*k2*l2*k4 + 4*i1*i2*j2*k2*l2*k4 + 5*i1*j2^2*k2*l2*k4 + 3*i1^3*k2^2*l2*k4 + 5*i1*i2*k2^2*l2*k4 +
            i1*j2*k2^2*l2*k4 + 5*i1*k2^3*l2*k4 + i1^5*l2^2*k4 + 6*i1^3*i2*l2^2*k4 + 4*i1*i2^2*l2^2*k4 + 6*i1*i2*j2*l2^2*k4
            + 4*i1*j2^2*l2^2*k4 + 2*i1^3*k2*l2^2*k4 + 6*i1*i2*k2*l2^2*k4 + 2*i1*j2*k2*l2^2*k4 + 4*i1^3*l2^3*k4 +
            3*i1*i2*l2^3*k4 + 3*i1*j2*l2^3*k4 + i1^6*i3*k4 + i1^4*k2*i3*k4 + i1^2*k2^2*i3*k4 + 4*k2^3*i3*k4 +
            2*i1^4*l2*i3*k4 + 5*i1^2*k2*l2*i3*k4 + 2*k2^2*l2*i3*k4 + 5*i1^2*l2^2*i3*k4 + 2*k2*l2^2*i3*k4 + 2*l2^3*i3*k4 +
            5*i1^4*i2*k3*k4 + 2*i1^4*j2*k3*k4 + 4*i1^4*k2*k3*k4 + 3*i1^2*i2*k2*k3*k4 + 4*i1^2*j2*k2*k3*k4 +
            4*i1^2*k2^2*k3*k4 + i2*k2^2*k3*k4 + 6*j2*k2^2*k3*k4 + 6*k2^3*k3*k4 + i1^4*l2*k3*k4 + 5*i1^2*i2*l2*k3*k4 +
            2*i1^2*j2*l2*k3*k4 + 4*i1^2*k2*l2*k3*k4 + 4*i2*k2*l2*k3*k4 + 3*j2*k2*l2*k3*k4 + 4*i1^2*l2^2*k3*k4 +
            5*i2*l2^2*k3*k4 + 2*j2*l2^2*k3*k4 + 6*k2*l2^2*k3*k4 + 5*l2^3*k3*k4 + 5*i1^3*k3^2*k4 + 4*i1*k2*k3^2*k4 +
            6*i1*l2*k3^2*k4 + 3*i1^9*l4 + 2*i1^7*i2*l4 + 6*i1^5*i2^2*l4 + 6*i1^3*i2^3*l4 + 4*i1*i2^4*l4 + 3*i1^7*j2*l4 +
            3*i1^5*i2*j2*l4 + 5*i1*i2^3*j2*l4 + 5*i1^5*j2^2*l4 + 3*i1^3*i2*j2^2*l4 + 3*i1*i2^2*j2^2*l4 + 5*i1^3*j2^3*l4 +
            5*i1*i2*j2^3*l4 + 4*i1*j2^4*l4 + 5*i1^7*k2*l4 + 3*i1^3*i2^2*k2*l4 + 3*i1*i2^3*k2*l4 + 5*i1^5*j2*k2*l4 +
            6*i1^3*i2*j2*k2*l4 + 4*i1*i2^2*j2*k2*l4 + 3*i1^3*j2^2*k2*l4 + 4*i1*i2*j2^2*k2*l4 + 3*i1*j2^3*k2*l4 +
            6*i1^5*k2^2*l4 + 4*i1^3*i2*k2^2*l4 + 5*i1*i2^2*k2^2*l4 + 6*i1*j2^2*k2^2*l4 + 2*i1^3*k2^3*l4 + 4*i1*i2*k2^3*l4
            + 6*i1*j2*k2^3*l4 + 4*i1*k2^4*l4 + 2*i1^7*l2*l4 + i1^5*i2*l2*l4 + i1^3*i2^2*l2*l4 + 5*i1*i2^3*l2*l4 +
            4*i1^5*j2*l2*l4 + 4*i1^3*i2*j2*l2*l4 + i1^3*j2^2*l2*l4 + 6*i1*i2*j2^2*l2*l4 + 3*i1*j2^3*l2*l4 +
            5*i1^5*k2*l2*l4 + 3*i1*i2^2*k2*l2*l4 + 2*i1^3*j2*k2*l2*l4 + i1*j2^2*k2*l2*l4 + 3*i1^3*k2^2*l2*l4 +
            5*i1*i2*k2^2*l2*l4 + 4*i1*j2*k2^2*l2*l4 + 6*i1*k2^3*l2*l4 + i1^5*l2^2*l4 + 6*i1^3*i2*l2^2*l4 +
            6*i1*i2^2*l2^2*l4 + 6*i1^3*j2*l2^2*l4 + 6*i1*i2*j2*l2^2*l4 + 3*i1*i2*k2*l2^2*l4 + 6*i1*j2*k2*l2^2*l4 +
            5*i1*k2^2*l2^2*l4 + 3*i1^3*l2^3*l4 + 6*i1*j2*l2^3*l4 + i1*k2*l2^3*l4 + 6*i1*l2^4*l4 + 3*i1^6*j3*l4 +
            i1^4*i2*j3*l4 + 5*i2^3*j3*l4 + i1^4*j2*j3*l4 + i1^2*i2*j2*j3*l4 + 6*i2^2*j2*j3*l4 + 6*i1^2*j2^2*j3*l4 +
            i2*j2^2*j3*l4 + 2*j2^3*j3*l4 + 2*i1^4*k2*j3*l4 + 6*i1^2*i2*k2*j3*l4 + 3*i2^2*k2*j3*l4 + 3*i1^2*j2*k2*j3*l4 +
            4*j2^2*k2*j3*l4 + 2*i1^2*k2^2*j3*l4 + 2*i2*k2^2*j3*l4 + 3*j2*k2^2*j3*l4 + 5*k2^3*j3*l4 + 4*i1^2*i2*l2*j3*l4 +
            4*i2^2*l2*j3*l4 + 2*i1^2*j2*l2*j3*l4 + 5*i2*j2*l2*j3*l4 + 5*j2^2*l2*j3*l4 + 4*i1^2*k2*l2*j3*l4 +
            2*i2*k2*l2*j3*l4 + 3*k2^2*l2*j3*l4 + 4*i2*l2^2*j3*l4 + 3*j2*l2^2*j3*l4 + 5*k2*l2^2*j3*l4 + l2^3*j3*l4 +
            3*i1^3*j3^2*l4 + 5*i1*i2*j3^2*l4 + i1*k2*j3^2*l4 + j3^3*l4 + 2*i1^6*k3*l4 + 3*i1^4*i2*k3*l4 +
            6*i1^2*i2^2*k3*l4 + 3*i1^4*j2*k3*l4 + i1^2*i2*j2*k3*l4 + i1^4*k2*k3*l4 + i1^2*i2*k2*k3*l4 + 3*i2^2*k2*k3*l4 +
            4*i2*j2*k2*k3*l4 + 2*i1^2*k2^2*k3*l4 + j2*k2^2*k3*l4 + 3*k2^3*k3*l4 + 5*i1^4*l2*k3*l4 + 3*i1^2*i2*l2*k3*l4 +
            i2^2*l2*k3*l4 + 4*i1^2*j2*l2*k3*l4 + 3*i2*j2*l2*k3*l4 + 3*j2^2*l2*k3*l4 + i1^2*k2*l2*k3*l4 + 6*i2*k2*l2*k3*l4
            + 6*k2^2*l2*k3*l4 + i1^2*l2^2*k3*l4 + 2*i2*l2^2*k3*l4 + j2*l2^2*k3*l4 + 2*k2*l2^2*k3*l4 + 5*l2^3*k3*l4 +
            4*i1*i2*j3*k3*l4 + 2*i1*l2*j3*k3*l4 + 4*j3^2*k3*l4 + i1^3*k3^2*l4 + i1*i2*k3^2*l4 + i1*j2*k3^2*l4 +
            5*i1*k2*k3^2*l4 + 3*i1*l2*k3^2*l4 + 4*j3*k3^2*l4 + k3^3*l4 + 6*i1^6*l3*l4 + 2*i1^4*i2*l3*l4 +
            2*i1^2*i2^2*l3*l4 + 5*i2^3*l3*l4 + 5*i1^4*j2*l3*l4 + 6*i2*j2^2*l3*l4 + 3*j2^3*l3*l4 + 5*i1^4*k2*l3*l4 +
            5*i1^2*i2*k2*l3*l4 + 5*i2^2*k2*l3*l4 + 5*i1^2*j2*k2*l3*l4 + 4*j2^2*k2*l3*l4 + i2*k2^2*l3*l4 + 3*k2^3*l3*l4 +
            i1^4*l2*l3*l4 + i1^2*i2*l2*l3*l4 + 4*i2^2*l2*l3*l4 + 2*i1^2*j2*l2*l3*l4 + 6*i2*j2*l2*l3*l4 + 3*j2^2*l2*l3*l4 +
            4*i1^2*k2*l2*l3*l4 + 4*i2*k2*l2*l3*l4 + 6*j2*k2*l2*l3*l4 + 6*k2^2*l2*l3*l4 + i1^2*l2^2*l3*l4 + 2*j2*l2^2*l3*l4
            + 2*k2*l2^2*l3*l4 + l2^3*l3*l4 + 5*i1^3*j3*l3*l4 + 5*i1*i2*j3*l3*l4 + 6*i1*j2*j3*l3*l4 + 5*i1*k2*j3*l3*l4 +
            4*i1*l2*j3*l3*l4 + 3*j3^2*l3*l4 + 4*i1^3*l3^2*l4 + 5*i1*i2*l3^2*l4 + i1*j2*l3^2*l4 + 4*i1*k2*l3^2*l4 +
            3*j3*l3^2*l4 + 2*l3^3*l4 + 5*i1^5*l4^2 + 4*i1^3*i2*l4^2 + 6*i1*i2^2*l4^2 + 4*i1^3*j2*l4^2 + 5*i1*i2*j2*l4^2 +
            3*i1*j2^2*l4^2 + i1^3*k2*l4^2 + 2*i1*i2*k2*l4^2 + i1*k2^2*l4^2 + 2*i1*i2*l2*l4^2 + 5*i1*j2*l2*l4^2 +
            2*i1*l2^2*l4^2 + i1^2*j3*l4^2 + 3*i2*j3*l4^2 + 4*j2*j3*l4^2 + 6*k2*j3*l4^2 + 4*l2*j3*l4^2 + 3*i2*k3*l4^2 +
            6*j2*k3*l4^2 + 4*k2*k3*l4^2 + i1^2*l3*l4^2 + 3*i2*l3*l4^2 + 2*j2*l3*l4^2 + 5*k2*l3*l4^2 + l2*l3*l4^2 +
            4*i1^9*m4 + 5*i1^7*i2*m4 + 5*i1^5*i2^2*m4 + i1^3*i2^3*m4 + 4*i1*i2^4*m4 + 6*i1^7*j2*m4 + i1^5*i2*j2*m4 +
            4*i1^3*i2^2*j2*m4 + 5*i1*i2^3*j2*m4 + 6*i1^5*j2^2*m4 + 3*i1^3*i2*j2^2*m4 + 3*i1*i2^2*j2^2*m4 + 6*i1^3*j2^3*m4
            + 5*i1*i2*j2^3*m4 + 4*i1*j2^4*m4 + 4*i1^7*k2*m4 + 3*i1^5*i2*k2*m4 + 3*i1*i2^3*k2*m4 + 6*i1^5*j2*k2*m4 +
            3*i1^3*i2*j2*k2*m4 + 3*i1*i2^2*j2*k2*m4 + 4*i1^3*j2^2*k2*m4 + 6*i1*i2*j2^2*k2*m4 + 2*i1*j2^3*k2*m4 +
            2*i1^5*k2^2*m4 + 4*i1*i2^2*k2^2*m4 + 4*i1^3*j2*k2^2*m4 + 4*i1*i2*j2*k2^2*m4 + 2*i1*j2^2*k2^2*m4 +
            3*i1*i2*k2^3*m4 + 4*i1*j2*k2^3*m4 + 2*i1*k2^4*m4 + 4*i1^7*l2*m4 + i1^5*i2*l2*m4 + 4*i1^3*i2^2*l2*m4 +
            2*i1*i2^3*l2*m4 + 2*i1^5*j2*l2*m4 + 5*i1^3*i2*j2*l2*m4 + 4*i1*i2^2*j2*l2*m4 + 6*i1^3*j2^2*l2*m4 +
            i1*j2^3*l2*m4 + 2*i1^5*k2*l2*m4 + i1^3*i2*k2*l2*m4 + i1*i2^2*k2*l2*m4 + 4*i1^3*j2*k2*l2*m4 +
            4*i1*j2^2*k2*l2*m4 + 6*i1^3*k2^2*l2*m4 + 5*i1*i2*k2^2*l2*m4 + 4*i1*j2*k2^2*l2*m4 + 2*i1^5*l2^2*m4 +
            4*i1^3*i2*l2^2*m4 + i1*i2^2*l2^2*m4 + 5*i1^3*j2*l2^2*m4 + 6*i1*i2*j2*l2^2*m4 + 2*i1^3*k2*l2^2*m4 +
            4*i1*i2*k2*l2^2*m4 + 6*i1*j2*k2*l2^2*m4 + 4*i1*k2^2*l2^2*m4 + 4*i1^3*l2^3*m4 + 2*i1*i2*l2^3*m4 + i1*j2*l2^3*m4
            + 2*i1*k2*l2^3*m4 + 4*i1^4*i2*i3*m4 + 3*i1^2*i2^2*i3*m4 + 6*i1^4*j2*i3*m4 + i1^2*i2*j2*i3*m4 +
            3*i1^2*j2^2*i3*m4 + 2*i1^4*k2*i3*m4 + 2*i2^2*k2*i3*m4 + 5*i1^2*j2*k2*i3*m4 + 3*i2*j2*k2*i3*m4 +
            2*j2^2*k2*i3*m4 + 4*i1^2*k2^2*i3*m4 + 5*j2*k2^2*i3*m4 + k2^3*i3*m4 + 2*i1^4*l2*i3*m4 + 4*i1^2*i2*l2*i3*m4 +
            4*i2^2*l2*i3*m4 + 4*i1^2*j2*l2*i3*m4 + 6*i2*j2*l2*i3*m4 + 4*j2^2*l2*i3*m4 + 5*i2*k2*l2*i3*m4 +
            3*j2*k2*l2*i3*m4 + 3*k2^2*l2*i3*m4 + i1^2*l2^2*i3*m4 + 5*i2*l2^2*i3*m4 + 2*j2*l2^2*i3*m4 + 2*l2^3*i3*m4 +
            i1^3*i3^2*m4 + 3*i1*k2*i3^2*m4 + 4*i1*l2*i3^2*m4 + 2*i1^6*j3*m4 + i1^4*i2*j3*m4 + 6*i1^2*i2^2*j3*m4 +
            6*i2^3*j3*m4 + 2*i1^4*j2*j3*m4 + i1^2*i2*j2*j3*m4 + 3*i2^2*j2*j3*m4 + 4*i2*j2^2*j3*m4 + j2^3*j3*m4 +
            2*i1^4*k2*j3*m4 + 4*i1^2*i2*k2*j3*m4 + 2*i2^2*k2*j3*m4 + 2*i1^2*j2*k2*j3*m4 + 4*i2*j2*k2*j3*m4 + j2^2*k2*j3*m4
            + 2*i1^2*k2^2*j3*m4 + 6*i2*k2^2*j3*m4 + 4*j2*k2^2*j3*m4 + 6*k2^3*j3*m4 + 6*i1^4*l2*j3*m4 + 6*i1^2*i2*l2*j3*m4
            + 2*i1^2*j2*l2*j3*m4 + i2*j2*l2*j3*m4 + 6*j2^2*l2*j3*m4 + 5*i1^2*k2*l2*j3*m4 + 6*i2*k2*l2*j3*m4 +
            j2*k2*l2*j3*m4 + 4*k2^2*l2*j3*m4 + 2*i2*l2^2*j3*m4 + 5*j2*l2^2*j3*m4 + 6*k2*l2^2*j3*m4 + l2^3*j3*m4 +
            6*i1^3*j3^2*m4 + 6*i1*i2*j3^2*m4 + 5*i1*j2*j3^2*m4 + 3*i1*k2*j3^2*m4 + 4*i1*l2*j3^2*m4 + 3*i1^6*l3*m4 +
            3*i1^4*i2*l3*m4 + 2*i1^2*i2^2*l3*m4 + i2^3*l3*m4 + 3*i1^4*j2*l3*m4 + 3*i1^2*i2*j2*l3*m4 + 2*i2^2*j2*l3*m4 +
            i1^2*j2^2*l3*m4 + 4*j2^3*l3*m4 + 4*i1^4*k2*l3*m4 + 2*i1^2*i2*k2*l3*m4 + 3*i1^2*j2*k2*l3*m4 + 2*i2*j2*k2*l3*m4
            + 5*j2^2*k2*l3*m4 + 4*i1^2*k2^2*l3*m4 + 3*i2*k2^2*l3*m4 + 4*j2*k2^2*l3*m4 + 2*k2^3*l3*m4 + 5*i1^4*l2*l3*m4 +
            2*i2^2*l2*l3*m4 + 5*i2*j2*l2*l3*m4 + j2^2*l2*l3*m4 + 3*i1^2*k2*l2*l3*m4 + 5*i2*k2*l2*l3*m4 + 5*j2*k2*l2*l3*m4
            + 4*k2^2*l2*l3*m4 + 4*j2*l2^2*l3*m4 + 4*k2*l2^2*l3*m4 + 4*l2^3*l3*m4 + 2*i1^3*i3*l3*m4 + 5*i1*i2*i3*l3*m4 +
            6*i1*j2*i3*l3*m4 + 5*i1*k2*i3*l3*m4 + 6*i1*l2*i3*l3*m4 + 2*i3^2*l3*m4 + 5*i1^3*j3*l3*m4 + 2*i1*i2*j3*l3*m4 +
            5*i1*j2*j3*l3*m4 + i1*k2*j3*l3*m4 + 2*i1*l2*j3*l3*m4 + 2*j3^2*l3*m4 + 3*i1*i2*l3^2*m4 + 6*i1*k2*l3^2*m4 +
            6*i1*l2*l3^2*m4 + 3*i3*l3^2*m4 + 5*j3*l3^2*m4 + 5*l3^3*m4 + i1^8*i5 + i1^6*i2*i5 + 6*i1^6*j2*i5 + 6*i1^6*k2*i5
            + 2*i1^4*i2*k2*i5 + 5*i1^4*j2*k2*i5 + 5*i1^4*k2^2*i5 + 4*i1^2*i2*k2^2*i5 + 3*i1^2*j2*k2^2*i5 + 3*i1^2*k2^3*i5
            + 3*k2^4*i5 + 5*i1^6*l2*i5 + 2*i1^4*i2*l2*i5 + 5*i1^4*j2*l2*i5 + 3*i1^4*k2*l2*i5 + 6*i1^2*i2*k2*l2*i5 +
            i1^2*j2*k2*l2*i5 + 5*i1^2*k2^2*l2*i5 + 6*i2*k2^2*l2*i5 + j2*k2^2*l2*i5 + 6*k2^3*l2*i5 + 5*i1^4*l2^2*i5 +
            3*i1^2*i2*l2^2*i5 + 4*i1^2*j2*l2^2*i5 + 5*i1^2*k2*l2^2*i5 + 2*i2*k2*l2^2*i5 + 5*j2*k2*l2^2*i5 + 4*k2^2*l2^2*i5
            + 4*i1^2*l2^3*i5 + 5*i2*l2^3*i5 + 2*j2*l2^3*i5 + 3*k2*l2^3*i5 + 4*l2^4*i5 + 2*i1^3*k2*j3*i5 + 6*i1*k2^2*j3*i5
            + 6*i1^3*l2*j3*i5 + 4*i1*k2*l2*j3*i5 + 2*i1*l2^2*j3*i5 + 4*i1^8*j5 + i1^6*i2*j5 + 3*i1^4*i2^2*j5 +
            6*i1^2*i2^3*j5 + i2^4*j5 + 5*i1^6*j2*j5 + 2*i1^4*i2*j2*j5 + 6*i1^2*i2^2*j2*j5 + 3*i2^3*j2*j5 + 4*i1^4*j2^2*j5
            + 5*i1^2*i2*j2^2*j5 + 6*i2^2*j2^2*j5 + 4*i1^2*j2^3*j5 + 3*i2*j2^3*j5 + j2^4*j5 + 2*i1^6*k2*j5 +
            5*i1^4*i2*k2*j5 + i1^2*i2^2*k2*j5 + 5*i2^3*k2*j5 + 6*i1^4*j2*k2*j5 + 5*i1^2*i2*j2*k2*j5 + 4*i2^2*j2*k2*j5 +
            6*i1^2*j2^2*k2*j5 + 5*i2*j2^2*k2*j5 + 6*i1^4*k2^2*j5 + 6*i1^2*i2*k2^2*j5 + 6*i2^2*k2^2*j5 + 6*i1^2*j2*k2^2*j5
            + 5*i2*j2*k2^2*j5 + 4*j2^2*k2^2*j5 + 6*i1^2*k2^3*j5 + 3*j2*k2^3*j5 + 4*k2^4*j5 + 2*i1^6*l2*j5 +
            5*i1^4*i2*l2*j5 + i1^2*i2^2*l2*j5 + 6*i2^3*l2*j5 + 5*i1^4*j2*l2*j5 + 6*i1^2*i2*j2*l2*j5 + 6*i2^2*j2*l2*j5 +
            6*i1^2*j2^2*l2*j5 + 5*i2*j2^2*l2*j5 + 4*j2^3*l2*j5 + 5*i1^4*k2*l2*j5 + 5*i1^2*i2*k2*l2*j5 + 2*i2^2*k2*l2*j5 +
            3*i1^2*j2*k2*l2*j5 + 2*i2*j2*k2*l2*j5 + 2*j2^2*k2*l2*j5 + 4*i1^2*k2^2*l2*j5 + i2*k2^2*l2*j5 + j2*k2^2*l2*j5 +
            4*k2^3*l2*j5 + 6*i1^4*l2^2*j5 + 2*i1^2*i2*l2^2*j5 + 3*i2^2*l2^2*j5 + 5*i2*j2*l2^2*j5 + 5*j2^2*l2^2*j5 +
            6*i1^2*k2*l2^2*j5 + i2*k2*l2^2*j5 + 4*j2*k2*l2^2*j5 + 3*k2^2*l2^2*j5 + 3*i1^2*l2^3*j5 + 4*i2*l2^3*j5 +
            6*j2*l2^3*j5 + 3*k2*l2^3*j5 + 2*l2^4*j5 + 6*i1^5*i3*j5 + 6*i1^3*i2*i3*j5 + 4*i1^3*j2*i3*j5 + 6*i1^3*k2*i3*j5 +
            3*i1*i2*k2*i3*j5 + 2*i1*k2^2*i3*j5 + i1^3*l2*i3*j5 + 5*i1*i2*l2*i3*j5 + 2*i1*k2*l2*i3*j5 + 2*i1*l2^2*i3*j5 +
            i1^2*i3^2*j5 + 6*k2*i3^2*j5 + l2*i3^2*j5 + 3*i1^5*j3*j5 + i1^3*i2*j3*j5 + i1*i2^2*j3*j5 + 5*i1^3*j2*j3*j5 +
            6*i1*j2^2*j3*j5 + i1^3*k2*j3*j5 + 5*i1*k2^2*j3*j5 + 5*i1*i2*l2*j3*j5 + 6*i1*j2*l2*j3*j5 + 4*i1*l2^2*j3*j5 +
            2*i1^2*i3*j3*j5 + 5*i2*i3*j3*j5 + 2*j2*i3*j3*j5 + 3*k2*i3*j3*j5 + 2*l2*i3*j3*j5 + 5*i1^2*j3^2*j5 +
            3*i2*j3^2*j5 + j2*j3^2*j5 + 5*k2*j3^2*j5 + 5*l2*j3^2*j5 + 6*i1^5*k3*j5 + 2*i1^3*i2*k3*j5 + 3*i1*i2^2*k3*j5 +
            4*i1^3*j2*k3*j5 + 5*i1*i2*j2*k3*j5 + 6*i1*j2^2*k3*j5 + 6*i1^3*k2*k3*j5 + 2*i1*j2*k2*k3*j5 + 4*i1*k2^2*k3*j5 +
            i1*i2*l2*k3*j5 + 6*i1*j2*l2*k3*j5 + 4*i1*k2*l2*k3*j5 + i1*l2^2*k3*j5 + 6*i1^2*j3*k3*j5 + 3*i2*j3*k3*j5 +
            6*j2*j3*k3*j5 + 5*k2*j3*k3*j5 + 6*l2*j3*k3*j5 + 4*i2*k3^2*j5 + 5*j2*k3^2*j5 + 3*k2*k3^2*j5 + 4*l2*k3^2*j5 +
            6*i1^5*l3*j5 + 3*i1^3*i2*l3*j5 + i1*i2^2*l3*j5 + i1^3*j2*l3*j5 + i1*i2*j2*l3*j5 + 3*i1*j2^2*l3*j5 +
            6*i1^3*k2*l3*j5 + i1*j2*k2*l3*j5 + 6*i1*k2^2*l3*j5 + 4*i1^3*l2*l3*j5 + 5*i1*j2*l2*l3*j5 + 5*i1*l2^2*l3*j5 +
            6*i1^2*i3*l3*j5 + 5*i2*i3*l3*j5 + 6*j2*i3*l3*j5 + 4*k2*i3*l3*j5 + 6*l2*i3*l3*j5 + 6*i1^2*l3^2*j5 +
            4*i2*l3^2*j5 + 4*j2*l3^2*j5 + 6*k2*l3^2*j5 + 6*l2*l3^2*j5 + 2*i1^4*i4*j5 + 5*i1^2*i2*i4*j5 + 4*i2^2*i4*j5 +
            4*i2*j2*i4*j5 + 6*j2^2*i4*j5 + i1^2*k2*i4*j5 + 5*i2*k2*i4*j5 + 2*j2*k2*i4*j5 + 5*k2^2*i4*j5 + 4*i1^2*l2*i4*j5
            + 5*i2*l2*i4*j5 + k2*l2*i4*j5 + 2*l2^2*i4*j5 + i1^4*l4*j5 + 6*i1^2*i2*l4*j5 + 2*i2^2*l4*j5 + 4*i1^2*j2*l4*j5 +
            6*i2*j2*l4*j5 + 6*j2^2*l4*j5 + 3*j2*k2*l4*j5 + k2^2*l4*j5 + 2*i1^2*l2*l4*j5 + 6*i2*l2*l4*j5 + 4*j2*l2*l4*j5 +
            l2^2*l4*j5 + 5*i1^8*k5 + 5*i1^6*i2*k5 + 6*i1^4*i2^2*k5 + 2*i1^4*i2*j2*k5 + 6*i1^4*j2^2*k5 + 5*i1^6*k2*k5 +
            2*i1^4*i2*k2*k5 + i1^2*i2^2*k2*k5 + 6*i1^4*j2*k2*k5 + 5*i1^2*i2*j2*k2*k5 + i1^2*j2^2*k2*k5 + 4*i1^4*k2^2*k5 +
            i1^2*i2*k2^2*k5 + 2*i2^2*k2^2*k5 + 6*i1^2*j2*k2^2*k5 + 3*i2*j2*k2^2*k5 + 2*j2^2*k2^2*k5 + i1^2*k2^3*k5 +
            4*j2*k2^3*k5 + 5*k2^4*k5 + 3*i1^6*l2*k5 + 6*i1^4*i2*l2*k5 + 4*i1^2*i2^2*l2*k5 + 2*i1^4*j2*l2*k5 +
            6*i1^2*i2*j2*l2*k5 + 4*i1^2*j2^2*l2*k5 + 6*i1^4*k2*l2*k5 + 4*i1^2*j2*k2*l2*k5 + 4*i1^2*k2^2*l2*k5 +
            2*i2*k2^2*l2*k5 + 2*j2*k2^2*l2*k5 + 5*k2^3*l2*k5 + i1^4*l2^2*k5 + 4*i1^2*i2*l2^2*k5 + 3*i2^2*l2^2*k5 +
            2*i1^2*j2*l2^2*k5 + i2*j2*l2^2*k5 + 3*j2^2*l2^2*k5 + i1^2*k2*l2^2*k5 + 2*i2*k2*l2^2*k5 + 4*j2*k2*l2^2*k5 +
            3*k2^2*l2^2*k5 + i1^2*l2^3*k5 + 4*i2*l2^3*k5 + 2*j2*l2^3*k5 + l2^4*k5 + i1^5*j3*k5 + 3*i1^3*i2*j3*k5 +
            4*i1^3*j2*j3*k5 + 3*i1^3*k2*j3*k5 + 5*i1*i2*k2*j3*k5 + 2*i1*j2*k2*j3*k5 + 4*i1*k2^2*j3*k5 + 2*i1^3*l2*j3*k5 +
            6*i1*i2*l2*j3*k5 + i1*j2*l2*j3*k5 + 3*i1*k2*l2*j3*k5 + 6*i1*l2^2*j3*k5 + 6*k2*j3^2*k5 + 4*l2*j3^2*k5 +
            3*i1^3*i2*k3*k5 + 4*i1^3*j2*k3*k5 + 6*i1^3*k2*k3*k5 + 2*i1*i2*k2*k3*k5 + 5*i1*j2*k2*k3*k5 + 2*i1*k2^2*k3*k5 +
            2*i1^3*l2*k3*k5 + 2*i1*i2*l2*k3*k5 + 5*i1*j2*l2*k3*k5 + 2*i1*k2*l2*k3*k5 + 5*i1*l2^2*k3*k5 + 5*i1^2*k3^2*k5 +
            2*k2*k3^2*k5 + 2*l2*k3^2*k5 + i1^5*l3*k5 + 4*i1*i2^2*l3*k5 + 6*i1^3*j2*l3*k5 + 6*i1*i2*j2*l3*k5 +
            4*i1*j2^2*l3*k5 + i1^3*k2*l3*k5 + i1*i2*k2*l3*k5 + 5*i1*j2*k2*l3*k5 + 2*i1*k2^2*l3*k5 + 6*i1^3*l2*l3*k5 +
            3*i1*i2*l2*l3*k5 + 3*i1*j2*l2*l3*k5 + 5*i1*k2*l2*l3*k5 + 5*i1*l2^2*l3*k5 + 6*i2*j3*l3*k5 + j2*j3*l3*k5 +
            5*k2*j3*l3*k5 + 6*l2*j3*l3*k5 + 2*i1^2*i2*l4*k5 + 5*i1^2*j2*l4*k5 + 4*i1^2*k2*l4*k5 + i2*k2*l4*k5 +
            6*j2*k2*l4*k5 + 3*k2^2*l4*k5 + 5*i1^2*l2*l4*k5 + 3*i2*l2*l4*k5 + 4*j2*l2*l4*k5 + k2*l2*l4*k5 + l2^2*l4*k5 +
            3*i1*k3*l4*k5 + 4*i1*l3*l4*k5 + 3*l4^2*k5 + 4*i1^8*l5 + 4*i1^4*i2^2*l5 + i1^2*i2^3*l5 + 6*i1^6*j2*l5 +
            4*i1^4*i2*j2*l5 + 4*i1^2*i2^2*j2*l5 + 6*i1^4*j2^2*l5 + 3*i1^2*i2*j2^2*l5 + 6*i1^2*j2^3*l5 + i1^6*k2*l5 +
            3*i1^4*i2*k2*l5 + 3*i1^2*i2^2*k2*l5 + i2^3*k2*l5 + i1^4*j2*k2*l5 + 3*i1^2*i2*j2*k2*l5 + 4*i2^2*j2*k2*l5 +
            i1^2*j2^2*k2*l5 + 3*i2*j2^2*k2*l5 + 6*j2^3*k2*l5 + i1^4*k2^2*l5 + i1^2*i2*k2^2*l5 + 4*i2^2*k2^2*l5 +
            3*j2^2*k2^2*l5 + i1^2*k2^3*l5 + 3*j2*k2^3*l5 + 4*k2^4*l5 + 6*i1^6*l2*l5 + i1^2*i2^2*l2*l5 + 2*i2^3*l2*l5 +
            i1^4*j2*l2*l5 + 3*i1^2*i2*j2*l2*l5 + i2^2*j2*l2*l5 + 3*i1^2*j2^2*l2*l5 + 6*i2*j2^2*l2*l5 + 5*j2^3*l2*l5 +
            2*i1^4*k2*l2*l5 + 6*i2^2*k2*l2*l5 + 3*i1^2*j2*k2*l2*l5 + 3*i2*j2*k2*l2*l5 + 5*j2^2*k2*l2*l5 + i1^2*k2^2*l2*l5
            + i2*k2^2*l2*l5 + 6*j2*k2^2*l2*l5 + 6*k2^3*l2*l5 + 6*i1^2*i2*l2^2*l5 + 3*i2^2*l2^2*l5 + 2*i1^2*j2*l2^2*l5 +
            2*i2*j2*l2^2*l5 + 2*j2^2*l2^2*l5 + 2*i1^2*k2*l2^2*l5 + 3*i2*k2*l2^2*l5 + 5*j2*k2*l2^2*l5 + 4*k2^2*l2^2*l5 +
            5*i1^2*l2^3*l5 + 2*i2*l2^3*l5 + 5*j2*l2^3*l5 + 5*k2*l2^3*l5 + 5*l2^4*l5 + 2*i1^5*i3*l5 + 5*i1^3*k2*i3*l5 +
            2*i1*i2*k2*i3*l5 + 5*i1*j2*k2*i3*l5 + 5*i1*k2^2*i3*l5 + 2*i1^3*l2*i3*l5 + 3*i1*i2*l2*i3*l5 + 4*i1*j2*l2*i3*l5
            + 3*i1*k2*l2*i3*l5 + i1*l2^2*i3*l5 + 6*i1^5*j3*l5 + i1*i2^2*j3*l5 + 5*i1^3*j2*j3*l5 + 5*i1*i2*j2*j3*l5 +
            i1*j2^2*j3*l5 + i1^3*k2*j3*l5 + 3*i1*i2*k2*j3*l5 + 3*i1*j2*k2*j3*l5 + i1*k2^2*j3*l5 + 3*i1^3*l2*j3*l5 +
            4*i1*i2*l2*j3*l5 + 2*i1*j2*l2*j3*l5 + i1*k2*l2*j3*l5 + 6*i1*l2^2*j3*l5 + i1^2*j3^2*l5 + 5*i2*j3^2*l5 +
            2*j2*j3^2*l5 + 3*k2*j3^2*l5 + 2*l2*j3^2*l5 + 4*i1^5*k3*l5 + 4*i1^3*i2*k3*l5 + 5*i1*i2^2*k3*l5 +
            6*i1^3*j2*k3*l5 + 4*i1*i2*j2*k3*l5 + 5*i1*j2^2*k3*l5 + 2*i1^3*k2*k3*l5 + i1*i2*k2*k3*l5 + 2*i1*j2*k2*k3*l5 +
            i1*k2^2*k3*l5 + 5*i1*i2*l2*k3*l5 + 6*i1*j2*l2*k3*l5 + 6*i1*l2^2*k3*l5 + 2*i1^2*j3*k3*l5 + 3*i2*j3*k3*l5 +
            4*j2*j3*k3*l5 + 5*k2*j3*k3*l5 + 6*l2*j3*k3*l5 + 6*i1^2*k3^2*l5 + k2*k3^2*l5 + 3*l2*k3^2*l5 + i1^5*l3*l5 +
            2*i1^3*i2*l3*l5 + 3*i1*i2^2*l3*l5 + 3*i1^3*j2*l3*l5 + 4*i1*i2*j2*l3*l5 + 4*i1^3*k2*l3*l5 + i1*i2*k2*l3*l5 +
            2*i1*j2*k2*l3*l5 + i1*k2^2*l3*l5 + i1^3*l2*l3*l5 + 3*i1*j2*l2*l3*l5 + 2*i1^2*i3*l3*l5 + 3*i2*i3*l3*l5 +
            4*j2*i3*l3*l5 + 3*k2*i3*l3*l5 + 6*l2*i3*l3*l5 + 2*i1^2*j3*l3*l5 + 3*i2*j3*l3*l5 + j2*j3*l3*l5 + 3*k2*j3*l3*l5
            + 6*l2*j3*l3*l5 + 6*i1^2*l3^2*l5 + i2*l3^2*l5 + 3*j2*l3^2*l5 + 5*k2*l3^2*l5 + 6*l2*l3^2*l5 + i1^4*k4*l5 +
            4*i1^2*k2*k4*l5 + 3*k2^2*k4*l5 + 2*i1^2*l2*k4*l5 + 3*l2^2*k4*l5 + 3*i1^4*m4*l5 + 5*i1^2*i2*m4*l5 + i2^2*m4*l5
            + 2*i1^2*j2*m4*l5 + 5*i2*j2*m4*l5 + j2^2*m4*l5 + 4*i1^2*k2*m4*l5 + 2*i2*k2*m4*l5 + j2*k2*m4*l5 + 6*k2^2*m4*l5
            + i1^2*l2*m4*l5 + 6*i2*l2*m4*l5 + j2*l2*m4*l5 + 4*i1*i3*m4*l5 + 3*i1*l3*m4*l5 + 6*i1^3*l5^2 + 5*i1*i2*l5^2 +
            2*i1*j2*l5^2 + 4*i1*k2*l5^2 + i1*l2*l5^2 + 6*j3*l5^2 + 6*k3*l5^2 + 3*l3*l5^2 + 6*i1^5*k2*i6 + 3*i1^3*k2^2*i6 +
            4*i1*k2^3*i6 + 4*i1^5*l2*i6 + 6*i1^3*k2*l2*i6 + 2*i1*k2^2*l2*i6 + 4*i1^3*l2^2*i6 + 5*i1*k2*l2^2*i6 +
            5*i1*l2^3*i6 + 3*i1^2*k2*l3*i6 + 5*k2^2*l3*i6 + 2*i1^2*l2*l3*i6 + 4*k2*l2*l3*i6 + 6*l2^2*l3*i6;
        Append(~JI, J13); Append(~Wght, 13);
    end if;
    if degmax le 13 then return JI, Wght; end if;

    /* Degree 14 */
    if degmin le 14 then
        J14:= 2*i1^14 + i1^12*i2 + 2*i1^10*i2^2 + 2*i1^8*i2^3 + 2*i1^6*i2^4 + 2*i1^4*i2^5 + i1^2*i2^6 + 6*i2^7 +
            4*i1^12*j2 + 2*i1^10*i2*j2 + 2*i1^8*i2^2*j2 + 4*i1^6*i2^3*j2 + 6*i1^4*i2^4*j2 + i1^2*i2^5*j2 + i1^10*j2^2 +
            i1^8*i2*j2^2 + 5*i1^6*i2^2*j2^2 + 5*i1^4*i2^3*j2^2 + i1^2*i2^4*j2^2 + 5*i1^6*i2*j2^3 + 6*i1^4*i2^2*j2^3 +
            i1^2*i2^3*j2^3 + 5*i1^6*j2^4 + 2*i1^4*i2*j2^4 + i1^2*i2^2*j2^4 + i1^2*i2*j2^5 + i1^2*j2^6 + j2^7 + 6*i1^12*k2
            + 5*i1^10*i2*k2 + i1^8*i2^2*k2 + 5*i1^6*i2^3*k2 + 2*i1^4*i2^4*k2 + 5*i2^6*k2 + 2*i1^10*j2*k2 + 3*i1^8*i2*j2*k2
            + 6*i1^6*i2^2*j2*k2 + 4*i1^4*i2^3*j2*k2 + i1^2*i2^4*j2*k2 + 5*i2^5*j2*k2 + 3*i1^8*j2^2*k2 + 4*i1^6*i2*j2^2*k2
            + i1^4*i2^2*j2^2*k2 + 3*i1^2*i2^3*j2^2*k2 + 5*i2^4*j2^2*k2 + 6*i1^6*j2^3*k2 + 6*i1^4*i2*j2^3*k2 +
            6*i1^2*i2^2*j2^3*k2 + 5*i2^3*j2^3*k2 + i1^4*j2^4*k2 + 3*i1^2*i2*j2^4*k2 + 5*i2^2*j2^4*k2 + i1^2*j2^5*k2 +
            5*i2*j2^5*k2 + 5*j2^6*k2 + 6*i1^10*k2^2 + i1^8*i2*k2^2 + 3*i1^6*i2^2*k2^2 + 5*i1^4*i2^3*k2^2 +
            6*i1^2*i2^4*k2^2 + 2*i2^5*k2^2 + 2*i1^8*j2*k2^2 + 6*i1^6*i2*j2*k2^2 + 5*i1^4*i2^2*j2*k2^2 + i1^2*i2^3*j2*k2^2
            + 6*i2^4*j2*k2^2 + i1^4*i2*j2^2*k2^2 + 5*i1^2*i2^2*j2^2*k2^2 + 5*i2^3*j2^2*k2^2 + 3*i1^4*j2^3*k2^2 +
            5*i1^2*i2*j2^3*k2^2 + 6*i2^2*j2^3*k2^2 + 4*i1^2*j2^4*k2^2 + 2*i2*j2^4*k2^2 + 4*i1^8*k2^3 + 6*i1^6*i2*k2^3 +
            i1^4*i2^2*k2^3 + 5*i1^2*i2^3*k2^3 + 4*i2^4*k2^3 + 5*i1^6*j2*k2^3 + 3*i1^4*i2*j2*k2^3 + 4*i1^2*i2^2*j2*k2^3 +
            3*i2^3*j2*k2^3 + 5*i1^4*j2^2*k2^3 + 3*i1^2*i2*j2^2*k2^3 + i2^2*j2^2*k2^3 + 5*i1^2*j2^3*k2^3 + i2*j2^3*k2^3 +
            5*j2^4*k2^3 + 3*i1^6*k2^4 + i1^4*i2*k2^4 + 2*i1^2*i2^2*k2^4 + 6*i2^3*k2^4 + 2*i1^4*j2*k2^4 + 3*i1^2*i2*j2*k2^4
            + 6*i2^2*j2*k2^4 + 3*i1^2*j2^2*k2^4 + 5*i2*j2^2*k2^4 + 2*j2^3*k2^4 + 5*i1^2*i2*k2^5 + 5*i2^2*k2^5 +
            6*i1^2*j2*k2^5 + i2*j2*k2^5 + j2^2*k2^5 + 5*i2*k2^6 + 5*k2^7 + 5*i1^10*i2*l2 + 6*i1^8*i2^2*l2 + 5*i1^6*i2^3*l2
            + 5*i1^4*i2^4*l2 + 5*i1^2*i2^5*l2 + i2^6*l2 + 2*i1^10*j2*l2 + 5*i1^8*i2*j2*l2 + i1^6*i2^2*j2*l2 +
            2*i1^2*i2^4*j2*l2 + i2^5*j2*l2 + 6*i1^6*i2*j2^2*l2 + 2*i1^4*i2^2*j2^2*l2 + 5*i1^2*i2^3*j2^2*l2 + i2^4*j2^2*l2
            + 4*i1^6*j2^3*l2 + 4*i1^4*i2*j2^3*l2 + i2^3*j2^3*l2 + 3*i1^4*j2^4*l2 + i1^2*i2*j2^4*l2 + i2^2*j2^4*l2 +
            i1^2*j2^5*l2 + i2*j2^5*l2 + j2^6*l2 + 4*i1^8*i2*k2*l2 + 6*i1^6*i2^2*k2*l2 + 6*i1^2*i2^4*k2*l2 + 6*i2^5*k2*l2 +
            3*i1^8*j2*k2*l2 + 4*i1^6*i2*j2*k2*l2 + 5*i1^4*i2^2*j2*k2*l2 + 6*i1^2*i2^3*j2*k2*l2 + 3*i2^4*j2*k2*l2 +
            i1^6*j2^2*k2*l2 + 2*i1^4*i2*j2^2*k2*l2 + 4*i1^2*i2^2*j2^2*k2*l2 + 5*i2^3*j2^2*k2*l2 + i1^4*j2^3*k2*l2 +
            6*i1^2*i2*j2^3*k2*l2 + 5*i2^2*j2^3*k2*l2 + 6*i1^2*j2^4*k2*l2 + 3*i2*j2^4*k2*l2 + 6*j2^5*k2*l2 + 3*i1^8*k2^2*l2
            + 6*i1^4*i2^2*k2^2*l2 + 5*i1^2*i2^3*k2^2*l2 + 6*i1^6*j2*k2^2*l2 + 3*i1^2*i2^2*j2*k2^2*l2 + 6*i2^3*j2*k2^2*l2 +
            3*i1^4*j2^2*k2^2*l2 + 4*i1^2*i2*j2^2*k2^2*l2 + 3*i2^2*j2^2*k2^2*l2 + 5*i1^2*j2^3*k2^2*l2 + 4*i2*j2^3*k2^2*l2 +
            j2^4*k2^2*l2 + 2*i1^4*i2*k2^3*l2 + 4*i1^2*i2^2*k2^3*l2 + 6*i2^3*k2^3*l2 + i1^4*j2*k2^3*l2 +
            4*i1^2*i2*j2*k2^3*l2 + i2^2*j2*k2^3*l2 + 5*i1^2*j2^2*k2^3*l2 + i2*j2^2*k2^3*l2 + 3*j2^3*k2^3*l2 +
            6*i1^4*k2^4*l2 + 6*i1^2*i2*k2^4*l2 + i2^2*k2^4*l2 + i1^2*j2*k2^4*l2 + 4*i2*j2*k2^4*l2 + 4*j2^2*k2^4*l2 +
            3*i1^2*k2^5*l2 + 3*i2*k2^5*l2 + 5*i1^10*l2^2 + 6*i1^8*i2*l2^2 + 2*i1^4*i2^3*l2^2 + i2^5*l2^2 + 5*i1^8*j2*l2^2
            + 3*i1^6*i2*j2*l2^2 + 5*i1^2*i2^3*j2*l2^2 + 2*i2^4*j2*l2^2 + 2*i1^6*j2^2*l2^2 + 4*i1^4*i2*j2^2*l2^2 +
            6*i1^2*i2^2*j2^2*l2^2 + 3*i2^3*j2^2*l2^2 + i1^4*j2^3*l2^2 + i1^2*i2*j2^3*l2^2 + 4*i2^2*j2^3*l2^2 +
            2*i1^2*j2^4*l2^2 + 5*i2*j2^4*l2^2 + 6*j2^5*l2^2 + 4*i1^6*i2*k2*l2^2 + 6*i1^4*i2^2*k2*l2^2 +
            5*i1^2*i2^3*k2*l2^2 + 4*i2^4*k2*l2^2 + i1^6*j2*k2*l2^2 + 5*i1^2*i2^2*j2*k2*l2^2 + 6*i2^3*j2*k2*l2^2 +
            2*i1^2*i2*j2^2*k2*l2^2 + 2*i2^2*j2^2*k2*l2^2 + 6*i1^2*j2^3*k2*l2^2 + 4*i2*j2^3*k2*l2^2 + 5*j2^4*k2*l2^2 +
            2*i1^6*k2^2*l2^2 + 6*i2^3*k2^2*l2^2 + i1^4*j2*k2^2*l2^2 + 4*i1^2*i2*j2*k2^2*l2^2 + 2*i2^2*j2*k2^2*l2^2 +
            2*i1^2*j2^2*k2^2*l2^2 + 5*i2*j2^2*k2^2*l2^2 + 3*j2^3*k2^2*l2^2 + 2*i1^2*i2*k2^3*l2^2 + 2*i1^2*j2*k2^3*l2^2 +
            2*i2*j2*k2^3*l2^2 + 5*i2*k2^4*l2^2 + 3*j2*k2^4*l2^2 + 2*k2^5*l2^2 + 5*i1^8*l2^3 + 6*i1^6*i2*l2^3 +
            4*i1^4*i2^2*l2^3 + 4*i1^2*i2^3*l2^3 + 4*i2^4*l2^3 + i1^6*j2*l2^3 + 3*i1^4*i2*j2*l2^3 + i2^3*j2*l2^3 +
            2*i1^4*j2^2*l2^3 + 4*i1^2*i2*j2^2*l2^3 + i2^2*j2^2*l2^3 + 6*i1^2*j2^3*l2^3 + j2^4*l2^3 + 5*i1^6*k2*l2^3 +
            4*i1^4*i2*k2*l2^3 + 6*i1^2*i2^2*k2*l2^3 + 6*i2^3*k2*l2^3 + i1^4*j2*k2*l2^3 + 6*i1^2*i2*j2*k2*l2^3 +
            4*i2^2*j2*k2*l2^3 + 6*i1^2*j2^2*k2*l2^3 + 2*i2*j2^2*k2*l2^3 + 3*j2^3*k2*l2^3 + 3*i1^4*k2^2*l2^3 +
            5*i1^2*i2*k2^2*l2^3 + 2*i2^2*k2^2*l2^3 + 6*i1^2*j2*k2^2*l2^3 + 2*j2^2*k2^2*l2^3 + 4*i1^2*k2^3*l2^3 +
            2*i2*k2^3*l2^3 + 6*k2^4*l2^3 + 3*i1^6*l2^4 + i1^4*i2*l2^4 + i1^2*i2^2*l2^4 + 3*i2^3*l2^4 + 4*i1^4*j2*l2^4 +
            5*i2^2*j2*l2^4 + 2*i1^2*j2^2*l2^4 + 3*i2*j2^2*l2^4 + j2^3*l2^4 + 3*i1^4*k2*l2^4 + 5*i1^2*i2*k2*l2^4 +
            6*i2^2*k2*l2^4 + 3*i1^2*j2*k2*l2^4 + 6*j2^2*k2*l2^4 + 5*i1^2*k2^2*l2^4 + 2*j2*k2^2*l2^4 + 5*k2^3*l2^4 +
            4*i1^4*l2^5 + 2*i1^2*i2*l2^5 + 6*i2^2*l2^5 + 6*i1^2*j2*l2^5 + i2*j2*l2^5 + 2*j2^2*l2^5 + i1^2*k2*l2^5 +
            3*i2*k2*l2^5 + 2*j2*k2*l2^5 + 3*k2^2*l2^5 + 4*i1^2*l2^6 + 5*i2*l2^6 + j2*l2^6 + 4*k2*l2^6 + 5*l2^7 +
            4*i1^11*i3 + 2*i1^9*i2*i3 + 5*i1^7*i2^2*i3 + 4*i1^5*i2^3*i3 + 3*i1^3*i2^4*i3 + 6*i1^9*j2*i3 + 2*i1^7*i2*j2*i3
            + 5*i1^5*i2^2*j2*i3 + 2*i1^3*i2^3*j2*i3 + 5*i1^7*j2^2*i3 + 6*i1^5*i2*j2^2*i3 + 4*i1^3*i2^2*j2^2*i3 +
            6*i1^5*j2^3*i3 + 2*i1^3*i2*j2^3*i3 + 3*i1^3*j2^4*i3 + 5*i1^9*k2*i3 + 2*i1^7*i2*k2*i3 + i1^5*i2^2*k2*i3 +
            3*i1*i2^4*k2*i3 + 5*i1^7*j2*k2*i3 + 4*i1^5*i2*j2*k2*i3 + i1^3*i2^2*j2*k2*i3 + 2*i1*i2^3*j2*k2*i3 +
            2*i1^5*j2^2*k2*i3 + 5*i1^3*i2*j2^2*k2*i3 + 4*i1*i2^2*j2^2*k2*i3 + i1^3*j2^3*k2*i3 + 2*i1*i2*j2^3*k2*i3 +
            3*i1*j2^4*k2*i3 + 5*i1^7*k2^2*i3 + 4*i1^3*i2^2*k2^2*i3 + 2*i1*i2^3*k2^2*i3 + 6*i1^5*j2*k2^2*i3 +
            5*i1^3*i2*j2*k2^2*i3 + 3*i1*i2^2*j2*k2^2*i3 + 5*i1^3*j2^2*k2^2*i3 + 2*i1*i2*j2^2*k2^2*i3 + i1^5*k2^3*i3 +
            4*i1^3*i2*k2^3*i3 + 3*i1*i2^2*k2^3*i3 + 2*i1*i2*j2*k2^3*i3 + 5*i1*j2^2*k2^3*i3 + 5*i1^3*k2^4*i3 +
            5*i1*i2*k2^4*i3 + 4*i1*j2*k2^4*i3 + 4*i1*k2^5*i3 + i1^9*l2*i3 + 6*i1^5*i2^2*l2*i3 + 6*i1^3*i2^3*l2*i3 +
            i1*i2^4*l2*i3 + 4*i1^7*j2*l2*i3 + 4*i1^5*i2*j2*l2*i3 + 6*i1^3*i2^2*j2*l2*i3 + 3*i1*i2^3*j2*l2*i3 +
            3*i1^5*j2^2*l2*i3 + 5*i1^3*i2*j2^2*l2*i3 + 6*i1*i2^2*j2^2*l2*i3 + 4*i1^3*j2^3*l2*i3 + 3*i1*i2*j2^3*l2*i3 +
            i1*j2^4*l2*i3 + 2*i1^7*k2*l2*i3 + 2*i1^5*i2*k2*l2*i3 + 2*i1*i2^3*k2*l2*i3 + i1^5*j2*k2*l2*i3 +
            4*i1^3*i2*j2*k2*l2*i3 + 2*i1*i2^2*j2*k2*l2*i3 + 2*i1^3*j2^2*k2*l2*i3 + 4*i1*i2*j2^2*k2*l2*i3 +
            6*i1*j2^3*k2*l2*i3 + 5*i1^5*k2^2*l2*i3 + 4*i1^3*i2*k2^2*l2*i3 + 3*i1*i2^2*k2^2*l2*i3 + i1^3*j2*k2^2*l2*i3 +
            i1*i2*j2*k2^2*l2*i3 + i1*j2^2*k2^2*l2*i3 + 3*i1^3*k2^3*l2*i3 + 6*i1*i2*k2^3*l2*i3 + 4*i1*k2^4*l2*i3 +
            i1^5*i2*l2^2*i3 + 6*i1^3*i2^2*l2^2*i3 + 5*i1*i2^3*l2^2*i3 + i1^5*j2*l2^2*i3 + 4*i1^3*i2*j2*l2^2*i3 +
            2*i1*i2^2*j2*l2^2*i3 + 5*i1^3*j2^2*l2^2*i3 + 2*i1*i2*j2^2*l2^2*i3 + 5*i1*j2^3*l2^2*i3 + 6*i1^5*k2*l2^2*i3 +
            i1^3*i2*k2*l2^2*i3 + 6*i1^3*j2*k2*l2^2*i3 + 2*i1*i2*j2*k2*l2^2*i3 + 3*i1^3*k2^2*l2^2*i3 + 5*i1*i2*k2^2*l2^2*i3
            + 3*i1*j2*k2^2*l2^2*i3 + 2*i1*k2^3*l2^2*i3 + 3*i1^3*i2*l2^3*i3 + 4*i1*i2^2*l2^3*i3 + 2*i1^3*j2*l2^3*i3 +
            3*i1*j2^2*l2^3*i3 + 6*i1*i2*k2*l2^3*i3 + 2*i1*j2*k2*l2^3*i3 + 6*i1*k2^2*l2^3*i3 + 3*i1^3*l2^4*i3 +
            i1*i2*l2^4*i3 + 4*i1*j2*l2^4*i3 + 2*i1*k2*l2^4*i3 + i1*l2^5*i3 + 4*i1^8*i3^2 + 2*i1^6*i2*i3^2 + i1^4*i2^2*i3^2
            + i1^6*j2*i3^2 + 5*i1^4*i2*j2*i3^2 + i1^4*j2^2*i3^2 + 5*i1^6*k2*i3^2 + 2*i1^4*i2*k2*i3^2 + 2*i1^2*i2^2*k2*i3^2
            + 3*i1^4*j2*k2*i3^2 + 3*i1^2*i2*j2*k2*i3^2 + 2*i1^2*j2^2*k2*i3^2 + 5*i1^4*k2^2*i3^2 + 6*i1^2*i2*k2^2*i3^2 +
            6*i2^2*k2^2*i3^2 + i1^2*j2*k2^2*i3^2 + 2*i2*j2*k2^2*i3^2 + 6*j2^2*k2^2*i3^2 + 5*i1^2*k2^3*i3^2 +
            6*i2*k2^3*i3^2 + 4*j2*k2^3*i3^2 + 6*k2^4*i3^2 + 2*i1^6*l2*i3^2 + i1^4*i2*l2*i3^2 + 2*i1^2*i2^2*l2*i3^2 +
            6*i1^4*j2*l2*i3^2 + 3*i1^2*i2*j2*l2*i3^2 + 2*i1^2*j2^2*l2*i3^2 + 5*i1^4*k2*l2*i3^2 + 2*i1^2*i2*k2*l2*i3^2 +
            3*i2^2*k2*l2*i3^2 + 6*i1^2*j2*k2*l2*i3^2 + i2*j2*k2*l2*i3^2 + 3*j2^2*k2*l2*i3^2 + 5*i1^2*k2^2*l2*i3^2 +
            3*i2*k2^2*l2*i3^2 + 5*j2*k2^2*l2*i3^2 + 2*k2^3*l2*i3^2 + 4*i1^4*l2^2*i3^2 + 2*i2^2*l2^2*i3^2 +
            3*i2*j2*l2^2*i3^2 + 2*j2^2*l2^2*i3^2 + i1^2*k2*l2^2*i3^2 + 5*i2*k2*l2^2*i3^2 + j2*k2*l2^2*i3^2 +
            2*k2^2*l2^2*i3^2 + 4*i1^2*l2^3*i3^2 + 4*i2*l2^3*i3^2 + 6*j2*l2^3*i3^2 + 4*l2^4*i3^2 + 2*i1^5*i3^3 +
            2*i1^3*k2*i3^3 + 4*i1*k2^2*i3^3 + 5*i1^3*l2*i3^3 + 5*i1*k2*l2*i3^3 + 5*i1*l2^2*i3^3 + 6*i1^9*i2*j3 +
            3*i1^5*i2^3*j3 + 4*i1^3*i2^4*j3 + 5*i1*i2^5*j3 + i1^9*j2*j3 + 3*i1^7*i2*j2*j3 + i1^5*i2^2*j2*j3 +
            i1^3*i2^3*j2*j3 + 3*i1*i2^4*j2*j3 + 5*i1^7*j2^2*j3 + 4*i1^5*i2*j2^2*j3 + i1^3*i2^2*j2^2*j3 + i1*i2^3*j2^2*j3 +
            6*i1^5*j2^3*j3 + 6*i1*i2^2*j2^3*j3 + i1^3*j2^4*j3 + 4*i1*i2*j2^4*j3 + 2*i1*j2^5*j3 + 4*i1^9*k2*j3 +
            4*i1^7*i2*k2*j3 + 2*i1^5*i2^2*k2*j3 + 6*i1^3*i2^3*k2*j3 + i1*i2^4*k2*j3 + 2*i1^7*j2*k2*j3 +
            2*i1^3*i2^2*j2*k2*j3 + 3*i1*i2^3*j2*k2*j3 + 3*i1^5*j2^2*k2*j3 + 6*i1^3*i2*j2^2*k2*j3 + 6*i1*i2^2*j2^2*k2*j3 +
            3*i1*i2*j2^3*k2*j3 + i1*j2^4*k2*j3 + 6*i1^7*k2^2*j3 + 4*i1^3*i2^2*k2^2*j3 + 6*i1*i2^3*k2^2*j3 +
            4*i1^5*j2*k2^2*j3 + 4*i1^3*i2*j2*k2^2*j3 + i1*i2^2*j2*k2^2*j3 + 5*i1*i2*j2^2*k2^2*j3 + 2*i1*j2^3*k2^2*j3 +
            i1*i2^2*k2^3*j3 + i1^3*j2*k2^3*j3 + 4*i1*i2*j2*k2^3*j3 + 3*i1^3*k2^4*j3 + 3*i1*i2*k2^4*j3 + 5*i1*j2*k2^4*j3 +
            2*i1*k2^5*j3 + 4*i1^9*l2*j3 + 4*i1^7*i2*l2*j3 + 4*i1^5*i2^2*l2*j3 + 2*i1^3*i2^3*l2*j3 + 4*i1*i2^4*l2*j3 +
            4*i1^7*j2*l2*j3 + 6*i1^5*i2*j2*l2*j3 + 2*i1^3*i2^2*j2*l2*j3 + 6*i1*i2^3*j2*l2*j3 + 4*i1^5*j2^2*l2*j3 +
            4*i1^3*i2*j2^2*l2*j3 + 6*i1^3*j2^3*l2*j3 + i1*i2*j2^3*l2*j3 + 3*i1*j2^4*l2*j3 + 4*i1^7*k2*l2*j3 +
            2*i1^5*i2*k2*l2*j3 + 2*i1^3*i2^2*k2*l2*j3 + 6*i1^5*j2*k2*l2*j3 + 3*i1^3*i2*j2*k2*l2*j3 + 4*i1*i2*j2^2*k2*l2*j3
            + 3*i1*j2^3*k2*l2*j3 + 6*i1^3*i2*k2^2*l2*j3 + 2*i1*i2^2*k2^2*l2*j3 + 6*i1^3*j2*k2^2*l2*j3 +
            6*i1*i2*j2*k2^2*l2*j3 + 5*i1*j2^2*k2^2*l2*j3 + 4*i1^3*k2^3*l2*j3 + 2*i1*i2*k2^3*l2*j3 + 6*i1*j2*k2^3*l2*j3 +
            5*i1*k2^4*l2*j3 + 4*i1^7*l2^2*j3 + i1^3*i2^2*l2^2*j3 + 4*i1*i2^3*l2^2*j3 + 6*i1*i2^2*j2*l2^2*j3 +
            6*i1^3*j2^2*l2^2*j3 + 2*i1*i2*j2^2*l2^2*j3 + 2*i1*j2^3*l2^2*j3 + 6*i1^5*k2*l2^2*j3 + 5*i1^3*i2*k2*l2^2*j3 +
            i1*i2^2*k2*l2^2*j3 + 4*i1^3*j2*k2*l2^2*j3 + 6*i1*i2*j2*k2*l2^2*j3 + i1*j2^2*k2*l2^2*j3 + 5*i1^3*k2^2*l2^2*j3 +
            4*i1*i2*k2^2*l2^2*j3 + 6*i1*j2*k2^2*l2^2*j3 + i1*k2^3*l2^2*j3 + i1^5*l2^3*j3 + 5*i1^3*i2*l2^3*j3 +
            2*i1^3*j2*l2^3*j3 + 3*i1*i2*j2*l2^3*j3 + 3*i1^3*k2*l2^3*j3 + 2*i1*i2*k2*l2^3*j3 + 3*i1*j2*k2*l2^3*j3 +
            2*i1*k2^2*l2^3*j3 + 3*i1^3*l2^4*j3 + i1*i2*l2^4*j3 + 2*i1*j2*l2^4*j3 + 6*i1*k2*l2^4*j3 + 6*i1*l2^5*j3 +
            4*i1^8*i3*j3 + i1^6*i2*i3*j3 + i1^2*i2^3*i3*j3 + 5*i1^6*j2*i3*j3 + 4*i1^4*i2*j2*i3*j3 + 4*i1^2*i2^2*j2*i3*j3 +
            3*i1^4*j2^2*i3*j3 + 3*i1^2*i2*j2^2*i3*j3 + 6*i1^2*j2^3*i3*j3 + 5*i1^4*i2*k2*i3*j3 + 4*i1^2*i2^2*k2*i3*j3 +
            4*i2^3*k2*i3*j3 + 2*i1^4*j2*k2*i3*j3 + i1^2*i2*j2*k2*i3*j3 + 2*i2^2*j2*k2*i3*j3 + 2*i1^2*j2^2*k2*i3*j3 +
            5*i2*j2^2*k2*i3*j3 + 3*j2^3*k2*i3*j3 + 4*i1^4*k2^2*i3*j3 + 2*i2^2*k2^2*i3*j3 + 6*i1^2*j2*k2^2*i3*j3 +
            5*i2*j2*k2^2*i3*j3 + 6*i1^2*k2^3*i3*j3 + 2*i2*k2^3*i3*j3 + 6*j2*k2^3*i3*j3 + 3*k2^4*i3*j3 + 5*i1^6*l2*i3*j3 +
            3*i1^4*i2*l2*i3*j3 + 6*i1^2*i2^2*l2*i3*j3 + 4*i2^3*l2*i3*j3 + 6*i1^4*j2*l2*i3*j3 + 3*i1^2*i2*j2*l2*i3*j3 +
            2*i2^2*j2*l2*i3*j3 + 5*i1^2*j2^2*l2*i3*j3 + 5*i2*j2^2*l2*i3*j3 + 3*j2^3*l2*i3*j3 + 5*i1^4*k2*l2*i3*j3 +
            2*i1^2*i2*k2*l2*i3*j3 + 3*i2^2*k2*l2*i3*j3 + 2*i1^2*j2*k2*l2*i3*j3 + 4*i2*j2*k2*l2*i3*j3 +
            2*i1^2*k2^2*l2*i3*j3 + 6*j2*k2^2*l2*i3*j3 + 4*k2^3*l2*i3*j3 + 2*i1^4*l2^2*i3*j3 + i1^2*i2*l2^2*i3*j3 +
            4*i1^2*j2*l2^2*i3*j3 + 3*i2*j2*l2^2*i3*j3 + 4*j2^2*l2^2*i3*j3 + 6*i1^2*k2*l2^2*i3*j3 + 4*i2*k2*l2^2*i3*j3 +
            6*j2*k2*l2^2*i3*j3 + 3*k2^2*l2^2*i3*j3 + 4*i1^2*l2^3*i3*j3 + 2*i2*l2^3*i3*j3 + 2*j2*l2^3*i3*j3 + 2*l2^4*i3*j3
            + i1^5*i3^2*j3 + i1^3*i2*i3^2*j3 + 6*i1^3*j2*i3^2*j3 + 5*i1*i2*k2*i3^2*j3 + 2*i1*j2*k2*i3^2*j3 +
            6*i1*k2^2*i3^2*j3 + i1^3*l2*i3^2*j3 + i1*i2*l2*i3^2*j3 + 6*i1*j2*l2*i3^2*j3 + 5*i1*k2*l2*i3^2*j3 +
            4*i1*l2^2*i3^2*j3 + 6*i1^8*j3^2 + 5*i1^6*i2*j3^2 + 2*i1^4*i2^2*j3^2 + 2*i1^2*i2^3*j3^2 + 3*i2^4*j3^2 +
            2*i1^4*i2*j2*j3^2 + 4*i1^2*i2^2*j2*j3^2 + 2*i2^3*j2*j3^2 + 3*i1^4*j2^2*j3^2 + 4*i2^2*j2^2*j3^2 +
            i1^2*j2^3*j3^2 + 2*i2*j2^3*j3^2 + 3*j2^4*j3^2 + 4*i1^4*i2*k2*j3^2 + 2*i1^2*i2^2*k2*j3^2 + 6*i2^3*k2*j3^2 +
            6*i1^4*j2*k2*j3^2 + 4*i1^2*i2*j2*k2*j3^2 + 5*i2^2*j2*k2*j3^2 + 2*i1^2*j2^2*k2*j3^2 + 3*j2^3*k2*j3^2 +
            2*i1^4*k2^2*j3^2 + 3*i1^2*i2*k2^2*j3^2 + 3*i2^2*k2^2*j3^2 + i1^2*j2*k2^2*j3^2 + 6*i2*j2*k2^2*j3^2 +
            4*j2^2*k2^2*j3^2 + i1^2*k2^3*j3^2 + 5*i2*k2^3*j3^2 + j2*k2^3*j3^2 + 4*i1^6*l2*j3^2 + 4*i1^4*i2*l2*j3^2 +
            3*i1^2*i2^2*l2*j3^2 + 6*i2^3*l2*j3^2 + i1^4*j2*l2*j3^2 + 3*i1^2*i2*j2*l2*j3^2 + 3*i2^2*j2*l2*j3^2 +
            3*i1^2*j2^2*l2*j3^2 + 4*i2*j2^2*l2*j3^2 + j2^3*l2*j3^2 + 3*i1^4*k2*l2*j3^2 + 2*i1^2*i2*k2*l2*j3^2 +
            3*i2^2*k2*l2*j3^2 + 3*i1^2*j2*k2*l2*j3^2 + 5*i1^2*k2^2*l2*j3^2 + 4*i2*k2^2*l2*j3^2 + 3*j2*k2^2*l2*j3^2 +
            5*k2^3*l2*j3^2 + 5*i1^4*l2^2*j3^2 + 3*i2^2*l2^2*j3^2 + 4*i1^2*j2*l2^2*j3^2 + 6*i2*j2*l2^2*j3^2 +
            j2^2*l2^2*j3^2 + 3*i1^2*k2*l2^2*j3^2 + 3*i2*k2*l2^2*j3^2 + j2*k2*l2^2*j3^2 + k2^2*l2^2*j3^2 + 3*i1^2*l2^3*j3^2
            + 5*i2*l2^3*j3^2 + 6*j2*l2^3*j3^2 + 5*k2*l2^3*j3^2 + l2^4*j3^2 + 5*i1^5*i3*j3^2 + i1^3*i2*i3*j3^2 +
            i1*i2^2*i3*j3^2 + 3*i1^3*j2*i3*j3^2 + 5*i1*i2*j2*i3*j3^2 + i1*j2^2*i3*j3^2 + 5*i1^3*k2*i3*j3^2 +
            6*i1*i2*k2*i3*j3^2 + 4*i1*j2*k2*i3*j3^2 + 4*i1*k2^2*i3*j3^2 + 3*i1^3*l2*i3*j3^2 + i1*i2*l2*i3*j3^2 +
            2*i1*j2*l2*i3*j3^2 + 5*i1*l2^2*i3*j3^2 + 4*i1^2*i3^2*j3^2 + 4*k2*i3^2*j3^2 + 6*l2*i3^2*j3^2 + 2*i1^5*j3^3 +
            2*i1^3*i2*j3^3 + 6*i1*i2^2*j3^3 + i1^3*j2*j3^3 + i1*i2*j2*j3^3 + i1^3*k2*j3^3 + 4*i1*i2*k2*j3^3 +
            6*i1*k2^2*j3^3 + 3*i1^3*l2*j3^3 + 6*i1*i2*l2*j3^3 + 2*i1*j2*l2*j3^3 + 3*i1*k2*l2*j3^3 + i1*l2^2*j3^3 +
            2*i1^2*i3*j3^3 + 2*i2*i3*j3^3 + 5*j2*i3*j3^3 + 6*k2*i3*j3^3 + 5*l2*i3*j3^3 + 3*j2*j3^4 + k2*j3^4 + 3*l2*j3^4 +
            5*i1^11*k3 + 4*i1^9*i2*k3 + 6*i1^7*i2^2*k3 + 6*i1^5*i2^3*k3 + 5*i1^3*i2^4*k3 + 3*i1*i2^5*k3 + 5*i1^9*j2*k3 +
            6*i1^5*i2^2*j2*k3 + 3*i1^3*i2^3*j2*k3 + 6*i1*i2^4*j2*k3 + 4*i1^7*j2^2*k3 + 6*i1^5*i2*j2^2*k3 +
            3*i1^3*i2^2*j2^2*k3 + 2*i1*i2^3*j2^2*k3 + 3*i1^5*j2^3*k3 + 5*i1*i2^2*j2^3*k3 + 3*i1^3*j2^4*k3 + i1*i2*j2^4*k3
            + 4*i1*j2^5*k3 + 5*i1^9*k2*k3 + 5*i1^7*i2*k2*k3 + 4*i1^5*i2^2*k2*k3 + 6*i1^3*i2^3*k2*k3 + i1*i2^4*k2*k3 +
            5*i1^7*j2*k2*k3 + 3*i1^3*i2^2*j2*k2*k3 + 6*i1*i2^3*j2*k2*k3 + i1^3*i2*j2^2*k2*k3 + 4*i1*i2^2*j2^2*k2*k3 +
            4*i1^3*j2^3*k2*k3 + 5*i1*i2*j2^3*k2*k3 + 5*i1*j2^4*k2*k3 + 6*i1^7*k2^2*k3 + i1^5*i2*k2^2*k3 +
            2*i1^3*i2^2*k2^2*k3 + 5*i1*i2^3*k2^2*k3 + 5*i1^5*j2*k2^2*k3 + 2*i1*i2^2*j2*k2^2*k3 + i1^3*j2^2*k2^2*k3 +
            3*i1*i2*j2^2*k2^2*k3 + 4*i1*j2^3*k2^2*k3 + 3*i1^5*k2^3*k3 + 5*i1*i2^2*k2^3*k3 + 3*i1^3*j2*k2^3*k3 +
            4*i1*j2^2*k2^3*k3 + 3*i1^3*k2^4*k3 + 5*i1*i2*k2^4*k3 + 4*i1*j2*k2^4*k3 + 4*i1^9*l2*k3 + 4*i1^7*i2*l2*k3 +
            i1^5*i2^2*l2*k3 + i1^3*i2^3*l2*k3 + 4*i1^7*j2*l2*k3 + 6*i1^3*i2^2*j2*l2*k3 + 3*i1^5*j2^2*l2*k3 +
            6*i1^3*i2*j2^2*l2*k3 + i1^3*j2^3*l2*k3 + 4*i1^7*k2*l2*k3 + 5*i1^5*i2*k2*l2*k3 + 6*i1*i2^3*k2*l2*k3 +
            4*i1^5*j2*k2*l2*k3 + 6*i1^3*i2*j2*k2*l2*k3 + 5*i1*i2^2*j2*k2*l2*k3 + 2*i1^3*j2^2*k2*l2*k3 +
            6*i1*i2*j2^2*k2*l2*k3 + 4*i1*j2^3*k2*l2*k3 + i1^5*k2^2*l2*k3 + 6*i1*i2^2*k2^2*l2*k3 + i1^3*j2*k2^2*l2*k3 +
            4*i1*i2*j2*k2^2*l2*k3 + 4*i1*j2^2*k2^2*l2*k3 + 5*i1^3*k2^3*l2*k3 + i1*i2*k2^3*l2*k3 + 5*i1*j2*k2^3*l2*k3 +
            6*i1*k2^4*l2*k3 + 4*i1^5*i2*l2^2*k3 + i1^3*i2^2*l2^2*k3 + 3*i1*i2^3*l2^2*k3 + i1^5*j2*l2^2*k3 +
            4*i1^3*i2*j2*l2^2*k3 + i1*i2^2*j2*l2^2*k3 + 3*i1^3*j2^2*l2^2*k3 + 3*i1*i2*j2^2*l2^2*k3 + 3*i1^5*k2*l2^2*k3 +
            4*i1^3*i2*k2*l2^2*k3 + 3*i1*i2^2*k2*l2^2*k3 + 3*i1^3*j2*k2*l2^2*k3 + 2*i1*j2^2*k2*l2^2*k3 +
            6*i1^3*k2^2*l2^2*k3 + 3*i1*i2*k2^2*l2^2*k3 + 2*i1*k2^3*l2^2*k3 + 6*i1^3*i2*l2^3*k3 + i1*i2^2*l2^3*k3 +
            i1^3*j2*l2^3*k3 + 4*i1*i2*j2*l2^3*k3 + 2*i1*j2^2*l2^3*k3 + 4*i1^3*k2*l2^3*k3 + 6*i1*j2*k2*l2^3*k3 +
            4*i1*k2^2*l2^3*k3 + 2*i1^3*l2^4*k3 + 6*i1*i2*l2^4*k3 + i1*j2*l2^4*k3 + 3*i1*k2*l2^4*k3 + 2*i1^8*j3*k3 +
            3*i1^6*i2*j3*k3 + 2*i1^4*i2^2*j3*k3 + 4*i1^2*i2^3*j3*k3 + 6*i2^4*j3*k3 + 5*i1^6*j2*j3*k3 + 6*i1^4*i2*j2*j3*k3
            + 2*i1^2*i2^2*j2*j3*k3 + 4*i2^3*j2*j3*k3 + 5*i1^2*i2*j2^2*j3*k3 + i2^2*j2^2*j3*k3 + 3*i1^2*j2^3*j3*k3 +
            4*i2*j2^3*j3*k3 + 6*j2^4*j3*k3 + 5*i1^6*k2*j3*k3 + 3*i1^4*i2*k2*j3*k3 + 3*i1^2*i2^2*k2*j3*k3 + 4*i2^3*k2*j3*k3
            + 2*i1^4*j2*k2*j3*k3 + 4*i1^2*i2*j2*k2*j3*k3 + i2^2*j2*k2*j3*k3 + 5*i1^2*j2^2*k2*j3*k3 + 2*j2^3*k2*j3*k3 +
            2*i1^4*k2^2*j3*k3 + 6*i1^2*i2*k2^2*j3*k3 + 2*i2^2*k2^2*j3*k3 + 5*i1^2*j2*k2^2*j3*k3 + i2*j2*k2^2*j3*k3 +
            6*j2^2*k2^2*j3*k3 + 4*i1^2*k2^3*j3*k3 + 6*i2*k2^3*j3*k3 + 5*j2*k2^3*j3*k3 + k2^4*j3*k3 + 4*i1^4*i2*l2*j3*k3 +
            4*i1^2*i2^2*l2*j3*k3 + 3*i2^3*l2*j3*k3 + 3*i1^4*j2*l2*j3*k3 + 6*i1^2*i2*j2*l2*j3*k3 + 5*i2^2*j2*l2*j3*k3 +
            4*i1^2*j2^2*l2*j3*k3 + 2*i2*j2^2*l2*j3*k3 + 4*j2^3*l2*j3*k3 + i1^4*k2*l2*j3*k3 + 4*i1^2*i2*k2*l2*j3*k3 +
            5*i2^2*k2*l2*j3*k3 + 3*i1^2*j2*k2*l2*j3*k3 + 3*i2*j2*k2*l2*j3*k3 + 3*j2^2*k2*l2*j3*k3 + 4*i1^2*k2^2*l2*j3*k3 +
            i2*k2^2*l2*j3*k3 + 4*j2*k2^2*l2*j3*k3 + k2^3*l2*j3*k3 + 5*i1^4*l2^2*j3*k3 + 2*i2^2*l2^2*j3*k3 +
            i1^2*j2*l2^2*j3*k3 + 6*i2*j2*l2^2*j3*k3 + 2*j2^2*l2^2*j3*k3 + 5*i1^2*k2*l2^2*j3*k3 + 4*i2*k2*l2^2*j3*k3 +
            3*j2*k2*l2^2*j3*k3 + 2*k2^2*l2^2*j3*k3 + 6*i2*l2^3*j3*k3 + 2*k2*l2^3*j3*k3 + 2*l2^4*j3*k3 + 4*i1^3*i2*j3^2*k3
            + 5*i1^3*j2*j3^2*k3 + i1*i2*j2*j3^2*k3 + 6*i1*j2^2*j3^2*k3 + 2*i1*i2*k2*j3^2*k3 + 3*i1*k2^2*j3^2*k3 +
            4*i1^3*l2*j3^2*k3 + 2*i1*i2*l2*j3^2*k3 + 3*i1*j2*l2*j3^2*k3 + i1*k2*l2*j3^2*k3 + 4*i1*l2^2*j3^2*k3 +
            4*i1^2*j3^3*k3 + 3*i2*j3^3*k3 + 5*j2*j3^3*k3 + 6*k2*j3^3*k3 + 3*i1^4*i2^2*k3^2 + 2*i1^2*i2^3*k3^2 +
            4*i2^4*k3^2 + 3*i1^6*j2*k3^2 + 6*i1^2*i2^2*j2*k3^2 + 5*i2^3*j2*k3^2 + 4*i1^4*j2^2*k3^2 + 3*i1^2*i2*j2^2*k3^2 +
            3*i2^2*j2^2*k3^2 + 3*i1^2*j2^3*k3^2 + 5*i2*j2^3*k3^2 + 4*j2^4*k3^2 + 2*i1^4*i2*k2*k3^2 + 4*i2^3*k2*k3^2 +
            3*i1^4*j2*k2*k3^2 + 2*i1^2*i2*j2*k2*k3^2 + 3*i2^2*j2*k2*k3^2 + i1^2*j2^2*k2*k3^2 + 3*i2*j2^2*k2*k3^2 +
            4*j2^3*k2*k3^2 + 4*i1^4*k2^2*k3^2 + 3*i1^2*i2*k2^2*k3^2 + 4*i2^2*k2^2*k3^2 + 4*i1^2*j2*k2^2*k3^2 +
            2*i2*j2*k2^2*k3^2 + 4*j2^2*k2^2*k3^2 + 5*i1^2*k2^3*k3^2 + 4*i2*k2^3*k3^2 + 3*k2^4*k3^2 + 2*i1^6*l2*k3^2 +
            2*i1^2*i2^2*l2*k3^2 + i2^3*l2*k3^2 + 4*i1^4*j2*l2*k3^2 + 5*i1^2*i2*j2*l2*k3^2 + 2*i2^2*j2*l2*k3^2 +
            4*i1^2*j2^2*l2*k3^2 + 4*j2^3*l2*k3^2 + 6*i1^4*k2*l2*k3^2 + 5*i1^2*i2*k2*l2*k3^2 + i2^2*k2*l2*k3^2 +
            4*i1^2*j2*k2*l2*k3^2 + 3*i2*j2*k2*l2*k3^2 + 6*j2^2*k2*l2*k3^2 + 6*i1^2*k2^2*l2*k3^2 + i2*k2^2*l2*k3^2 +
            j2*k2^2*l2*k3^2 + i1^4*l2^2*k3^2 + 3*i1^2*i2*l2^2*k3^2 + 5*i2^2*l2^2*k3^2 + 5*i1^2*j2*l2^2*k3^2 +
            4*i2*j2*l2^2*k3^2 + j2^2*l2^2*k3^2 + 2*i2*k2*l2^2*k3^2 + 5*k2^2*l2^2*k3^2 + 2*i1^2*l2^3*k3^2 + i2*l2^3*k3^2 +
            6*j2*l2^3*k3^2 + 5*k2*l2^3*k3^2 + 6*l2^4*k3^2 + 2*i1^3*i2*j3*k3^2 + i1*i2^2*j3*k3^2 + 3*i1^3*j2*j3*k3^2 +
            4*i1*i2*j2*j3*k3^2 + 2*i1*j2^2*j3*k3^2 + 4*i1^3*k2*j3*k3^2 + 6*i1*i2*k2*j3*k3^2 + 6*i1*j2*k2*j3*k3^2 +
            4*i1*k2^2*j3*k3^2 + i1^3*l2*j3*k3^2 + 6*i1*i2*l2*j3*k3^2 + i1*j2*l2*j3*k3^2 + 3*i1*k2*l2*j3*k3^2 +
            4*i1*l2^2*j3*k3^2 + 4*i2*j3^2*k3^2 + 6*j2*j3^2*k3^2 + 3*k2*j3^2*k3^2 + l2*j3^2*k3^2 + 5*i1^5*k3^3 +
            5*i1^3*i2*k3^3 + 4*i1*i2^2*k3^3 + 5*i1^3*j2*k3^3 + i1*i2*j2*k3^3 + 2*i1*j2^2*k3^3 + 3*i1^3*k2*k3^3 +
            6*i1*i2*k2*k3^3 + 5*i1*k2^2*k3^3 + 3*i1^3*l2*k3^3 + 3*i1*i2*l2*k3^3 + 6*i1*j2*l2*k3^3 + 3*i1*k2*l2*k3^3 +
            3*i1*l2^2*k3^3 + i1^2*j3*k3^3 + 2*j2*j3*k3^3 + k2*j3*k3^3 + 6*l2*j3*k3^3 + 2*i1^2*k3^4 + 3*i2*k3^4 + 3*j2*k3^4
            + 4*l2*k3^4 + 3*i1^11*l3 + 3*i1^9*i2*l3 + 4*i1^7*i2^2*l3 + i1^5*i2^3*l3 + 2*i1^3*i2^4*l3 + 6*i1*i2^5*l3 +
            2*i1^9*j2*l3 + 6*i1^7*i2*j2*l3 + 2*i1^5*i2^2*j2*l3 + 2*i1^3*i2^3*j2*l3 + 3*i1*i2^4*j2*l3 + 5*i1^7*j2^2*l3 +
            3*i1^5*i2*j2^2*l3 + 3*i1^3*i2^2*j2^2*l3 + 5*i1*i2^3*j2^2*l3 + 3*i1^5*j2^3*l3 + i1^3*i2*j2^3*l3 +
            5*i1*i2^2*j2^3*l3 + 6*i1^3*j2^4*l3 + 3*i1*i2*j2^4*l3 + 6*i1*j2^5*l3 + 3*i1^7*i2*k2*l3 + 5*i1^5*i2^2*k2*l3 +
            5*i1^3*i2^3*k2*l3 + 4*i1*i2^4*k2*l3 + 2*i1^7*j2*k2*l3 + i1^5*i2*j2*k2*l3 + 2*i1^3*i2^2*j2*k2*l3 +
            5*i1*i2^3*j2*k2*l3 + i1^3*i2*j2^2*k2*l3 + 4*i1*i2^2*j2^2*k2*l3 + 3*i1^3*j2^3*k2*l3 + 3*i1*i2*j2^3*k2*l3 +
            5*i1*j2^4*k2*l3 + 2*i1^7*k2^2*l3 + i1^5*i2*k2^2*l3 + i1^3*i2^2*k2^2*l3 + 5*i1*i2^3*k2^2*l3 + 6*i1^5*j2*k2^2*l3
            + 3*i1^3*j2^2*k2^2*l3 + 6*i1*i2*j2^2*k2^2*l3 + i1*j2^3*k2^2*l3 + 4*i1^5*k2^3*l3 + 2*i1*i2^2*k2^3*l3 +
            2*i1^3*j2*k2^3*l3 + 2*i1*i2*j2*k2^3*l3 + 6*i1*j2^2*k2^3*l3 + 2*i1^3*k2^4*l3 + 3*i1*i2*k2^4*l3 + i1*j2*k2^4*l3
            + 4*i1^9*l2*l3 + i1^7*i2*l2*l3 + 6*i1^3*i2^3*l2*l3 + i1*i2^4*l2*l3 + 4*i1^7*j2*l2*l3 + 4*i1^5*i2*j2*l2*l3 +
            3*i1^3*i2^2*j2*l2*l3 + 2*i1*i2^3*j2*l2*l3 + 2*i1^5*j2^2*l2*l3 + 4*i1^3*i2*j2^2*l2*l3 + 4*i1*i2^2*j2^2*l2*l3 +
            5*i1^3*j2^3*l2*l3 + 3*i1*i2*j2^3*l2*l3 + 4*i1*j2^4*l2*l3 + 5*i1^7*k2*l2*l3 + 4*i1^5*i2*k2*l2*l3 +
            2*i1^3*i2^2*k2*l2*l3 + 5*i1^3*i2*j2*k2*l2*l3 + 2*i1*i2^2*j2*k2*l2*l3 + 2*i1^3*j2^2*k2*l2*l3 +
            5*i1*i2*j2^2*k2*l2*l3 + i1*j2^3*k2*l2*l3 + 6*i1^3*i2*k2^2*l2*l3 + 4*i1*i2^2*k2^2*l2*l3 + 3*i1^3*j2*k2^2*l2*l3
            + 3*i1*i2*j2*k2^2*l2*l3 + 3*i1*j2^2*k2^2*l2*l3 + 5*i1^3*k2^3*l2*l3 + 6*i1*i2*k2^3*l2*l3 + 2*i1*j2*k2^3*l2*l3 +
            i1*k2^4*l2*l3 + 6*i1^7*l2^2*l3 + 2*i1^3*i2^2*l2^2*l3 + 2*i1*i2^3*l2^2*l3 + 6*i1^5*j2*l2^2*l3 +
            3*i1^3*i2*j2*l2^2*l3 + 6*i1*i2^2*j2*l2^2*l3 + i1^3*j2^2*l2^2*l3 + i1*i2*j2^2*l2^2*l3 + 4*i1*j2^3*l2^2*l3 +
            2*i1^5*k2*l2^2*l3 + i1^3*i2*k2*l2^2*l3 + 5*i1*i2^2*k2*l2^2*l3 + 2*i1*i2*j2*k2*l2^2*l3 + 2*i1*j2^2*k2*l2^2*l3 +
            6*i1^3*k2^2*l2^2*l3 + 3*i1*j2*k2^2*l2^2*l3 + i1*k2^3*l2^2*l3 + 4*i1^5*l2^3*l3 + 3*i1*i2^2*l2^3*l3 +
            i1*i2*j2*l2^3*l3 + 2*i1*j2^2*l2^3*l3 + 6*i1^3*k2*l2^3*l3 + 6*i1*i2*k2*l2^3*l3 + 4*i1*j2*k2*l2^3*l3 +
            6*i1*k2^2*l2^3*l3 + 2*i1^3*l2^4*l3 + 2*i1*i2*l2^4*l3 + 2*i1*j2*l2^4*l3 + 4*i1*k2*l2^4*l3 + 3*i1*l2^5*l3 +
            6*i1^8*i3*l3 + 3*i1^6*i2*i3*l3 + i1^4*i2^2*i3*l3 + 2*i1^2*i2^3*i3*l3 + 3*i2^4*i3*l3 + i1^6*j2*i3*l3 +
            2*i2^3*j2*i3*l3 + i1^2*i2*j2^2*i3*l3 + 4*i2^2*j2^2*i3*l3 + 4*i1^2*j2^3*i3*l3 + 2*i2*j2^3*i3*l3 + 3*j2^4*i3*l3
            + 6*i1^6*k2*i3*l3 + 3*i1^4*i2*k2*i3*l3 + 6*i1^2*i2^2*k2*i3*l3 + 4*i2^3*k2*i3*l3 + 4*i1^4*j2*k2*i3*l3 +
            4*i1^2*i2*j2*k2*i3*l3 + 6*i1^2*j2^2*k2*i3*l3 + 2*i2*j2^2*k2*i3*l3 + j2^3*k2*i3*l3 + 4*i1^4*k2^2*i3*l3 +
            6*i2^2*k2^2*i3*l3 + 2*i1^2*j2*k2^2*i3*l3 + 6*i2*j2*k2^2*i3*l3 + 4*j2^2*k2^2*i3*l3 + i1^2*k2^3*i3*l3 +
            4*i2*k2^3*i3*l3 + 2*j2*k2^3*i3*l3 + k2^4*i3*l3 + 4*i1^4*i2*l2*i3*l3 + 2*i1^2*i2^2*l2*i3*l3 + i2^3*l2*i3*l3 +
            2*i1^4*j2*l2*i3*l3 + 2*i2^2*j2*l2*i3*l3 + 3*i1^2*j2^2*l2*i3*l3 + 4*j2^3*l2*i3*l3 + i1^4*k2*l2*i3*l3 +
            3*i1^2*i2*k2*l2*i3*l3 + 3*i2^2*k2*l2*i3*l3 + 5*i1^2*j2*k2*l2*i3*l3 + 3*i2*j2*k2*l2*i3*l3 + 6*j2^2*k2*l2*i3*l3
            + 5*i2*k2^2*l2*i3*l3 + 5*k2^3*l2*i3*l3 + 4*i2^2*l2^2*i3*l3 + 2*i1^2*j2*l2^2*i3*l3 + 6*i2*j2*l2^2*i3*l3 +
            2*j2^2*l2^2*i3*l3 + 2*i1^2*k2*l2^2*i3*l3 + 6*i2*k2*l2^2*i3*l3 + j2*k2*l2^2*i3*l3 + 6*k2^2*l2^2*i3*l3 +
            3*i1^2*l2^3*i3*l3 + 3*i2*l2^3*i3*l3 + 4*j2*l2^3*i3*l3 + 4*l2^4*i3*l3 + 2*i1^5*i3^2*l3 + 2*i1*i2^2*i3^2*l3 +
            5*i1^3*j2*i3^2*l3 + 3*i1*i2*j2*i3^2*l3 + 2*i1*j2^2*i3^2*l3 + 2*i1^3*k2*i3^2*l3 + 3*i1*i2*k2*i3^2*l3 +
            5*i1*j2*k2*i3^2*l3 + 3*i1*k2^2*i3^2*l3 + 4*i1^3*l2*i3^2*l3 + 6*i1*j2*l2*i3^2*l3 + i1*k2*l2*i3^2*l3 +
            5*i1*l2^2*i3^2*l3 + 6*i1^2*i3^3*l3 + k2*i3^3*l3 + l2*i3^3*l3 + 6*i1^8*j3*l3 + 5*i1^6*i2*j3*l3 +
            6*i1^4*i2^2*j3*l3 + 2*i1^2*i2^3*j3*l3 + 6*i2^4*j3*l3 + 4*i1^6*j2*j3*l3 + 6*i1^4*i2*j2*j3*l3 +
            3*i1^2*i2^2*j2*j3*l3 + i2^3*j2*j3*l3 + 2*i1^4*j2^2*j3*l3 + 6*i1^2*i2*j2^2*j3*l3 + 3*i2^2*j2^2*j3*l3 +
            3*i1^2*j2^3*j3*l3 + 2*i2*j2^3*j3*l3 + 2*j2^4*j3*l3 + i1^6*k2*j3*l3 + 4*i1^4*i2*k2*j3*l3 + i1^2*i2^2*k2*j3*l3 +
            2*i2^3*k2*j3*l3 + 3*i1^4*j2*k2*j3*l3 + i1^2*i2*j2*k2*j3*l3 + 3*i2^2*j2*k2*j3*l3 + 6*i1^2*j2^2*k2*j3*l3 +
            6*i2*j2^2*k2*j3*l3 + 3*j2^3*k2*j3*l3 + i1^4*k2^2*j3*l3 + 5*i2^2*k2^2*j3*l3 + 2*i2*j2*k2^2*j3*l3 +
            2*j2^2*k2^2*j3*l3 + i1^2*k2^3*j3*l3 + 6*i2*k2^3*j3*l3 + j2*k2^3*j3*l3 + 5*k2^4*j3*l3 + 3*i1^6*l2*j3*l3 +
            6*i1^4*i2*l2*j3*l3 + 6*i1^2*i2^2*l2*j3*l3 + 2*i2^3*l2*j3*l3 + 2*i1^4*j2*l2*j3*l3 + 2*i1^2*i2*j2*l2*j3*l3 +
            2*i2^2*j2*l2*j3*l3 + i1^2*j2^2*l2*j3*l3 + 4*i2*j2^2*l2*j3*l3 + 6*j2^3*l2*j3*l3 + 4*i2^2*k2*l2*j3*l3 +
            4*i1^2*j2*k2*l2*j3*l3 + 6*i2*j2*k2*l2*j3*l3 + 6*j2^2*k2*l2*j3*l3 + 4*i1^2*k2^2*l2*j3*l3 + 5*j2*k2^2*l2*j3*l3 +
            i1^2*i2*l2^2*j3*l3 + 2*i2^2*l2^2*j3*l3 + 4*i1^2*j2*l2^2*j3*l3 + 5*i2*j2*l2^2*j3*l3 + 4*i1^2*k2*l2^2*j3*l3 +
            3*i2*k2*l2^2*j3*l3 + 5*j2*k2*l2^2*j3*l3 + 2*k2^2*l2^2*j3*l3 + 3*i2*l2^3*j3*l3 + 2*j2*l2^3*j3*l3 +
            k2*l2^3*j3*l3 + 3*l2^4*j3*l3 + 4*i1^5*j3^2*l3 + 3*i1^3*i2*j3^2*l3 + 5*i1*i2^2*j3^2*l3 + 2*i1^3*j2*j3^2*l3 +
            3*i1*i2*j2*j3^2*l3 + 4*i1*j2^2*j3^2*l3 + i1^3*k2*j3^2*l3 + 6*i1*i2*k2*j3^2*l3 + 5*i1*j2*k2*j3^2*l3 +
            4*i1*k2^2*j3^2*l3 + 4*i1^3*l2*j3^2*l3 + 5*i1*i2*l2*j3^2*l3 + 2*i1*k2*l2*j3^2*l3 + 4*i2*j3^3*l3 + 2*i1^8*l3^2 +
            5*i1^6*i2*l3^2 + 6*i1^4*i2^2*l3^2 + 6*i1^6*j2*l3^2 + i1^4*i2*j2*l3^2 + i1^2*i2^2*j2*l3^2 + 5*i2^3*j2*l3^2 +
            i1^4*j2^2*l3^2 + i1^2*i2*j2^2*l3^2 + 3*i1^2*j2^3*l3^2 + 6*i2*j2^3*l3^2 + 3*j2^4*l3^2 + 4*i1^6*k2*l3^2 +
            i1^4*i2*k2*l3^2 + 2*i1^2*i2^2*k2*l3^2 + 5*i2^3*k2*l3^2 + 5*i1^2*i2*j2*k2*l3^2 + 4*i2^2*j2*k2*l3^2 +
            6*i1^2*j2^2*k2*l3^2 + 4*i2*j2^2*k2*l3^2 + 4*j2^3*k2*l3^2 + i1^4*k2^2*l3^2 + 6*i1^2*i2*k2^2*l3^2 +
            5*i2*j2*k2^2*l3^2 + 2*j2^2*k2^2*l3^2 + 3*i2*k2^3*l3^2 + 2*k2^4*l3^2 + 6*i1^6*l2*l3^2 + 5*i1^4*i2*l2*l3^2 +
            6*i1^2*i2^2*l2*l3^2 + 5*i2^3*l2*l3^2 + i1^4*j2*l2*l3^2 + 2*i2^2*j2*l2*l3^2 + 6*i1^2*j2^2*l2*l3^2 +
            6*i2*j2^2*l2*l3^2 + 3*j2^3*l2*l3^2 + 3*i1^4*k2*l2*l3^2 + 4*i1^2*i2*k2*l2*l3^2 + 6*i2^2*k2*l2*l3^2 +
            i2*j2*k2*l2*l3^2 + 4*j2^2*k2*l2*l3^2 + 5*i1^2*k2^2*l2*l3^2 + 6*i2*k2^2*l2*l3^2 + 4*i1^4*l2^2*l3^2 +
            3*i1^2*i2*l2^2*l3^2 + i2^2*l2^2*l3^2 + 2*i1^2*j2*l2^2*l3^2 + 3*i2*j2*l2^2*l3^2 + 3*i1^2*k2*l2^2*l3^2 +
            5*i2*k2*l2^2*l3^2 + 5*j2*k2*l2^2*l3^2 + 2*k2^2*l2^2*l3^2 + 2*i1^2*l2^3*l3^2 + 3*i2*l2^3*l3^2 + 4*j2*l2^3*l3^2
            + 3*k2*l2^3*l3^2 + 6*l2^4*l3^2 + 3*i1^5*i3*l3^2 + 6*i1^3*i2*i3*l3^2 + i1*i2^2*i3*l3^2 + i1^3*j2*i3*l3^2 +
            4*i1*i2*j2*i3*l3^2 + 4*i1*i2*k2*i3*l3^2 + 4*i1*j2*k2*i3*l3^2 + 3*i1*k2^2*i3*l3^2 + 4*i1^3*l2*i3*l3^2 +
            3*i1*i2*l2*i3*l3^2 + 5*i1*j2*l2*i3*l3^2 + 4*i1*k2*l2*i3*l3^2 + 2*i1*l2^2*i3*l3^2 + i1^2*i3^2*l3^2 +
            j2*i3^2*l3^2 + 6*l2*i3^2*l3^2 + 2*i1^3*i2*j3*l3^2 + 2*i1*i2^2*j3*l3^2 + 4*i1^3*j2*j3*l3^2 + 2*i1*i2*j2*j3*l3^2
            + 4*i1*j2^2*j3*l3^2 + 5*i1^3*k2*j3*l3^2 + 6*i1*i2*k2*j3*l3^2 + i1*j2*k2*j3*l3^2 + 2*i1*k2^2*j3*l3^2 +
            6*i1*i2*l2*j3*l3^2 + 2*i1*j2*l2*j3*l3^2 + 3*i1*k2*l2*j3*l3^2 + 4*i1*l2^2*j3*l3^2 + 4*i1^2*j3^2*l3^2 +
            6*i2*j3^2*l3^2 + 4*j2*j3^2*l3^2 + 6*k2*j3^2*l3^2 + 4*l2*j3^2*l3^2 + 5*i1^5*l3^3 + i1^3*i2*l3^3 +
            4*i1*i2^2*l3^3 + 3*i1^3*j2*l3^3 + 2*i1*i2*j2*l3^3 + 6*i1*j2^2*l3^3 + 2*i1*i2*k2*l3^3 + 6*i1*j2*k2*l3^3 +
            4*i1^3*l2*l3^3 + 5*i1*i2*l2*l3^3 + 2*i1*j2*l2*l3^3 + 2*i1*k2*l2*l3^3 + 4*i1*l2^2*l3^3 + 3*i1^2*i3*l3^3 +
            6*i2*i3*l3^3 + 4*k2*i3*l3^3 + l2*i3*l3^3 + 2*i1^2*j3*l3^3 + 5*i2*j3*l3^3 + k2*j3*l3^3 + 3*l2*j3*l3^3 +
            i1^2*l3^4 + 6*i2*l3^4 + 4*j2*l3^4 + 2*k2*l3^4 + 6*l2*l3^4 + i1^10*i4 + 5*i1^8*i2*i4 + 5*i1^6*i2^2*i4 +
            3*i1^4*i2^3*i4 + 2*i2^5*i4 + 5*i1^8*j2*i4 + 3*i1^6*i2*j2*i4 + 2*i1^4*i2^2*j2*i4 + 2*i1^2*i2^3*j2*i4 +
            4*i2^4*j2*i4 + i1^6*j2^2*i4 + i1^2*i2^2*j2^2*i4 + 6*i2^3*j2^2*i4 + 2*i1^4*j2^3*i4 + 6*i1^2*i2*j2^3*i4 +
            i2^2*j2^3*i4 + 5*i1^2*j2^4*i4 + 3*i2*j2^4*i4 + 5*j2^5*i4 + 3*i1^8*k2*i4 + 2*i1^6*i2*k2*i4 + 3*i1^4*i2^2*k2*i4
            + i2^4*k2*i4 + 6*i1^6*j2*k2*i4 + 6*i1^4*i2*j2*k2*i4 + 4*i1^2*i2^2*j2*k2*i4 + 2*i2^3*j2*k2*i4 +
            4*i1^4*j2^2*k2*i4 + 5*i1^2*i2*j2^2*k2*i4 + 2*i2^2*j2^2*k2*i4 + 5*i1^2*j2^3*k2*i4 + 2*j2^4*k2*i4 +
            4*i1^6*k2^2*i4 + 6*i1^4*i2*k2^2*i4 + i1^2*i2^2*k2^2*i4 + 4*i2^3*k2^2*i4 + i1^4*j2*k2^2*i4 +
            5*i1^2*i2*j2*k2^2*i4 + 5*i2^2*j2*k2^2*i4 + 4*i1^2*j2^2*k2^2*i4 + 2*i2*j2^2*k2^2*i4 + 3*j2^3*k2^2*i4 +
            2*i1^4*k2^3*i4 + 6*i1^2*i2*k2^3*i4 + 4*i2^2*k2^3*i4 + 6*i1^2*j2*k2^3*i4 + 2*i2*j2*k2^3*i4 + 4*j2^2*k2^3*i4 +
            4*i1^2*k2^4*i4 + 2*i2*k2^4*i4 + 4*j2*k2^4*i4 + k2^5*i4 + 5*i1^8*l2*i4 + 4*i1^6*i2*l2*i4 + 4*i1^4*i2^2*l2*i4 +
            6*i1^2*i2^3*l2*i4 + i2^4*l2*i4 + 5*i1^6*j2*l2*i4 + 6*i1^4*i2*j2*l2*i4 + 4*i1^2*i2^2*j2*l2*i4 + i2^3*j2*l2*i4 +
            3*i1^2*i2*j2^2*l2*i4 + 5*i2^2*j2^2*l2*i4 + i1^2*j2^3*l2*i4 + 4*i2*j2^3*l2*i4 + 3*j2^4*l2*i4 + 2*i1^6*k2*l2*i4
            + 4*i1^4*i2*k2*l2*i4 + 4*i1^2*i2^2*k2*l2*i4 + 5*i2^3*k2*l2*i4 + 6*i1^4*j2*k2*l2*i4 + 6*i1^2*i2*j2*k2*l2*i4 +
            6*i2^2*j2*k2*l2*i4 + 4*i1^2*j2^2*k2*l2*i4 + i2*j2^2*k2*l2*i4 + 2*j2^3*k2*l2*i4 + 2*i1^2*i2*k2^2*l2*i4 +
            6*i2^2*k2^2*l2*i4 + 5*i2*j2*k2^2*l2*i4 + 5*j2^2*k2^2*l2*i4 + 2*i2*k2^3*l2*i4 + 5*j2*k2^3*l2*i4 + 5*k2^4*l2*i4
            + i1^6*l2^2*i4 + 2*i1^4*i2*l2^2*i4 + i1^2*i2^2*l2^2*i4 + 5*i2^3*l2^2*i4 + 2*i1^4*j2*l2^2*i4 +
            4*i1^2*i2*j2*l2^2*i4 + i2^2*j2*l2^2*i4 + i1^2*j2^2*l2^2*i4 + j2^3*l2^2*i4 + 4*i1^4*k2*l2^2*i4 +
            6*i1^2*i2*k2*l2^2*i4 + 2*i2^2*k2*l2^2*i4 + 3*i1^2*j2*k2*l2^2*i4 + 5*j2^2*k2*l2^2*i4 + 6*i1^2*k2^2*l2^2*i4 +
            5*i2*k2^2*l2^2*i4 + 2*j2*k2^2*l2^2*i4 + 3*i1^4*l2^3*i4 + 4*i2^2*l2^3*i4 + 5*i1^2*j2*l2^3*i4 + 4*j2^2*l2^3*i4 +
            3*i1^2*k2*l2^3*i4 + 4*i2*k2*l2^3*i4 + 3*j2*k2*l2^3*i4 + 2*k2^2*l2^3*i4 + 4*i1^2*l2^4*i4 + 3*i2*l2^4*i4 +
            3*k2*l2^4*i4 + 5*l2^5*i4 + 4*i1^7*i3*i4 + 3*i1^5*i2*i3*i4 + 4*i1^3*i2^2*i3*i4 + 4*i1*i2^3*i3*i4 +
            3*i1^5*j2*i3*i4 + i1^3*i2*j2*i3*i4 + 2*i1*i2^2*j2*i3*i4 + 2*i1^3*j2^2*i3*i4 + 5*i1*i2*j2^2*i3*i4 +
            3*i1*j2^3*i3*i4 + 2*i1^5*k2*i3*i4 + 5*i1^3*i2*k2*i3*i4 + i1*i2^2*k2*i3*i4 + 2*i1*i2*j2*k2*i3*i4 +
            4*i1*j2^2*k2*i3*i4 + 6*i1^3*k2^2*i3*i4 + 5*i1*i2*k2^2*i3*i4 + 5*i1*j2*k2^2*i3*i4 + 2*i1*k2^3*i3*i4 +
            i1^5*l2*i3*i4 + 2*i1^3*i2*l2*i3*i4 + 4*i1*i2^2*l2*i3*i4 + 5*i1^3*j2*l2*i3*i4 + 3*i1*j2^2*l2*i3*i4 +
            3*i1^3*k2*l2*i3*i4 + 6*i1*i2*k2*l2*i3*i4 + 2*i1*j2*k2*l2*i3*i4 + 5*i1*k2^2*l2*i3*i4 + 3*i1^3*l2^2*i3*i4 +
            5*i1*i2*l2^2*i3*i4 + 5*i1*j2*l2^2*i3*i4 + i1*k2*l2^2*i3*i4 + 6*i1*l2^3*i3*i4 + 2*i1^4*i3^2*i4 +
            3*i1^2*i2*i3^2*i4 + 4*i1^2*j2*i3^2*i4 + 6*i1^2*k2*i3^2*i4 + 4*i2*k2*i3^2*i4 + 3*j2*k2*i3^2*i4 + 5*k2^2*i3^2*i4
            + 6*i1^2*l2*i3^2*i4 + 2*i2*l2*i3^2*i4 + 5*j2*l2*i3^2*i4 + k2*l2*i3^2*i4 + 6*l2^2*i3^2*i4 + 2*i1^7*j3*i4 +
            5*i1^5*i2*j3*i4 + 2*i1^3*i2^2*j3*i4 + 5*i1*i2^3*j3*i4 + 6*i1^5*j2*j3*i4 + 5*i1^3*i2*j2*j3*i4 +
            2*i1*i2^2*j2*j3*i4 + 4*i1^3*j2^2*j3*i4 + 2*i1*i2*j2^2*j3*i4 + 5*i1*j2^3*j3*i4 + 2*i1^5*k2*j3*i4 +
            2*i1^3*i2*k2*j3*i4 + i1^3*j2*k2*j3*i4 + 4*i1*i2*j2*k2*j3*i4 + 2*i1*j2^2*k2*j3*i4 + i1^3*k2^2*j3*i4 +
            i1*j2*k2^2*j3*i4 + i1^5*l2*j3*i4 + 4*i1*i2*j2*l2*j3*i4 + 4*i1*j2^2*l2*j3*i4 + 3*i1^3*k2*l2*j3*i4 +
            5*i1*i2*k2*l2*j3*i4 + i1*j2*k2*l2*j3*i4 + 4*i1*k2^2*l2*j3*i4 + 6*i1^3*l2^2*j3*i4 + i1*i2*l2^2*j3*i4 +
            3*i1*j2*l2^2*j3*i4 + 2*i1*k2*l2^2*j3*i4 + 2*i1*l2^3*j3*i4 + 6*i1^4*i3*j3*i4 + 2*i1^2*i2*i3*j3*i4 +
            4*i2^2*i3*j3*i4 + 2*i1^2*j2*i3*j3*i4 + 6*i2*j2*i3*j3*i4 + 4*j2^2*i3*j3*i4 + 5*i1^2*k2*i3*j3*i4 +
            2*j2*k2*i3*j3*i4 + 2*k2^2*i3*j3*i4 + 4*i1^2*l2*i3*j3*i4 + 6*i2*l2*i3*j3*i4 + 3*j2*l2*i3*j3*i4 + k2*l2*i3*j3*i4
            + l2^2*i3*j3*i4 + 2*i1*i3^2*j3*i4 + 3*i1^4*j3^2*i4 + 2*i2^2*j3^2*i4 + 5*i1^2*j2*j3^2*i4 + 2*i2*j2*j3^2*i4 +
            3*j2^2*j3^2*i4 + 5*i2*k2*j3^2*i4 + j2*k2*j3^2*i4 + 6*i1^2*l2*j3^2*i4 + 2*i2*l2*j3^2*i4 + 4*j2*l2*j3^2*i4 +
            5*k2*l2*j3^2*i4 + 5*l2^2*j3^2*i4 + 4*i1*i3*j3^2*i4 + 6*i1*j3^3*i4 + 5*i1^7*k3*i4 + 5*i1^5*i2*k3*i4 +
            3*i1^3*i2^2*k3*i4 + i1*i2^3*k3*i4 + 6*i1^5*j2*k3*i4 + 2*i1^3*i2*j2*k3*i4 + 3*i1*i2^2*j2*k3*i4 +
            i1^3*j2^2*k3*i4 + 5*i1*i2*j2^2*k3*i4 + 5*i1*j2^3*k3*i4 + 4*i1^3*i2*k2*k3*i4 + 4*i1*i2^2*k2*k3*i4 +
            6*i1^3*j2*k2*k3*i4 + 3*i1*i2*j2*k2*k3*i4 + 5*i1*j2^2*k2*k3*i4 + 4*i1^3*k2^2*k3*i4 + 4*i1*i2*k2^2*k3*i4 +
            6*i1*j2*k2^2*k3*i4 + 4*i1*k2^3*k3*i4 + i1^5*l2*k3*i4 + 6*i1^3*i2*l2*k3*i4 + i1*i2^2*l2*k3*i4 +
            5*i1*i2*j2*l2*k3*i4 + 5*i1*j2^2*l2*k3*i4 + 2*i1^3*k2*l2*k3*i4 + 3*i1*j2*k2*l2*k3*i4 + 3*i1*k2^2*l2*k3*i4 +
            4*i1^3*l2^2*k3*i4 + 4*i1*i2*l2^2*k3*i4 + 3*i1*j2*l2^2*k3*i4 + 6*i1*k2*l2^2*k3*i4 + i1*l2^3*k3*i4 +
            4*i1^4*j3*k3*i4 + 4*i1^2*i2*j3*k3*i4 + 3*i1^2*j2*j3*k3*i4 + 3*i2*j2*j3*k3*i4 + 4*j2^2*j3*k3*i4 +
            5*i1^2*k2*j3*k3*i4 + 5*i2*k2*j3*k3*i4 + 5*k2^2*j3*k3*i4 + 5*i1^2*l2*j3*k3*i4 + 6*i2*l2*j3*k3*i4 +
            5*j2*l2*j3*k3*i4 + 4*k2*l2*j3*k3*i4 + 6*l2^2*j3*k3*i4 + i1^10*j4 + 6*i1^8*i2*j4 + 3*i1^6*i2^2*j4 +
            4*i1^4*i2^3*j4 + 6*i1^2*i2^4*j4 + 4*i1^6*i2*j2*j4 + 4*i1^2*i2^3*j2*j4 + 4*i1^6*j2^2*j4 + 2*i1^4*i2*j2^2*j4 +
            i1^2*i2^2*j2^2*j4 + i1^4*j2^3*j4 + 4*i1^2*i2*j2^3*j4 + 6*i1^2*j2^4*j4 + 4*i1^8*k2*j4 + i1^6*i2*k2*j4 +
            4*i1^2*i2^3*k2*j4 + 2*i2^4*k2*j4 + 4*i1^4*i2*j2*k2*j4 + 2*i1^2*i2^2*j2*k2*j4 + 6*i2^3*j2*k2*j4 +
            i1^4*j2^2*k2*j4 + 5*i1^2*i2*j2^2*k2*j4 + 5*i2^2*j2^2*k2*j4 + 3*i1^2*j2^3*k2*j4 + 6*i2*j2^3*k2*j4 +
            2*j2^4*k2*j4 + 6*i1^6*k2^2*j4 + i1^4*i2*k2^2*j4 + 6*i1^2*i2^2*k2^2*j4 + 4*i2^3*k2^2*j4 + i1^2*i2*j2*k2^2*j4 +
            3*i2^2*j2*k2^2*j4 + i1^2*j2^2*k2^2*j4 + 3*i2*j2^2*k2^2*j4 + 4*j2^3*k2^2*j4 + 4*i1^4*k2^3*j4 +
            6*i1^2*i2*k2^3*j4 + i2^2*k2^3*j4 + 6*i1^2*j2*k2^3*j4 + 3*i2*j2*k2^3*j4 + 4*i1^2*k2^4*j4 + 2*i2*k2^4*j4 +
            j2*k2^4*j4 + 6*k2^5*j4 + i1^8*l2*j4 + 4*i1^4*i2^2*l2*j4 + i1^2*i2^3*l2*j4 + 5*i1^6*j2*l2*j4 + i1^4*i2*j2*l2*j4
            + 2*i1^2*i2^2*j2*l2*j4 + 2*i1^4*j2^2*l2*j4 + 4*i1^2*j2^3*l2*j4 + 2*i1^6*k2*l2*j4 + 2*i1^4*i2*k2*l2*j4 +
            6*i1^2*i2^2*k2*l2*j4 + 5*i2^3*k2*l2*j4 + 3*i1^2*i2*j2*k2*l2*j4 + i2^2*j2*k2*l2*j4 + 3*i1^2*j2^2*k2*l2*j4 +
            4*i2*j2^2*k2*l2*j4 + 4*j2^3*k2*l2*j4 + 3*i1^4*k2^2*l2*j4 + 2*i1^2*i2*k2^2*l2*j4 + 3*i2^2*k2^2*l2*j4 +
            5*i1^2*k2^3*l2*j4 + 6*i2*k2^3*l2*j4 + 5*j2*k2^3*l2*j4 + 2*k2^4*l2*j4 + 6*i1^4*i2*l2^2*j4 + 2*i1^2*i2^2*l2^2*j4
            + i2^3*l2^2*j4 + 2*i1^4*j2*l2^2*j4 + i1^2*i2*j2*l2^2*j4 + 2*i2^2*j2*l2^2*j4 + 2*i1^2*j2^2*l2^2*j4 +
            4*j2^3*l2^2*j4 + 6*i1^4*k2*l2^2*j4 + 3*i1^2*i2*k2*l2^2*j4 + 5*i2^2*k2*l2^2*j4 + 5*i1^2*j2*k2*l2^2*j4 +
            3*i2*j2*k2*l2^2*j4 + j2^2*k2*l2^2*j4 + i1^2*k2^2*l2^2*j4 + 3*i2*k2^2*l2^2*j4 + j2*k2^2*l2^2*j4 +
            4*k2^3*l2^2*j4 + 2*i1^4*l2^3*j4 + i1^2*i2*l2^3*j4 + 6*i2^2*l2^3*j4 + 3*i1^2*j2*l2^3*j4 + 6*i2*j2*l2^3*j4 +
            4*j2^2*l2^3*j4 + 3*i1^2*k2*l2^3*j4 + i2*k2*l2^3*j4 + 2*j2*k2*l2^3*j4 + 5*k2^2*l2^3*j4 + 5*i1^2*l2^4*j4 +
            2*i2*l2^4*j4 + 2*j2*l2^4*j4 + 4*k2*l2^4*j4 + 3*l2^5*j4 + i1^7*j3*j4 + 6*i1^5*i2*j3*j4 + 3*i1^3*i2^2*j3*j4 +
            2*i1*i2^3*j3*j4 + 2*i1^5*j2*j3*j4 + 6*i1^3*i2*j2*j3*j4 + i1*i2^2*j2*j3*j4 + 5*i1^3*j2^2*j3*j4 +
            6*i1*i2*j2^2*j3*j4 + 5*i1*j2^3*j3*j4 + 4*i1^3*i2*k2*j3*j4 + 4*i1*i2^2*k2*j3*j4 + 5*i1^3*j2*k2*j3*j4 +
            3*i1*j2^2*k2*j3*j4 + 2*i1^3*k2^2*j3*j4 + 4*i1*i2*k2^2*j3*j4 + 3*i1*j2*k2^2*j3*j4 + 4*i1*k2^3*j3*j4 +
            3*i1^5*l2*j3*j4 + 6*i1*i2^2*l2*j3*j4 + 2*i1^3*j2*l2*j3*j4 + 2*i1*i2*j2*l2*j3*j4 + 6*i1*j2^2*l2*j3*j4 +
            4*i1^3*k2*l2*j3*j4 + 6*i1*j2*k2*l2*j3*j4 + 6*i1*k2^2*l2*j3*j4 + 6*i1^3*l2^2*j3*j4 + i1*i2*l2^2*j3*j4 +
            6*i1*j2*l2^2*j3*j4 + 5*i1*k2*l2^2*j3*j4 + 2*i1*l2^3*j3*j4 + 3*i1^4*j3^2*j4 + 6*i1^2*i2*j3^2*j4 +
            2*i2^2*j3^2*j4 + 3*i2*j2*j3^2*j4 + 2*j2^2*j3^2*j4 + 2*i1^2*k2*j3^2*j4 + i2*k2*j3^2*j4 + 2*j2*k2*j3^2*j4 +
            4*k2^2*j3^2*j4 + i1^2*l2*j3^2*j4 + 4*j2*l2*j3^2*j4 + 4*k2*l2*j3^2*j4 + l2^2*j3^2*j4 + 4*i1*j3^3*j4 +
            4*i1^10*k4 + 5*i1^6*i2^2*k4 + 5*i1^4*i2^3*k4 + 5*i1^8*j2*k4 + 5*i1^6*i2*j2*k4 + 6*i1^4*i2^2*j2*k4 +
            4*i1^6*j2^2*k4 + i1^4*i2*j2^2*k4 + 2*i1^4*j2^3*k4 + 4*i1^8*k2*k4 + 6*i1^6*i2*k2*k4 + 5*i1^4*i2^2*k2*k4 +
            5*i1^2*i2^3*k2*k4 + 4*i1^6*j2*k2*k4 + 4*i1^4*i2*j2*k2*k4 + 6*i1^2*i2^2*j2*k2*k4 + 5*i1^4*j2^2*k2*k4 +
            i1^2*i2*j2^2*k2*k4 + 2*i1^2*j2^3*k2*k4 + 3*i1^6*k2^2*k4 + 2*i1^4*i2*k2^2*k4 + 6*i1^2*i2^2*k2^2*k4 +
            5*i2^3*k2^2*k4 + 5*i1^4*j2*k2^2*k4 + i1^2*i2*j2*k2^2*k4 + 6*i2^2*j2*k2^2*k4 + i2*j2^2*k2^2*k4 + 2*j2^3*k2^2*k4
            + 4*i1^4*k2^3*k4 + 3*i2^2*k2^3*k4 + 5*i1^2*j2*k2^3*k4 + 6*i2*j2*k2^3*k4 + j2^2*k2^3*k4 + 3*i1^2*k2^4*k4 +
            3*i2*k2^4*k4 + 5*j2*k2^4*k4 + 2*k2^5*k4 + 3*i1^8*l2*k4 + 3*i1^6*i2*l2*k4 + i1^4*i2^2*l2*k4 + 2*i1^2*i2^3*l2*k4
            + 6*i1^6*j2*l2*k4 + i1^2*i2^2*j2*l2*k4 + 6*i1^4*j2^2*l2*k4 + 6*i1^2*i2*j2^2*l2*k4 + 5*i1^2*j2^3*l2*k4 +
            3*i1^6*k2*l2*k4 + 3*i1^4*i2*k2*l2*k4 + 4*i1^2*i2^2*k2*l2*k4 + 5*i2^3*k2*l2*k4 + 5*i1^4*j2*k2*l2*k4 +
            5*i1^2*i2*j2*k2*l2*k4 + 6*i2^2*j2*k2*l2*k4 + 5*i1^2*j2^2*k2*l2*k4 + i2*j2^2*k2*l2*k4 + 2*j2^3*k2*l2*k4 +
            4*i1^4*k2^2*l2*k4 + 2*i1^2*i2*k2^2*l2*k4 + 6*i2^2*k2^2*l2*k4 + 4*i1^2*j2*k2^2*l2*k4 + 3*i2*j2*k2^2*l2*k4 +
            5*j2^2*k2^2*l2*k4 + 6*i1^2*k2^3*l2*k4 + 2*i2*k2^3*l2*k4 + j2*k2^3*l2*k4 + 6*k2^4*l2*k4 + 4*i1^6*l2^2*k4 +
            3*i1^4*i2*l2^2*k4 + 4*i1^2*i2^2*l2^2*k4 + 3*i2^3*l2^2*k4 + 5*i1^4*j2*l2^2*k4 + i1^2*i2*j2*l2^2*k4 +
            5*i2^2*j2*l2^2*k4 + 2*i1^2*j2^2*l2^2*k4 + 2*i2*j2^2*l2^2*k4 + 4*j2^3*l2^2*k4 + 3*i1^4*k2*l2^2*k4 +
            5*i1^2*i2*k2*l2^2*k4 + 3*i2^2*k2*l2^2*k4 + 3*i2*j2*k2*l2^2*k4 + j2^2*k2*l2^2*k4 + 6*i1^2*k2^2*l2^2*k4 +
            3*i2*k2^2*l2^2*k4 + j2*k2^2*l2^2*k4 + 5*k2^3*l2^2*k4 + 2*i1^4*l2^3*k4 + 2*i1^2*i2*l2^3*k4 + 2*i2^2*l2^3*k4 +
            5*j2^2*l2^3*k4 + 5*j2*k2*l2^3*k4 + k2^2*l2^3*k4 + i1^2*l2^4*k4 + i2*l2^4*k4 + 3*k2*l2^4*k4 + 3*l2^5*k4 +
            4*i1^7*i3*k4 + 6*i1^5*i2*i3*k4 + i1^5*j2*i3*k4 + 3*i1^5*k2*i3*k4 + 3*i1^3*i2*k2*i3*k4 + 4*i1^3*j2*k2*i3*k4 +
            5*i1^3*k2^2*i3*k4 + 4*i1*i2*k2^2*i3*k4 + 3*i1*j2*k2^2*i3*k4 + 6*i1*k2^3*i3*k4 + i1^5*l2*i3*k4 +
            i1^3*i2*l2*i3*k4 + 6*i1^3*j2*l2*i3*k4 + 2*i1^3*k2*l2*i3*k4 + 4*i1*i2*k2*l2*i3*k4 + 3*i1*j2*k2*l2*i3*k4 +
            5*i1*k2^2*l2*i3*k4 + 5*i1^3*l2^2*i3*k4 + 2*i1*i2*l2^2*i3*k4 + 5*i1*j2*l2^2*i3*k4 + 4*i1*k2*l2^2*i3*k4 +
            6*i1*l2^3*i3*k4 + i1^7*k3*k4 + 3*i1^5*i2*k3*k4 + 4*i1^3*i2^2*k3*k4 + 2*i1^5*j2*k3*k4 + 6*i1^3*i2*j2*k3*k4 +
            4*i1^3*j2^2*k3*k4 + i1^5*k2*k3*k4 + 2*i1^3*i2*k2*k3*k4 + 6*i1*i2^2*k2*k3*k4 + 5*i1^3*j2*k2*k3*k4 +
            2*i1*i2*j2*k2*k3*k4 + 6*i1*j2^2*k2*k3*k4 + 3*i1^3*k2^2*k3*k4 + 4*i1*i2*k2^2*k3*k4 + 4*i1*j2*k2^2*k3*k4 +
            4*i1*k2^3*k3*k4 + 6*i1^5*l2*k3*k4 + 6*i1*i2^2*l2*k3*k4 + 6*i1^3*j2*l2*k3*k4 + 2*i1*i2*j2*l2*k3*k4 +
            6*i1*j2^2*l2*k3*k4 + 2*i1^3*k2*l2*k3*k4 + i1*i2*k2*l2*k3*k4 + i1*j2*k2*l2*k3*k4 + 6*i1*k2^2*l2*k3*k4 +
            i1^3*l2^2*k3*k4 + 6*i1*i2*l2^2*k3*k4 + 5*i1*j2*l2^2*k3*k4 + 5*i1*k2*l2^2*k3*k4 + 5*i1*l2^3*k3*k4 +
            i1^4*k3^2*k4 + i1^2*i2*k3^2*k4 + 6*i1^2*j2*k3^2*k4 + 6*i1^2*k2*k3^2*k4 + 3*i2*k2*k3^2*k4 + 4*j2*k2*k3^2*k4 +
            6*k2^2*k3^2*k4 + i1^2*l2*k3^2*k4 + 5*i2*l2*k3^2*k4 + 2*j2*l2*k3^2*k4 + 3*k2*l2*k3^2*k4 + l2^2*k3^2*k4 +
            5*i1*k3^3*k4 + 2*i1^6*k4^2 + i1^4*k2*k4^2 + 2*i1^2*k2^2*k4^2 + 3*k2^3*k4^2 + 2*i1^4*l2*k4^2 +
            6*i1^2*k2*l2*k4^2 + 2*k2^2*l2*k4^2 + 3*i1^2*l2^2*k4^2 + 5*l2^3*k4^2 + 4*i1^10*l4 + i1^8*i2*l4 + 2*i1^6*i2^2*l4
            + 3*i1^4*i2^3*l4 + 3*i1^2*i2^4*l4 + i2^5*l4 + 5*i1^8*j2*l4 + i1^6*i2*j2*l4 + 3*i1^4*i2^2*j2*l4 + 2*i2^4*j2*l4
            + 2*i1^6*j2^2*l4 + 2*i1^4*i2*j2^2*l4 + 3*i1^2*i2^2*j2^2*l4 + 3*i2^3*j2^2*l4 + 6*i1^4*j2^3*l4 +
            3*i1^2*i2*j2^3*l4 + 4*i2^2*j2^3*l4 + 5*i1^2*j2^4*l4 + 5*i2*j2^4*l4 + 6*j2^5*l4 + i1^8*k2*l4 + 6*i1^6*i2*k2*l4
            + 5*i1^2*i2^3*k2*l4 + 2*i2^4*k2*l4 + 2*i1^6*j2*k2*l4 + 2*i1^4*i2*j2*k2*l4 + 2*i1^2*i2^2*j2*k2*l4 +
            i2^3*j2*k2*l4 + 6*i1^4*j2^2*k2*l4 + 5*i1^2*i2*j2^2*k2*l4 + 6*i2^2*j2^2*k2*l4 + 2*i1^2*j2^3*k2*l4 +
            5*i2*j2^3*k2*l4 + i1^6*k2^2*l4 + 6*i1^2*i2^2*k2^2*l4 + 5*i2^3*k2^2*l4 + 6*i1^4*j2*k2^2*l4 +
            4*i1^2*i2*j2*k2^2*l4 + 3*i2^2*j2*k2^2*l4 + 3*i1^2*j2^2*k2^2*l4 + i2*j2^2*k2^2*l4 + 5*j2^3*k2^2*l4 +
            5*i1^4*k2^3*l4 + 3*i1^2*i2*k2^3*l4 + 5*i1^2*j2*k2^3*l4 + i2*j2*k2^3*l4 + j2^2*k2^3*l4 + 4*i1^2*k2^4*l4 +
            5*i2*k2^4*l4 + j2*k2^4*l4 + 4*k2^5*l4 + 6*i1^8*l2*l4 + i1^6*i2*l2*l4 + 4*i1^2*i2^3*l2*l4 + 4*i2^4*l2*l4 +
            4*i1^6*j2*l2*l4 + 2*i1^4*i2*j2*l2*l4 + 4*i1^2*i2^2*j2*l2*l4 + 3*i2^3*j2*l2*l4 + 2*i1^4*j2^2*l2*l4 +
            4*i1^2*i2*j2^2*l2*l4 + 2*i2^2*j2^2*l2*l4 + 2*i1^2*j2^3*l2*l4 + 6*i2*j2^3*l2*l4 + 6*j2^4*l2*l4 +
            6*i1^6*k2*l2*l4 + 3*i1^4*i2*k2*l2*l4 + 5*i1^2*i2^2*k2*l2*l4 + i2^3*k2*l2*l4 + 4*i1^4*j2*k2*l2*l4 +
            2*i1^2*i2*j2*k2*l2*l4 + 5*i2^2*j2*k2*l2*l4 + i1^2*j2^2*k2*l2*l4 + 4*i2*j2^2*k2*l2*l4 + 4*j2^3*k2*l2*l4 +
            i1^2*i2*k2^2*l2*l4 + 3*i2^2*k2^2*l2*l4 + 6*i1^2*j2*k2^2*l2*l4 + 4*i2*j2*k2^2*l2*l4 + 4*j2^2*k2^2*l2*l4 +
            3*i2*k2^3*l2*l4 + j2*k2^3*l2*l4 + 2*k2^4*l2*l4 + i1^2*i2^2*l2^2*l4 + 6*i2^3*l2^2*l4 + 5*i1^2*i2*j2*l2^2*l4 +
            4*i2^2*j2*l2^2*l4 + 3*i1^2*j2^2*l2^2*l4 + 6*i2*j2^2*l2^2*l4 + j2^3*l2^2*l4 + i1^4*k2*l2^2*l4 +
            2*i1^2*i2*k2*l2^2*l4 + 5*i2^2*k2*l2^2*l4 + i1^2*j2*k2*l2^2*l4 + 3*i2*j2*k2*l2^2*l4 + i1^2*k2^2*l2^2*l4 +
            5*i2*k2^2*l2^2*l4 + 2*j2*k2^2*l2^2*l4 + 5*k2^3*l2^2*l4 + 6*i1^4*l2^3*l4 + 4*i1^2*i2*l2^3*l4 + 2*i2^2*l2^3*l4 +
            2*i1^2*j2*l2^3*l4 + 3*i2*j2*l2^3*l4 + 2*j2^2*l2^3*l4 + 2*i1^2*k2*l2^3*l4 + i2*k2*l2^3*l4 + 2*j2*k2*l2^3*l4 +
            6*k2^2*l2^3*l4 + 5*i1^2*l2^4*l4 + 5*i2*l2^4*l4 + 6*j2*l2^4*l4 + 6*l2^5*l4 + 2*i1^7*j3*l4 + i1^5*i2*j3*l4 +
            4*i1^3*i2^2*j3*l4 + 5*i1*i2^3*j3*l4 + 5*i1^3*i2*j2*j3*l4 + 3*i1*i2^2*j2*j3*l4 + 5*i1^3*j2^2*j3*l4 +
            6*i1*j2^3*j3*l4 + 5*i1^5*k2*j3*l4 + 2*i1^3*i2*k2*j3*l4 + 5*i1*i2^2*k2*j3*l4 + 5*i1*i2*j2*k2*j3*l4 +
            2*i1^3*k2^2*j3*l4 + 5*i1*i2*k2^2*j3*l4 + 4*i1*j2*k2^2*j3*l4 + 5*i1*k2^3*j3*l4 + 6*i1^5*l2*j3*l4 +
            i1*i2^2*l2*j3*l4 + 2*i1^3*j2*l2*j3*l4 + 6*i1*i2*j2*l2*j3*l4 + 5*i1*j2^2*l2*j3*l4 + i1^3*k2*l2*j3*l4 +
            6*i1*i2*k2*l2*j3*l4 + 5*i1*j2*k2*l2*j3*l4 + 4*i1*k2^2*l2*j3*l4 + 3*i1^3*l2^2*j3*l4 + 5*i1*i2*l2^2*j3*l4 +
            4*i1*j2*l2^2*j3*l4 + 2*i1*l2^3*j3*l4 + 2*i1^4*j3^2*l4 + 4*i2^2*j3^2*l4 + 3*i1^2*j2*j3^2*l4 + i2*j2*j3^2*l4 +
            2*j2^2*j3^2*l4 + 4*i1^2*k2*j3^2*l4 + 3*i2*k2*j3^2*l4 + 6*j2*k2*j3^2*l4 + 3*k2^2*j3^2*l4 + 2*i1^2*l2*j3^2*l4 +
            3*i2*l2*j3^2*l4 + j2*l2*j3^2*l4 + 5*k2*l2*j3^2*l4 + l2^2*j3^2*l4 + 5*i1*j3^3*l4 + 2*i1^7*k3*l4 +
            3*i1^5*i2*k3*l4 + 2*i1*i2^3*k3*l4 + i1^3*i2*j2*k3*l4 + 5*i1*i2^2*j2*k3*l4 + i1^3*j2^2*k3*l4 +
            5*i1*i2*j2^2*k3*l4 + 2*i1*j2^3*k3*l4 + 6*i1^5*k2*k3*l4 + i1^3*i2*k2*k3*l4 + 2*i1*i2^2*k2*k3*l4 +
            3*i1^3*j2*k2*k3*l4 + 3*i1*i2*j2*k2*k3*l4 + i1*j2^2*k2*k3*l4 + 5*i1*i2*k2^2*k3*l4 + 2*i1*j2*k2^2*k3*l4 +
            2*i1*k2^3*k3*l4 + 6*i1^5*l2*k3*l4 + 2*i1^3*i2*l2*k3*l4 + i1^3*j2*l2*k3*l4 + i1*i2*j2*l2*k3*l4 +
            i1*j2^2*l2*k3*l4 + 5*i1^3*k2*l2*k3*l4 + 3*i1*i2*k2*l2*k3*l4 + 4*i1*j2*k2*l2*k3*l4 + i1*k2^2*l2*k3*l4 +
            4*i1^3*l2^2*k3*l4 + i1*i2*l2^2*k3*l4 + 2*i1*j2*l2^2*k3*l4 + 5*i1*k2*l2^2*k3*l4 + i1*l2^3*k3*l4 +
            4*i1^4*j3*k3*l4 + 4*i2^2*j3*k3*l4 + 3*i1^2*j2*j3*k3*l4 + i2*j2*j3*k3*l4 + 2*j2^2*j3*k3*l4 + 3*i1^2*k2*j3*k3*l4
            + 2*i2*k2*j3*k3*l4 + j2*k2*j3*k3*l4 + 5*k2^2*j3*k3*l4 + 2*i2*l2*j3*k3*l4 + 4*j2*l2*j3*k3*l4 + 2*k2*l2*j3*k3*l4
            + 3*l2^2*j3*k3*l4 + 3*i1*j3^2*k3*l4 + 3*i1^4*k3^2*l4 + 2*i2^2*k3^2*l4 + 5*i1^2*j2*k3^2*l4 + 4*i2*j2*k3^2*l4 +
            j2^2*k3^2*l4 + 4*i1^2*k2*k3^2*l4 + 3*i2*k2*k3^2*l4 + 2*j2*k2*k3^2*l4 + 6*k2^2*k3^2*l4 + 3*i1^2*l2*k3^2*l4 +
            3*i2*l2*k3^2*l4 + 2*j2*l2*k3^2*l4 + 6*k2*l2*k3^2*l4 + 6*l2^2*k3^2*l4 + i1*j3*k3^2*l4 + 5*i1*k3^3*l4 +
            4*i1^7*l3*l4 + 2*i1^5*i2*l3*l4 + 2*i1^3*i2^2*l3*l4 + i1*i2^3*l3*l4 + 4*i1^5*j2*l3*l4 + 2*i1^3*i2*j2*l3*l4 +
            2*i1*i2^2*j2*l3*l4 + 3*i1^3*j2^2*l3*l4 + 4*i1*i2*j2^2*l3*l4 + 5*i1^5*k2*l3*l4 + 2*i1^3*i2*k2*l3*l4 +
            3*i1*i2^2*k2*l3*l4 + i1^3*j2*k2*l3*l4 + 2*i1*j2^2*k2*l3*l4 + i1^3*k2^2*l3*l4 + 3*i1*i2*k2^2*l3*l4 +
            i1*j2*k2^2*l3*l4 + i1^5*l2*l3*l4 + 2*i1^3*i2*l2*l3*l4 + 5*i1*i2^2*l2*l3*l4 + i1^3*j2*l2*l3*l4 +
            i1*i2*j2*l2*l3*l4 + 5*i1*j2^2*l2*l3*l4 + 5*i1^3*k2*l2*l3*l4 + 5*i1*i2*k2*l2*l3*l4 + 3*i1*j2*k2*l2*l3*l4 +
            i1*k2^2*l2*l3*l4 + 2*i1*i2*l2^2*l3*l4 + 5*i1*j2*l2^2*l3*l4 + 5*i1*k2*l2^2*l3*l4 + i1*l2^3*l3*l4 +
            6*i1^4*j3*l3*l4 + 6*i1^2*i2*j3*l3*l4 + 5*i2^2*j3*l3*l4 + 5*i1^2*j2*j3*l3*l4 + 6*i2*j2*j3*l3*l4 +
            5*j2^2*j3*l3*l4 + 5*i1^2*k2*j3*l3*l4 + j2*k2*j3*l3*l4 + 4*k2^2*j3*l3*l4 + 5*i1^2*l2*j3*l3*l4 +
            3*i2*l2*j3*l3*l4 + 4*j2*l2*j3*l3*l4 + 6*k2*l2*j3*l3*l4 + 4*l2^2*j3*l3*l4 + 3*i1*j3^2*l3*l4 + 3*i1^4*l3^2*l4 +
            4*i1^2*i2*l3^2*l4 + 2*i2^2*l3^2*l4 + 5*i2*j2*l3^2*l4 + 2*j2^2*l3^2*l4 + 3*i2*k2*l3^2*l4 + 3*j2*k2*l3^2*l4 +
            k2^2*l3^2*l4 + 4*i1^2*l2*l3^2*l4 + 2*i2*l2*l3^2*l4 + j2*l2*l3^2*l4 + 3*l2^2*l3^2*l4 + 3*i1*j3*l3^2*l4 +
            4*i1*l3^3*l4 + 5*i1^6*l4^2 + 5*i1^4*i2*l4^2 + 2*i2^3*l4^2 + 5*i1^4*j2*l4^2 + i1^2*j2^2*l4^2 + i2*j2^2*l4^2 +
            4*j2^3*l4^2 + 6*i1^4*k2*l4^2 + 5*i1^2*i2*k2*l4^2 + 5*i2^2*k2*l4^2 + i1^2*j2*k2*l4^2 + i2*j2*k2*l4^2 +
            2*j2^2*k2*l4^2 + 2*j2*k2^2*l4^2 + 4*k2^3*l4^2 + 4*i1^4*l2*l4^2 + i1^2*i2*l2*l4^2 + 6*i2^2*l2*l4^2 +
            3*i1^2*j2*l2*l4^2 + 4*i2*j2*l2*l4^2 + 6*j2^2*l2*l4^2 + 2*i1^2*k2*l2*l4^2 + 6*j2*k2*l2*l4^2 + 2*k2^2*l2*l4^2 +
            5*i1^2*l2^2*l4^2 + 4*i2*l2^2*l4^2 + j2*l2^2*l4^2 + 2*k2*l2^2*l4^2 + 5*l2^3*l4^2 + 3*i1^3*j3*l4^2 +
            5*i1*i2*j3*l4^2 + i1*j2*j3*l4^2 + 3*i1*l2*j3*l4^2 + 2*j3^2*l4^2 + 5*i1^3*k3*l4^2 + 5*i1*i2*k3*l4^2 +
            2*i1*j2*k3*l4^2 + 6*i1*k2*k3*l4^2 + 5*i1*l2*k3*l4^2 + 3*j3*k3*l4^2 + 3*i1*i2*l3*l4^2 + i1*j2*l3*l4^2 +
            4*i1*k2*l3*l4^2 + 4*i1*l2*l3*l4^2 + 6*j3*l3*l4^2 + 4*l3^2*l4^2 + 6*i1^2*l4^3 + 2*i2*l4^3 + 2*k2*l4^3 +
            6*l2*l4^3 + 3*i1^10*m4 + 6*i1^8*i2*m4 + 5*i1^6*i2^2*m4 + 2*i1^2*i2^4*m4 + 3*i2^5*m4 + i1^8*j2*m4 +
            6*i1^6*i2*j2*m4 + 5*i1^4*i2^2*j2*m4 + 4*i1^2*i2^3*j2*m4 + 6*i2^4*j2*m4 + 5*i1^6*j2^2*m4 + i1^4*i2*j2^2*m4 +
            4*i1^2*i2^2*j2^2*m4 + 2*i2^3*j2^2*m4 + i1^4*j2^3*m4 + 5*i2^2*j2^3*m4 + 4*i1^2*j2^4*m4 + i2*j2^4*m4 + 4*j2^5*m4
            + i1^6*i2*k2*m4 + i1^4*i2^2*k2*m4 + 6*i1^2*i2^3*k2*m4 + i2^4*k2*m4 + 2*i1^6*j2*k2*m4 + i1^4*i2*j2*k2*m4 +
            2*i2^3*j2*k2*m4 + 6*i1^4*j2^2*k2*m4 + 2*i1^2*i2*j2^2*k2*m4 + 2*i2^2*j2^2*k2*m4 + 6*i1^2*j2^3*k2*m4 +
            2*j2^4*k2*m4 + i1^6*k2^2*m4 + 5*i2^3*k2^2*m4 + i1^4*j2*k2^2*m4 + 4*i1^2*i2*j2*k2^2*m4 + 4*i2^2*j2*k2^2*m4 +
            2*i1^2*j2^2*k2^2*m4 + 6*i2*j2^2*k2^2*m4 + 6*j2^3*k2^2*m4 + 2*i1^4*k2^3*m4 + 4*i1^2*i2*k2^3*m4 + 4*i2^2*k2^3*m4
            + 3*i1^2*j2*k2^3*m4 + 6*i2*j2*k2^3*m4 + 3*j2^2*k2^3*m4 + 5*i1^2*k2^4*m4 + 6*i2*k2^4*m4 + 6*j2*k2^4*m4 +
            6*i1^2*i2^3*l2*m4 + 4*i1^6*j2*l2*m4 + 4*i1^4*i2*j2*l2*m4 + 3*i1^2*i2^2*j2*l2*m4 + 5*i1^4*j2^2*l2*m4 +
            4*i1^2*i2*j2^2*l2*m4 + i1^2*j2^3*l2*m4 + 3*i1^6*k2*l2*m4 + 4*i1^4*i2*k2*l2*m4 + 6*i1^2*i2^2*k2*l2*m4 +
            5*i1^4*j2*k2*l2*m4 + 4*i2^2*j2*k2*l2*m4 + 3*i1^2*j2^2*k2*l2*m4 + 3*i2*j2^2*k2*l2*m4 + 5*i1^4*k2^2*l2*m4 +
            2*i1^2*i2*k2^2*l2*m4 + 5*i2^2*k2^2*l2*m4 + 4*i1^2*j2*k2^2*l2*m4 + 2*i2*j2*k2^2*l2*m4 + 2*j2^2*k2^2*l2*m4 +
            2*i1^2*k2^3*l2*m4 + 4*i2*k2^3*l2*m4 + 4*j2*k2^3*l2*m4 + 6*k2^4*l2*m4 + 3*i1^6*l2^2*m4 + 4*i1^4*i2*l2^2*m4 +
            6*i1^2*i2^2*l2^2*m4 + 6*i1^4*j2*l2^2*m4 + 6*i1^2*i2*j2*l2^2*m4 + 4*i2^2*j2*l2^2*m4 + 6*i1^2*j2^2*l2^2*m4 +
            6*i2*j2^2*l2^2*m4 + 4*j2^3*l2^2*m4 + 5*i1^4*k2*l2^2*m4 + 5*i1^2*i2*k2*l2^2*m4 + i2^2*k2*l2^2*m4 +
            2*i1^2*j2*k2*l2^2*m4 + 3*i2*j2*k2*l2^2*m4 + 2*j2^2*k2*l2^2*m4 + 2*i1^2*k2^2*l2^2*m4 + 2*i2*k2^2*l2^2*m4 +
            3*k2^3*l2^2*m4 + 2*i1^4*l2^3*m4 + 3*i1^2*i2*l2^3*m4 + 4*i1^2*j2*l2^3*m4 + 4*i2*j2*l2^3*m4 + 3*j2^2*l2^3*m4 +
            4*i1^2*k2*l2^3*m4 + 4*i2*k2*l2^3*m4 + 6*j2*k2*l2^3*m4 + 5*i1^2*l2^4*m4 + 3*i2*l2^4*m4 + 3*j2*l2^4*m4 +
            5*k2*l2^4*m4 + 3*l2^5*m4 + i1^7*i3*m4 + 6*i1^5*i2*i3*m4 + i1^3*i2^2*i3*m4 + 2*i1*i2^3*i3*m4 + i1^5*j2*i3*m4 +
            i1^3*i2*j2*i3*m4 + i1*i2^2*j2*i3*m4 + 5*i1^3*j2^2*i3*m4 + 6*i1*i2*j2^2*i3*m4 + 5*i1*j2^3*i3*m4 +
            3*i1^5*k2*i3*m4 + 3*i1^3*i2*k2*i3*m4 + i1*i2^2*k2*i3*m4 + 4*i1^3*j2*k2*i3*m4 + 6*i1*j2^2*k2*i3*m4 +
            6*i1^3*k2^2*i3*m4 + 4*i1*j2*k2^2*i3*m4 + 2*i1*k2^3*i3*m4 + 3*i1^5*l2*i3*m4 + 6*i1^3*i2*l2*i3*m4 +
            4*i1*i2^2*l2*i3*m4 + 6*i1^3*j2*l2*i3*m4 + 4*i1*i2*j2*l2*i3*m4 + 6*i1*j2^2*l2*i3*m4 + 3*i1^3*k2*l2*i3*m4 +
            i1*i2*k2*l2*i3*m4 + 5*i1*j2*k2*l2*i3*m4 + 4*i1*k2^2*l2*i3*m4 + 5*i1^3*l2^2*i3*m4 + 3*i1*k2*l2^2*i3*m4 +
            6*i1*l2^3*i3*m4 + 2*i1^2*i2*i3^2*m4 + 5*i1^2*j2*i3^2*m4 + 4*i1^2*k2*i3^2*m4 + i2*k2*i3^2*m4 + 6*j2*k2*i3^2*m4
            + k2^2*i3^2*m4 + i1^2*l2*i3^2*m4 + 4*i2*l2*i3^2*m4 + 3*j2*l2*i3^2*m4 + k2*l2*i3^2*m4 + 5*l2^2*i3^2*m4 +
            i1^7*j3*m4 + 6*i1^5*i2*j3*m4 + 4*i1^3*i2^2*j3*m4 + 6*i1^5*j2*j3*m4 + i1^3*i2*j2*j3*m4 + 4*i1*i2^2*j2*j3*m4 +
            2*i1^3*j2^2*j3*m4 + 6*i1*i2*j2^2*j3*m4 + 4*i1*j2^3*j3*m4 + 3*i1^5*k2*j3*m4 + 2*i1^3*i2*k2*j3*m4 +
            6*i1*i2^2*k2*j3*m4 + i1^3*j2*k2*j3*m4 + 3*i1*i2*j2*k2*j3*m4 + 2*i1*j2^2*k2*j3*m4 + i1^3*k2^2*j3*m4 +
            6*i1*i2*k2^2*j3*m4 + i1*j2*k2^2*j3*m4 + 6*i1*k2^3*j3*m4 + 5*i1^5*l2*j3*m4 + 5*i1^3*i2*l2*j3*m4 +
            5*i1*i2^2*l2*j3*m4 + 3*i1^3*j2*l2*j3*m4 + 4*i1*i2*j2*l2*j3*m4 + 5*i1^3*k2*l2*j3*m4 + 4*i1*i2*k2*l2*j3*m4 +
            i1*j2*k2*l2*j3*m4 + 5*i1*k2^2*l2*j3*m4 + 2*i1^3*l2^2*j3*m4 + 5*i1*i2*l2^2*j3*m4 + i1*j2*l2^2*j3*m4 +
            i1*k2*l2^2*j3*m4 + 3*i1*l2^3*j3*m4 + 2*i1^4*j3^2*m4 + 2*i1^2*i2*j3^2*m4 + 6*i2^2*j3^2*m4 + 4*i1^2*j2*j3^2*m4 +
            3*i2*j2*j3^2*m4 + 5*j2^2*j3^2*m4 + 3*i1^2*k2*j3^2*m4 + 2*i2*k2*j3^2*m4 + 5*j2*k2*j3^2*m4 + 6*k2^2*j3^2*m4 +
            i1^2*l2*j3^2*m4 + 2*j2*l2*j3^2*m4 + 5*k2*l2*j3^2*m4 + i1*j3^3*m4 + 3*i1^7*l3*m4 + 6*i1^5*i2*l3*m4 +
            2*i1^3*i2^2*l3*m4 + 4*i1*i2^3*l3*m4 + 5*i1*i2^2*j2*l3*m4 + 4*i1^3*j2^2*l3*m4 + i1*i2*j2^2*l3*m4 +
            4*i1*j2^3*l3*m4 + 2*i1^5*k2*l3*m4 + i1^3*i2*k2*l3*m4 + i1^3*j2*k2*l3*m4 + i1*i2*j2*k2*l3*m4 + i1*j2^2*k2*l3*m4
            + i1^3*k2^2*l3*m4 + 2*i1*j2*k2^2*l3*m4 + 2*i1*k2^3*l3*m4 + 4*i1^5*l2*l3*m4 + 4*i1^3*i2*l2*l3*m4 +
            5*i1*i2^2*l2*l3*m4 + 4*i1^3*j2*l2*l3*m4 + 6*i1*i2*j2*l2*l3*m4 + 5*i1*j2^2*l2*l3*m4 + 4*i1^3*k2*l2*l3*m4 +
            i1*i2*k2*l2*l3*m4 + 2*i1*k2^2*l2*l3*m4 + 6*i1*i2*l2^2*l3*m4 + i1*j2*l2^2*l3*m4 + 6*i1*l2^3*l3*m4 +
            i1^4*i3*l3*m4 + 4*i1^2*i2*i3*l3*m4 + 6*i2^2*i3*l3*m4 + 3*i2*j2*i3*l3*m4 + 5*j2^2*i3*l3*m4 + 4*i1^2*k2*i3*l3*m4
            + 2*i2*k2*i3*l3*m4 + 2*j2*k2*i3*l3*m4 + 5*k2^2*i3*l3*m4 + 5*i1^2*l2*i3*l3*m4 + 6*i2*l2*i3*l3*m4 +
            3*j2*l2*i3*l3*m4 + 2*k2*l2*i3*l3*m4 + 3*l2^2*i3*l3*m4 + i1*i3^2*l3*m4 + 6*i1^4*j3*l3*m4 + 5*i1^2*i2*j3*l3*m4 +
            2*i2^2*j3*l3*m4 + 2*i1^2*j2*j3*l3*m4 + 6*i2*j2*j3*l3*m4 + 6*j2^2*j3*l3*m4 + 6*i1^2*k2*j3*l3*m4 +
            2*i2*k2*j3*l3*m4 + 5*k2^2*j3*l3*m4 + 4*i1^2*l2*j3*l3*m4 + i2*l2*j3*l3*m4 + 2*j2*l2*j3*l3*m4 + 3*l2^2*j3*l3*m4
            + 4*i1*j3^2*l3*m4 + 2*i1^4*l3^2*m4 + 4*i1^2*i2*l3^2*m4 + 3*i2^2*l3^2*m4 + 3*i1^2*j2*l3^2*m4 + 2*j2^2*l3^2*m4 +
            6*i1^2*k2*l3^2*m4 + 6*j2*k2*l3^2*m4 + 2*k2^2*l3^2*m4 + 6*i1^2*l2*l3^2*m4 + 3*i2*l2*l3^2*m4 + 4*j2*l2*l3^2*m4 +
            3*k2*l2*l3^2*m4 + 2*l2^2*l3^2*m4 + i1*j3*l3^2*m4 + 4*i1*l3^3*m4 + 3*i1^9*i5 + i1^5*i2^2*i5 + 5*i1^7*j2*i5 +
            5*i1^5*i2*j2*i5 + i1^5*j2^2*i5 + 4*i1^7*k2*i5 + 6*i1^5*i2*k2*i5 + 4*i1^3*i2^2*k2*i5 + 3*i1^5*j2*k2*i5 +
            6*i1^3*i2*j2*k2*i5 + 4*i1^3*j2^2*k2*i5 + 2*i1^5*k2^2*i5 + 5*i1^3*i2*k2^2*i5 + 3*i1*i2^2*k2^2*i5 +
            i1*i2*j2*k2^2*i5 + 3*i1*j2^2*k2^2*i5 + 3*i1^3*k2^3*i5 + i1*i2*k2^3*i5 + 2*i1*j2*k2^3*i5 + 3*i1*k2^4*i5 +
            5*i1^7*l2*i5 + 6*i1^5*i2*l2*i5 + 3*i1^3*i2^2*l2*i5 + 6*i1^5*j2*l2*i5 + i1^3*i2*j2*l2*i5 + 3*i1^3*j2^2*l2*i5 +
            2*i1^3*i2*k2*l2*i5 + 5*i1*i2^2*k2*l2*i5 + 3*i1^3*j2*k2*l2*i5 + 4*i1*i2*j2*k2*l2*i5 + 5*i1*j2^2*k2*l2*i5 +
            3*i1^3*k2^2*l2*i5 + 6*i1*i2*k2^2*l2*i5 + 3*i1*j2*k2^2*l2*i5 + 5*i1*k2^3*l2*i5 + 3*i1^5*l2^2*i5 +
            4*i1^3*i2*l2^2*i5 + 2*i1*i2^2*l2^2*i5 + 3*i1*i2*j2*l2^2*i5 + 2*i1*j2^2*l2^2*i5 + 2*i1^3*k2*l2^2*i5 +
            i1*i2*k2*l2^2*i5 + 2*i1*j2*k2*l2^2*i5 + 5*i1*k2^2*l2^2*i5 + 2*i1^3*l2^3*i5 + 5*i1*j2*l2^3*i5 + i1*k2*l2^3*i5 +
            i1*l2^4*i5 + 5*i1^6*i3*i5 + 2*i1^4*k2*i3*i5 + 3*i1^2*k2^2*i3*i5 + 3*k2^3*i3*i5 + 6*i1^4*l2*i3*i5 +
            5*i1^2*k2*l2*i3*i5 + 2*k2^2*l2*i3*i5 + 5*k2*l2^2*i3*i5 + 3*l2^3*i3*i5 + 6*i1^6*j3*i5 + i1^4*i2*j3*i5 +
            6*i1^4*j2*j3*i5 + 2*i1^4*k2*j3*i5 + 4*i1^2*i2*k2*j3*i5 + 3*i1^2*j2*k2*j3*i5 + 6*i1^2*k2^2*j3*i5 +
            6*i2*k2^2*j3*i5 + 5*j2*k2^2*j3*i5 + 3*k2^3*j3*i5 + i1^4*l2*j3*i5 + 6*i1^2*i2*l2*j3*i5 + i1^2*j2*l2*j3*i5 +
            6*i1^2*k2*l2*j3*i5 + 4*i2*k2*l2*j3*i5 + 3*j2*k2*l2*j3*i5 + 4*i1^2*l2^2*j3*i5 + 5*i2*l2^2*j3*i5 +
            2*j2*l2^2*j3*i5 + 2*l2^3*j3*i5 + 6*i1^3*j3^2*i5 + i1*k2*j3^2*i5 + 4*i1*l2*j3^2*i5 + j3^3*i5 + 6*i1^7*i2*j5 +
            2*i1^5*i2^2*j5 + i1^3*i2^3*j5 + 2*i1*i2^4*j5 + i1^7*j2*j5 + 3*i1^5*i2*j2*j5 + 5*i1*i2^3*j2*j5 + 2*i1^5*j2^2*j5
            + i1^3*i2*j2^2*j5 + i1*i2^2*j2^2*j5 + 5*i1^3*j2^3*j5 + 3*i1*i2*j2^3*j5 + 3*i1*j2^4*j5 + 4*i1^7*k2*j5 +
            3*i1^5*i2*k2*j5 + 2*i1^3*i2^2*k2*j5 + 3*i1*i2^3*k2*j5 + 4*i1^5*j2*k2*j5 + 2*i1^3*i2*j2*k2*j5 +
            6*i1*i2^2*j2*k2*j5 + 4*i1^3*j2^2*k2*j5 + 6*i1*i2*j2^2*k2*j5 + 6*i1*j2^3*k2*j5 + 5*i1^5*k2^2*j5 +
            2*i1^3*i2*k2^2*j5 + 3*i1*i2^2*k2^2*j5 + 5*i1^3*j2*k2^2*j5 + 2*i1*i2*j2*k2^2*j5 + 5*i1*j2^2*k2^2*j5 +
            3*i1^3*k2^3*j5 + 2*i1*i2*k2^3*j5 + 6*i1*j2*k2^3*j5 + 6*i1*k2^4*j5 + 4*i1^5*i2*l2*j5 + i1^3*i2^2*l2*j5 +
            3*i1*i2^3*l2*j5 + i1^5*j2*l2*j5 + 6*i1*i2^2*j2*l2*j5 + 4*i1^3*j2^2*l2*j5 + 5*i1*i2*j2^2*l2*j5 +
            3*i1^5*k2*l2*j5 + 5*i1^3*i2*k2*l2*j5 + 5*i1*i2^2*k2*l2*j5 + 2*i1^3*j2*k2*l2*j5 + 4*i1*i2*j2*k2*l2*j5 +
            i1*j2^2*k2*l2*j5 + 3*i1^3*k2^2*l2*j5 + i1*i2*k2^2*l2*j5 + 6*i1*j2*k2^2*l2*j5 + 6*i1*k2^3*l2*j5 +
            2*i1^5*l2^2*j5 + 4*i1^3*i2*l2^2*j5 + 5*i1*i2^2*l2^2*j5 + i1^3*j2*l2^2*j5 + 5*i1*i2*j2*l2^2*j5 +
            i1*j2^2*l2^2*j5 + 6*i1*i2*k2*l2^2*j5 + 5*i1*j2*k2*l2^2*j5 + 3*i1*k2^2*l2^2*j5 + 3*i1^3*l2^3*j5 +
            5*i1*i2*l2^3*j5 + 2*i1*j2*l2^3*j5 + 6*i1*k2*l2^3*j5 + 2*i1*l2^4*j5 + 5*i1^6*i3*j5 + 2*i1^4*i2*i3*j5 +
            i1^2*i2^2*i3*j5 + 4*i2^3*i3*j5 + 2*i1^4*j2*i3*j5 + 5*i1^2*i2*j2*i3*j5 + 2*i2^2*j2*i3*j5 + i1^2*j2^2*i3*j5 +
            5*i2*j2^2*i3*j5 + 3*j2^3*i3*j5 + 5*i1^4*k2*i3*j5 + 6*i1^2*i2*k2*i3*j5 + 4*i2^2*k2*i3*j5 + 5*i2*j2*k2*i3*j5 +
            5*j2^2*k2*i3*j5 + 6*i1^2*k2^2*i3*j5 + 2*i2*k2^2*i3*j5 + 3*j2*k2^2*i3*j5 + 5*k2^3*i3*j5 + i1^4*l2*i3*j5 +
            i1^2*i2*l2*i3*j5 + 6*i2^2*l2*i3*j5 + 3*i1^2*j2*l2*i3*j5 + 3*i2*j2*l2*i3*j5 + j2^2*l2*i3*j5 + 5*i2*k2*l2*i3*j5
            + j2*k2*l2*i3*j5 + 2*i2*l2^2*i3*j5 + j2*l2^2*i3*j5 + 4*k2*l2^2*i3*j5 + 4*l2^3*i3*j5 + 6*i1^3*i3^2*j5 +
            4*i1*i2*i3^2*j5 + 3*i1*j2*i3^2*j5 + 4*i1*k2*i3^2*j5 + 2*i1*l2*i3^2*j5 + 5*i3^3*j5 + i1^6*j3*j5 +
            3*i1^4*i2*j3*j5 + i1^2*i2^2*j3*j5 + 2*i2^3*j3*j5 + 4*i1^2*j2^2*j3*j5 + i2*j2^2*j3*j5 + 4*j2^3*j3*j5 +
            2*i1^4*k2*j3*j5 + 4*i1^2*i2*k2*j3*j5 + 3*i2^2*k2*j3*j5 + 4*i1^2*j2*k2*j3*j5 + 5*i2*j2*k2*j3*j5 + j2^2*k2*j3*j5
            + 6*i2*k2^2*j3*j5 + 4*j2*k2^2*j3*j5 + 3*k2^3*j3*j5 + 6*i1^4*l2*j3*j5 + 2*i1^2*i2*l2*j3*j5 + 3*i2^2*l2*j3*j5 +
            4*i1^2*j2*l2*j3*j5 + 2*i2*j2*l2*j3*j5 + 6*j2^2*l2*j3*j5 + i1^2*k2*l2*j3*j5 + 3*i2*k2*l2*j3*j5 + j2*k2*l2*j3*j5
            + 5*k2^2*l2*j3*j5 + 2*i1^2*l2^2*j3*j5 + 6*i2*l2^2*j3*j5 + 6*j2*l2^2*j3*j5 + 6*k2*l2^2*j3*j5 + 2*i1^3*i3*j3*j5
            + 4*i1*i2*i3*j3*j5 + i1*j2*i3*j3*j5 + 4*i1*k2*i3*j3*j5 + 3*i1*l2*i3*j3*j5 + 4*i3^2*j3*j5 + 3*i1*i2*j3^2*j5 +
            3*i1*j2*j3^2*j5 + 4*i1*k2*j3^2*j5 + 4*i3*j3^2*j5 + 4*j3^3*j5 + 6*i1^6*k3*j5 + 6*i1^4*i2*k3*j5 + i2^3*k3*j5 +
            5*i1^4*j2*k3*j5 + 5*i2^2*j2*k3*j5 + i1^2*j2^2*k3*j5 + i2*j2^2*k3*j5 + 3*i1^4*k2*k3*j5 + 5*i1^2*i2*k2*k3*j5 +
            5*i2^2*k2*k3*j5 + i1^2*j2*k2*k3*j5 + 2*i2*j2*k2*k3*j5 + 6*j2^2*k2*k3*j5 + i1^2*k2^2*k3*j5 + 2*i2*k2^2*k3*j5 +
            3*i1^4*l2*k3*j5 + 6*i1^2*i2*l2*k3*j5 + 3*j2^2*l2*k3*j5 + 5*i1^2*k2*l2*k3*j5 + 4*i2*k2*l2*k3*j5 +
            3*j2*k2*l2*k3*j5 + 4*k2^2*l2*k3*j5 + i1^2*l2^2*k3*j5 + 4*i2*l2^2*k3*j5 + j2*l2^2*k3*j5 + 3*k2*l2^2*k3*j5 +
            3*i1^3*j3*k3*j5 + 6*i1*i2*j3*k3*j5 + 6*i1*j2*j3*k3*j5 + 6*i1*k2*j3*k3*j5 + i1*l2*j3*k3*j5 + 4*j3^2*k3*j5 +
            2*i1^3*k3^2*j5 + 5*i1*i2*k3^2*j5 + 4*i1*j2*k3^2*j5 + i1*k2*k3^2*j5 + 4*i1*l2*k3^2*j5 + 6*j3*k3^2*j5 +
            4*i1^6*l3*j5 + 4*i1^2*i2*j2*l3*j5 + i2^2*j2*l3*j5 + 3*i1^2*j2^2*l3*j5 + 4*i2*j2^2*l3*j5 + 2*j2^3*l3*j5 +
            6*i1^4*k2*l3*j5 + 3*i1^2*i2*k2*l3*j5 + 4*i2^2*k2*l3*j5 + 6*i1^2*j2*k2*l3*j5 + 4*i2*j2*k2*l3*j5 +
            2*j2^2*k2*l3*j5 + 5*i1^2*k2^2*l3*j5 + 2*i2*k2^2*l3*j5 + 4*j2*k2^2*l3*j5 + 4*i1^4*l2*l3*j5 + 2*i1^2*i2*l2*l3*j5
            + 5*i2^2*l2*l3*j5 + 5*i2*j2*l2*l3*j5 + 6*j2^2*l2*l3*j5 + 2*i1^2*k2*l2*l3*j5 + i2*k2*l2*l3*j5 +
            2*j2*k2*l2*l3*j5 + k2^2*l2*l3*j5 + 2*i2*l2^2*l3*j5 + 3*k2*l2^2*l3*j5 + 6*l2^3*l3*j5 + 5*i1^3*i3*l3*j5 +
            4*i1*i2*i3*l3*j5 + 4*i1*j2*i3*l3*j5 + 4*i1*k2*i3*l3*j5 + 5*i1*l2*i3*l3*j5 + 6*i3^2*l3*j5 + i1^3*l3^2*j5 +
            3*i1*i2*l3^2*j5 + 4*i1*j2*l3^2*j5 + 6*i1*k2*l3^2*j5 + 3*i1*l2*l3^2*j5 + 6*l3^3*j5 + 2*i1^5*i4*j5 +
            2*i1*i2^2*i4*j5 + 6*i1^3*j2*i4*j5 + 5*i1*i2*j2*i4*j5 + 3*i1*j2^2*i4*j5 + 6*i1^3*k2*i4*j5 + 2*i1*i2*k2*i4*j5 +
            5*i1*j2*k2*i4*j5 + 4*i1*k2^2*i4*j5 + i1^3*l2*i4*j5 + 3*i1*i2*l2*i4*j5 + 6*i1*j2*l2*i4*j5 + 6*i1*l2^2*i4*j5 +
            i1^2*i3*i4*j5 + 2*i2*i3*i4*j5 + l2*i3*i4*j5 + 2*i2*j3*i4*j5 + 4*j2*j3*i4*j5 + 2*i1^5*l4*j5 + 4*i1^3*i2*l4*j5 +
            2*i1*i2^2*l4*j5 + 5*i1*i2*j2*l4*j5 + i1*j2^2*l4*j5 + 3*i1^3*k2*l4*j5 + 6*i1*k2^2*l4*j5 + 6*i1^3*l2*l4*j5 +
            2*i1*i2*l2*l4*j5 + i1*j2*l2*l4*j5 + i1*k2*l2*l4*j5 + 4*i1*l2^2*l4*j5 + 5*i1^2*j3*l4*j5 + 6*i2*j3*l4*j5 +
            5*l2*j3*l4*j5 + i1^2*k3*l4*j5 + 5*i2*k3*l4*j5 + 4*j2*k3*l4*j5 + 3*k2*k3*l4*j5 + l2*k3*l4*j5 + 4*i1^2*l3*l4*j5
            + 3*i2*l3*l4*j5 + 2*j2*l3*l4*j5 + k2*l3*l4*j5 + 4*l2*l3*l4*j5 + 5*i1^4*j5^2 + 6*i1^2*i2*j5^2 + 5*i2^2*j5^2 +
            5*i1^2*j2*j5^2 + 5*i2*j2*j5^2 + j2^2*j5^2 + 2*i1^2*k2*j5^2 + i2*k2*j5^2 + 4*j2*k2*j5^2 + 6*k2^2*j5^2 +
            3*i1^2*l2*j5^2 + 6*i2*l2*j5^2 + 5*j2*l2*j5^2 + 2*k2*l2*j5^2 + 5*l2^2*j5^2 + 2*i1^9*k5 + i1^5*i2^2*k5 +
            6*i1^3*i2^3*k5 + 6*i1^7*j2*k5 + 3*i1^5*i2*j2*k5 + 3*i1^3*i2^2*j2*k5 + 3*i1^5*j2^2*k5 + 4*i1^3*i2*j2^2*k5 +
            i1^3*j2^3*k5 + 5*i1^5*i2*k2*k5 + 3*i1^3*i2^2*k2*k5 + 6*i1*i2^3*k2*k5 + 4*i1^5*j2*k2*k5 + i1^3*i2*j2*k2*k5 +
            3*i1*i2^2*j2*k2*k5 + 3*i1^3*j2^2*k2*k5 + 4*i1*i2*j2^2*k2*k5 + i1*j2^3*k2*k5 + i1^5*k2^2*k5 + 4*i1^3*i2*k2^2*k5
            + 6*i1*i2^2*k2^2*k5 + 5*i1^3*j2*k2^2*k5 + 5*i1*i2*j2*k2^2*k5 + 3*i1*j2^2*k2^2*k5 + 4*i1^3*k2^3*k5 +
            2*i1*i2*k2^3*k5 + 4*i1*j2*k2^3*k5 + 4*i1*k2^4*k5 + 4*i1^7*l2*k5 + 4*i1^5*i2*l2*k5 + 3*i1^3*i2^2*l2*k5 +
            6*i1*i2^3*l2*k5 + 4*i1^5*j2*l2*k5 + 5*i1^3*i2*j2*l2*k5 + 3*i1*i2^2*j2*l2*k5 + 6*i1^3*j2^2*l2*k5 +
            4*i1*i2*j2^2*l2*k5 + i1*j2^3*l2*k5 + 4*i1^5*k2*l2*k5 + 5*i1^3*i2*k2*l2*k5 + 3*i1*i2^2*k2*l2*k5 +
            6*i1^3*j2*k2*l2*k5 + 5*i1*i2*j2*k2*l2*k5 + 6*i1*j2^2*k2*l2*k5 + i1^3*k2^2*l2*k5 + 2*i1*k2^3*l2*k5 +
            5*i1^5*l2^2*k5 + 6*i1^3*i2*l2^2*k5 + i1^3*j2*l2^2*k5 + 5*i1*i2*j2*l2^2*k5 + 2*i1*j2^2*l2^2*k5 +
            6*i1^3*k2*l2^2*k5 + 2*i1*j2*k2*l2^2*k5 + i1*k2^2*l2^2*k5 + 2*i1^3*l2^3*k5 + 5*i1*i2*l2^3*k5 + 6*i1*j2*l2^3*k5
            + 5*i1*k2*l2^3*k5 + 4*i1*l2^4*k5 + 5*i1^2*i2^2*j3*k5 + 5*i1^4*j2*j3*k5 + 4*i1^2*i2*j2*j3*k5 +
            5*i1^2*j2^2*j3*k5 + 5*i1^4*k2*j3*k5 + 4*i1^2*i2*k2*j3*k5 + 6*i2^2*k2*j3*k5 + 6*i1^2*j2*k2*j3*k5 +
            2*i2*j2*k2*j3*k5 + 6*j2^2*k2*j3*k5 + 2*i1^2*k2^2*j3*k5 + 3*i2*k2^2*j3*k5 + 2*j2*k2^2*j3*k5 + 2*i1^4*l2*j3*k5 +
            4*i1^2*i2*l2*j3*k5 + 4*i2^2*l2*j3*k5 + 3*i1^2*j2*l2*j3*k5 + 6*i2*j2*l2*j3*k5 + 4*j2^2*l2*j3*k5 +
            4*i1^2*k2*l2*j3*k5 + i2*k2*l2*j3*k5 + 2*j2*k2*l2*j3*k5 + 5*k2^2*l2*j3*k5 + 4*i1^2*l2^2*j3*k5 + 6*i2*l2^2*j3*k5
            + 5*k2*l2^2*j3*k5 + 2*l2^3*j3*k5 + 2*i1^3*j3^2*k5 + i1*i2*j3^2*k5 + 6*i1*j2*j3^2*k5 + 5*i1*k2*j3^2*k5 +
            6*i1*l2*j3^2*k5 + j3^3*k5 + 2*i1^6*k3*k5 + 6*i1^4*i2*k3*k5 + i1^2*i2^2*k3*k5 + i1^4*j2*k3*k5 +
            5*i1^2*i2*j2*k3*k5 + i1^2*j2^2*k3*k5 + 2*i1^4*k2*k3*k5 + 2*i1^2*i2*k2*k3*k5 + 6*i2^2*k2*k3*k5 +
            6*i1^2*j2*k2*k3*k5 + 2*i2*j2*k2*k3*k5 + 6*j2^2*k2*k3*k5 + 2*i1^2*k2^2*k3*k5 + 4*i2*k2^2*k3*k5 +
            6*j2*k2^2*k3*k5 + 6*k2^3*k3*k5 + 4*i1^4*l2*k3*k5 + 3*i1^2*i2*l2*k3*k5 + 2*i2^2*l2*k3*k5 + 6*i1^2*j2*l2*k3*k5 +
            3*i2*j2*l2*k3*k5 + 2*j2^2*l2*k3*k5 + i1^2*k2*l2*k3*k5 + 6*i2*k2*l2*k3*k5 + 2*j2*k2*l2*k3*k5 + 4*k2^2*l2*k3*k5
            + i1^2*l2^2*k3*k5 + 5*i2*l2^2*k3*k5 + 4*j2*l2^2*k3*k5 + 2*k2*l2^2*k3*k5 + 4*l2^3*k3*k5 + 5*i1^3*k3^2*k5 +
            i1*i2*k3^2*k5 + 6*i1*j2*k3^2*k5 + 4*i1*k2*k3^2*k5 + 5*i1*l2*k3^2*k5 + 2*k3^3*k5 + 2*i1^6*l3*k5 +
            2*i1^4*i2*l3*k5 + 6*i1^2*i2^2*l3*k5 + i2^3*l3*k5 + i1^4*j2*l3*k5 + 5*i1^2*i2*j2*l3*k5 + 4*i2^2*j2*l3*k5 +
            3*i1^2*j2^2*l3*k5 + 3*i2*j2^2*l3*k5 + 6*j2^3*l3*k5 + 2*i1^4*k2*l3*k5 + 2*i1^2*j2*k2*l3*k5 + i2*j2*k2*l3*k5 +
            6*j2^2*k2*l3*k5 + i2*k2^2*l3*k5 + 5*k2^3*l3*k5 + 6*i1^4*l2*l3*k5 + 6*i1^2*i2*l2*l3*k5 + 2*i2^2*l2*l3*k5 +
            4*i1^2*j2*l2*l3*k5 + 5*i2*j2*l2*l3*k5 + 2*i2*k2*l2*l3*k5 + 4*j2*k2*l2*l3*k5 + 6*k2^2*l2*l3*k5 +
            5*i1^2*l2^2*l3*k5 + 2*i2*l2^2*l3*k5 + k2*l2^2*l3*k5 + l2^3*l3*k5 + 4*i1^3*j3*l3*k5 + 5*i1*i2*j3*l3*k5 +
            5*i1*j2*j3*l3*k5 + 4*i1*k2*j3*l3*k5 + 6*i1*l2*j3*l3*k5 + 6*j3^2*l3*k5 + 6*i1^5*l4*k5 + 3*i1*i2^2*l4*k5 +
            i1^3*j2*l4*k5 + i1*i2*j2*l4*k5 + 3*i1*j2^2*l4*k5 + 3*i1^3*k2*l4*k5 + 3*i1*i2*k2*l4*k5 + 3*i1*j2*k2*l4*k5 +
            2*i1^3*l2*l4*k5 + i1*j2*l2*l4*k5 + i1*k2*l2*l4*k5 + 6*i1*l2^2*l4*k5 + 4*i1^2*j3*l4*k5 + 6*i2*j3*l4*k5 +
            j2*j3*l4*k5 + i1^2*k3*l4*k5 + 2*i2*k3*l4*k5 + 5*j2*k3*l4*k5 + 4*k2*k3*l4*k5 + 3*l2*k3*l4*k5 + 3*i1^2*l3*l4*k5
            + 4*i2*l3*l4*k5 + 4*k2*l3*l4*k5 + i1*l4^2*k5 + 2*i1^9*l5 + 5*i1^7*i2*l5 + 2*i1^5*i2^2*l5 + 5*i1^3*i2^3*l5 +
            6*i1*i2^4*l5 + 4*i1^7*j2*l5 + 5*i1^5*i2*j2*l5 + 6*i1^3*i2^2*j2*l5 + 4*i1*i2^3*j2*l5 + 2*i1^5*j2^2*l5 +
            i1^3*i2*j2^2*l5 + i1*i2^2*j2^2*l5 + 2*i1^3*j2^3*l5 + 4*i1*i2*j2^3*l5 + 6*i1*j2^4*l5 + 6*i1^7*k2*l5 +
            4*i1^5*i2*k2*l5 + i1^3*i2^2*k2*l5 + 3*i1*i2^3*k2*l5 + i1^5*j2*k2*l5 + 4*i1^3*i2*j2*k2*l5 + 6*i1*i2^2*j2*k2*l5
            + 5*i1*j2^3*k2*l5 + 4*i1^5*k2^2*l5 + 6*i1^3*i2*k2^2*l5 + 5*i1*i2^2*k2^2*l5 + i1^3*j2*k2^2*l5 +
            5*i1*i2*j2*k2^2*l5 + 5*i1*j2^2*k2^2*l5 + 2*i1^3*k2^3*l5 + i1*i2*k2^3*l5 + 6*i1*j2*k2^3*l5 + 3*i1*k2^4*l5 +
            5*i1^7*l2*l5 + 3*i1^5*i2*l2*l5 + 5*i1^3*i2^2*l2*l5 + i1*i2^3*l2*l5 + 6*i1^5*j2*l2*l5 + 5*i1^3*i2*j2*l2*l5 +
            5*i1*i2^2*j2*l2*l5 + i1*i2*j2^2*l2*l5 + 4*i1^3*i2*k2*l2*l5 + i1*i2^2*k2*l2*l5 + i1*j2^2*k2*l2*l5 +
            i1*j2*k2^2*l2*l5 + i1*k2^3*l2*l5 + 6*i1^5*l2^2*l5 + 4*i1^3*i2*l2^2*l5 + 5*i1*i2^2*l2^2*l5 + 2*i1^3*j2*l2^2*l5
            + 3*i1*i2*j2*l2^2*l5 + i1*j2^2*l2^2*l5 + 5*i1^3*k2*l2^2*l5 + 2*i1*i2*k2*l2^2*l5 + 5*i1*j2*k2*l2^2*l5 +
            i1*j2*l2^3*l5 + i1*k2*l2^3*l5 + i1*l2^4*l5 + 5*i1^6*i3*l5 + i1^4*i2*i3*l5 + i1^2*i2^2*i3*l5 +
            5*i1^2*i2*j2*i3*l5 + i1^2*j2^2*i3*l5 + i1^4*k2*i3*l5 + 6*i1^2*i2*k2*i3*l5 + 3*i2^2*k2*i3*l5 + i1^2*j2*k2*i3*l5
            + i2*j2*k2*i3*l5 + 3*j2^2*k2*i3*l5 + i1^2*k2^2*i3*l5 + 4*i2*k2^2*i3*l5 + 6*j2*k2^2*i3*l5 + k2^3*i3*l5 +
            4*i1^4*l2*i3*l5 + 3*i1^2*i2*l2*i3*l5 + 4*i2^2*l2*i3*l5 + 3*i1^2*j2*l2*i3*l5 + 6*i2*j2*l2*i3*l5 +
            4*j2^2*l2*i3*l5 + 3*i1^2*k2*l2*i3*l5 + 2*j2*k2*l2*i3*l5 + 5*k2^2*l2*i3*l5 + j2*l2^2*i3*l5 + 6*k2*l2^2*i3*l5 +
            5*l2^3*i3*l5 + 5*i1^3*i3^2*l5 + 3*i1*k2*i3^2*l5 + 5*i1*l2*i3^2*l5 + 4*i1^4*i2*j3*l5 + 3*i1^2*i2^2*j3*l5 +
            5*i2^3*j3*l5 + 2*i1^4*j2*j3*l5 + 6*i1^2*i2*j2*j3*l5 + 6*i2^2*j2*j3*l5 + 5*i1^2*j2^2*j3*l5 + i2*j2^2*j3*l5 +
            2*j2^3*j3*l5 + 2*i1^4*k2*j3*l5 + 5*i2^2*k2*j3*l5 + i1^2*j2*k2*j3*l5 + i2*j2*k2*j3*l5 + 2*j2^2*k2*j3*l5 +
            5*i1^2*k2^2*j3*l5 + 2*i2*k2^2*j3*l5 + 5*j2*k2^2*j3*l5 + 2*k2^3*j3*l5 + 4*i1^2*i2*l2*j3*l5 + 2*i2^2*l2*j3*l5 +
            3*i2*j2*l2*j3*l5 + 2*j2^2*l2*j3*l5 + 6*i1^2*k2*l2*j3*l5 + 2*i2*k2*l2*j3*l5 + 5*j2*k2*l2*j3*l5 +
            4*k2^2*l2*j3*l5 + 6*i1^2*l2^2*j3*l5 + 2*j2*l2^2*j3*l5 + k2*l2^2*j3*l5 + 4*l2^3*j3*l5 + 5*i1^3*j3^2*l5 +
            5*i1*j2*j3^2*l5 + 6*i1*k2*j3^2*l5 + i1*l2*j3^2*l5 + 4*j3^3*l5 + 5*i1^6*k3*l5 + 2*i1^4*i2*k3*l5 +
            i1^2*i2^2*k3*l5 + 5*i2^3*k3*l5 + 3*i1^2*i2*j2*k3*l5 + 6*i2^2*j2*k3*l5 + 3*i1^2*j2^2*k3*l5 + i2*j2^2*k3*l5 +
            2*j2^3*k3*l5 + 2*i1^4*k2*k3*l5 + 2*i1^2*i2*k2*k3*l5 + 2*i2^2*k2*k3*l5 + i1^2*j2*k2*k3*l5 + 6*i2*j2*k2*k3*l5 +
            6*j2^2*k2*k3*l5 + 3*i1^2*k2^2*k3*l5 + 6*i2*k2^2*k3*l5 + 6*k2^3*k3*l5 + 3*i1^4*l2*k3*l5 + 5*i1^2*i2*l2*k3*l5 +
            3*i1^2*j2*l2*k3*l5 + 3*i2*j2*l2*k3*l5 + 4*j2^2*l2*k3*l5 + 5*i1^2*k2*l2*k3*l5 + 5*i2*k2*l2*k3*l5 +
            j2*k2*l2*k3*l5 + 5*k2^2*l2*k3*l5 + 5*i1^2*l2^2*k3*l5 + 3*i2*l2^2*k3*l5 + 3*k2*l2^2*k3*l5 + 3*l2^3*k3*l5 +
            6*i1^3*j3*k3*l5 + 2*i1*j2*j3*k3*l5 + i1*l2*j3*k3*l5 + 5*j3^2*k3*l5 + 2*i1*i2*k3^2*l5 + 2*i1*j2*k3^2*l5 +
            2*j3*k3^2*l5 + 2*k3^3*l5 + 5*i1^4*i2*l3*l5 + 3*i1^2*i2^2*l3*l5 + 4*i1^4*j2*l3*l5 + 5*i1^2*i2*j2*l3*l5 +
            2*i2^2*j2*l3*l5 + 5*i1^2*j2^2*l3*l5 + 3*i2*j2^2*l3*l5 + 2*j2^3*l3*l5 + 4*i1^4*k2*l3*l5 + 4*i1^2*i2*k2*l3*l5 +
            6*i2^2*k2*l3*l5 + 5*i1^2*j2*k2*l3*l5 + 6*i2*j2*k2*l3*l5 + 4*j2^2*k2*l3*l5 + 2*i1^2*k2^2*l3*l5 +
            4*i2*k2^2*l3*l5 + 5*j2*k2^2*l3*l5 + 2*k2^3*l3*l5 + i2^2*l2*l3*l5 + 3*i1^2*j2*l2*l3*l5 + 4*i2*j2*l2*l3*l5 +
            4*j2^2*l2*l3*l5 + 2*i1^2*k2*l2*l3*l5 + 3*i2*k2*l2*l3*l5 + 2*j2*k2*l2*l3*l5 + 4*k2^2*l2*l3*l5 + i1^2*l2^2*l3*l5
            + 2*i2*l2^2*l3*l5 + 2*k2*l2^2*l3*l5 + 2*l2^3*l3*l5 + 4*i1^3*i3*l3*l5 + 6*i1*i2*i3*l3*l5 + i1*j2*i3*l3*l5 +
            i1*k2*i3*l3*l5 + 6*i1*l2*i3*l3*l5 + 5*i3^2*l3*l5 + 4*i1^3*j3*l3*l5 + 2*i1*i2*j3*l3*l5 + 6*i1*j2*j3*l3*l5 +
            4*i1*k2*j3*l3*l5 + 6*i1*l2*j3*l3*l5 + j3^2*l3*l5 + 5*i1^3*l3^2*l5 + 2*i1*i2*l3^2*l5 + 3*i1*k2*l3^2*l5 +
            6*i1*l2*l3^2*l5 + 6*i3*l3^2*l5 + l3^3*l5 + 6*i1^5*k4*l5 + 3*i1^3*i2*k4*l5 + 4*i1^3*j2*k4*l5 + 6*i1^3*k2*k4*l5
            + 2*i1*i2*k2*k4*l5 + 5*i1*j2*k2*k4*l5 + i1^3*l2*k4*l5 + 2*i1*i2*l2*k4*l5 + 5*i1*j2*l2*k4*l5 + 4*i1*k2*l2*k4*l5
            + 4*i1*l2^2*k4*l5 + 5*i1^2*k3*k4*l5 + 3*k2*k3*k4*l5 + 6*l2*k3*k4*l5 + 3*i1^5*m4*l5 + 5*i1^3*i2*m4*l5 +
            3*i1^3*j2*m4*l5 + 4*i1*i2*j2*m4*l5 + 3*i1*j2^2*m4*l5 + i1*i2*k2*m4*l5 + 6*i1*j2*k2*m4*l5 + i1*k2^2*m4*l5 +
            2*i1^3*l2*m4*l5 + 2*i1*i2*l2*m4*l5 + 4*i1*j2*l2*m4*l5 + 2*i1*l2^2*m4*l5 + 3*i1^2*i3*m4*l5 + 6*i2*i3*m4*l5 +
            j2*i3*m4*l5 + 6*k2*i3*m4*l5 + 6*l2*i3*m4*l5 + 6*i1^2*j3*m4*l5 + 4*i2*j3*m4*l5 + 5*k2*j3*m4*l5 + 3*l2*j3*m4*l5
            + i1^2*l3*m4*l5 + 2*i2*l3*m4*l5 + 5*k2*l3*m4*l5 + l2*l3*m4*l5 + 4*i1^2*i2*l5^2 + 4*i2^2*l5^2 + 4*i1^2*j2*l5^2
            + 6*i2*j2*l5^2 + 4*j2^2*l5^2 + i1^2*k2*l5^2 + 4*i2*k2*l5^2 + 2*j2*k2*l5^2 + 4*k2^2*l5^2 + 2*i1^2*l2*l5^2 +
            6*i2*l2*l5^2 + 6*j2*l2*l5^2 + 6*k2*l2*l5^2 + i1*i3*l5^2 + 3*i1*j3*l5^2 + 5*i1*k3*l5^2 + 5*k4*l5^2 + 4*m4*l5^2
            + 4*i1^6*i2*i6 + 3*i1^6*j2*i6 + 2*i1^6*k2*i6 + 6*i1^4*i2*k2*i6 + i1^4*j2*k2*i6 + i1^4*k2^2*i6 +
            3*i1^2*i2*k2^2*i6 + 4*i1^2*j2*k2^2*i6 + 3*i1^2*k2^3*i6 + 4*i2*k2^3*i6 + 3*j2*k2^3*i6 + 2*i1^6*l2*i6 +
            2*i1^4*i2*l2*i6 + 5*i1^4*j2*l2*i6 + i1^4*k2*l2*i6 + 4*i1^2*i2*k2*l2*i6 + 3*i1^2*j2*k2*l2*i6 + i1^2*k2^2*l2*i6
            + 5*i2*k2^2*l2*i6 + 2*j2*k2^2*l2*i6 + 6*k2^3*l2*i6 + 3*i1^4*l2^2*i6 + 2*i1^2*i2*l2^2*i6 + 5*i1^2*j2*l2^2*i6 +
            6*i1^2*k2*l2^2*i6 + 3*i2*k2*l2^2*i6 + 4*j2*k2*l2^2*i6 + 3*k2^2*l2^2*i6 + 3*i1^2*l2^3*i6 + i2*l2^3*i6 +
            2*j2*l2^3*i6 + k2*l2^3*i6 + l2^4*i6 + 3*i1^5*j3*i6 + i1^3*k2*j3*i6 + i1*k2^2*j3*i6 + i1^3*l2*j3*i6 +
            5*i1*k2*l2*j3*i6 + 4*i1*l2^2*j3*i6 + 3*i1^5*k3*i6 + 5*i1^3*k2*k3*i6 + 2*i1*k2^2*k3*i6 + 3*i1^3*l2*k3*i6 +
            5*i1*k2*l2*k3*i6 + 6*i1*l2^2*k3*i6 + 2*i1^5*l3*i6 + i1^3*i2*l3*i6 + 6*i1^3*j2*l3*i6 + 4*i1^3*k2*l3*i6 +
            4*i1*i2*k2*l3*i6 + 3*i1*j2*k2*l3*i6 + 5*i1*k2^2*l3*i6 + i1^3*l2*l3*i6 + 6*i1*i2*l2*l3*i6 + i1*j2*l2*l3*i6 +
            5*i1*k2*l2*l3*i6 + 3*i1*l2^2*l3*i6 + i1^2*j3*l3*i6 + 3*k2*j3*l3*i6 + 2*l2*j3*l3*i6 + 6*i1^4*l4*i6 +
            5*i1^2*k2*l4*i6 + k2^2*l4*i6 + 3*i1^2*l2*l4*i6 + k2*l2*l4*i6 + 2*l2^2*l4*i6 + 3*i1*l3*l4*i6 + 3*l4^2*i6 +
            3*i1^7*i7;
        Append(~JI, J14); Append(~Wght, 14);
    end if;
    if degmax le 14 then return JI, Wght; end if;

    /* Degree 15 */
    if not PrimaryOnly and degmin le 15 then
        J15:=
            i1^15 + 2*i1^13*i2 + 3*i1^11*i2^2 + i1^9*i2^3 + 3*i1^7*i2^4 + 3*i1^5*i2^5 + 6*i1^3*i2^6 + 2*i1^13*j2 +
            6*i1^11*i2*j2 + i1^9*i2^2*j2 + 4*i1^7*i2^3*j2 + 2*i1^5*i2^4*j2 + 6*i1^3*i2^5*j2 + i1^11*j2^2 +
            4*i1^9*i2*j2^2 + 2*i1^7*i2^2*j2^2 + 4*i1^5*i2^3*j2^2 + 6*i1^3*i2^4*j2^2 + 2*i1^9*j2^3 + 2*i1^5*i2^2*j2^3 +
            6*i1^3*i2^3*j2^3 + 5*i1^7*j2^4 + 3*i1^5*i2*j2^4 + 6*i1^3*i2^2*j2^4 + 6*i1^3*i2*j2^5 + 6*i1^3*j2^6 +
            4*i1^13*k2 + 4*i1^11*i2*k2 + 5*i1^9*i2^2*k2 + 5*i1^7*i2^3*k2 + 4*i1^5*i2^4*k2 + 6*i1^3*i2^5*k2 +
            i1*i2^6*k2 + 2*i1^11*j2*k2 + 5*i1^7*i2^2*j2*k2 + 4*i1^5*i2^3*j2*k2 + 5*i1^3*i2^4*j2*k2 + i1*i2^5*j2*k2 +
            2*i1^9*j2^2*k2 + 5*i1^7*i2*j2^2*k2 + 4*i1^3*i2^3*j2^2*k2 + i1*i2^4*j2^2*k2 + 3*i1^7*j2^3*k2 +
            3*i1^3*i2^2*j2^3*k2 + i1*i2^3*j2^3*k2 + 6*i1^5*j2^4*k2 + 2*i1^3*i2*j2^4*k2 + i1*i2^2*j2^4*k2 +
            i1^3*j2^5*k2 + i1*i2*j2^5*k2 + i1*j2^6*k2 + 2*i1^9*i2*k2^2 + 5*i1^7*i2^2*k2^2 + 2*i1^5*i2^3*k2^2 +
            4*i1^3*i2^4*k2^2 + 3*i1*i2^5*k2^2 + 6*i1^7*i2*j2*k2^2 + 5*i1^5*i2^2*j2*k2^2 + 5*i1*i2^4*j2*k2^2 +
            2*i1^7*j2^2*k2^2 + 5*i1^5*i2*j2^2*k2^2 + i1^3*i2^2*j2^2*k2^2 + 6*i1*i2^3*j2^2*k2^2 + 2*i1^5*j2^3*k2^2 +
            3*i1^3*i2*j2^3*k2^2 + 6*i1*i2^2*j2^3*k2^2 + 6*i1^3*j2^4*k2^2 + 5*i1*i2*j2^4*k2^2 + 3*i1*j2^5*k2^2 +
            4*i1^9*k2^3 + 6*i1^7*i2*k2^3 + i1^5*i2^2*k2^3 + i1*i2^4*k2^3 + 3*i1^7*j2*k2^3 + i1^3*i2^2*j2*k2^3 +
            6*i1*i2^3*j2*k2^3 + 3*i1^5*j2^2*k2^3 + 5*i1^3*i2*j2^2*k2^3 + i1*i2^2*j2^2*k2^3 + i1^3*j2^3*k2^3 +
            4*i1*i2*j2^3*k2^3 + 2*i1*j2^4*k2^3 + 2*i1^7*k2^4 + 3*i1^5*i2*k2^4 + 6*i1^3*i2^2*k2^4 + 2*i1*i2^3*k2^4 +
            5*i1^3*i2*j2*k2^4 + 2*i1^3*j2^2*k2^4 + 4*i1*i2*j2^2*k2^4 + 4*i1*j2^3*k2^4 + 3*i1^5*k2^5 + 6*i1^3*i2*k2^5 +
            i1*i2^2*k2^5 + i1^3*j2*k2^5 + 2*i1*i2*j2*k2^5 + 5*i1*j2^2*k2^5 + 4*i1^3*k2^6 + 4*i1*i2*k2^6 + 2*i1*j2*k2^6
            + i1*k2^7 + 4*i1^13*l2 + i1^11*i2*l2 + 5*i1^9*i2^2*l2 + 3*i1^7*i2^3*l2 + i1^5*i2^4*l2 + 2*i1^3*i2^5*l2 +
            2*i1*i2^6*l2 + i1^11*j2*l2 + 6*i1^7*i2^2*j2*l2 + 4*i1^3*i2^4*j2*l2 + 2*i1*i2^5*j2*l2 + 6*i1^9*j2^2*l2 +
            4*i1^7*i2*j2^2*l2 + 6*i1^3*i2^3*j2^2*l2 + 2*i1*i2^4*j2^2*l2 + 3*i1^5*i2*j2^3*l2 + i1^3*i2^2*j2^3*l2 +
            2*i1*i2^3*j2^3*l2 + 3*i1^5*j2^4*l2 + 3*i1^3*i2*j2^4*l2 + 2*i1*i2^2*j2^4*l2 + 5*i1^3*j2^5*l2 +
            2*i1*i2*j2^5*l2 + 2*i1*j2^6*l2 + 5*i1^11*k2*l2 + 3*i1^7*i2^2*k2*l2 + 5*i1^5*i2^3*k2*l2 + 4*i1^3*i2^4*k2*l2
            + 4*i1*i2^5*k2*l2 + i1^7*i2*j2*k2*l2 + 2*i1^5*i2^2*j2*k2*l2 + i1*i2^4*j2*k2*l2 + i1^7*j2^2*k2*l2 +
            4*i1^5*i2*j2^2*k2*l2 + 6*i1^3*i2^2*j2^2*k2*l2 + 5*i1*i2^3*j2^2*k2*l2 + 2*i1^5*j2^3*k2*l2 +
            2*i1*i2^2*j2^3*k2*l2 + 4*i1^3*j2^4*k2*l2 + 6*i1*i2*j2^4*k2*l2 + 3*i1*j2^5*k2*l2 + i1^9*k2^2*l2 +
            5*i1^7*i2*k2^2*l2 + 6*i1^5*i2^2*k2^2*l2 + 5*i1^3*i2^3*k2^2*l2 + i1*i2^4*k2^2*l2 + 5*i1^7*j2*k2^2*l2 +
            6*i1^5*i2*j2*k2^2*l2 + 3*i1^3*i2^2*j2*k2^2*l2 + 5*i1*i2^3*j2*k2^2*l2 + 6*i1^5*j2^2*k2^2*l2 +
            4*i1^3*i2*j2^2*k2^2*l2 + 2*i1*i2^2*j2^2*k2^2*l2 + 5*i1*i2*j2^3*k2^2*l2 + i1*j2^4*k2^2*l2 + 4*i1^7*k2^3*l2
            + i1^5*i2*k2^3*l2 + 3*i1^3*i2^2*k2^3*l2 + 4*i1*i2^3*k2^3*l2 + 4*i1^5*j2*k2^3*l2 + 4*i1*i2^2*j2*k2^3*l2 +
            5*i1^3*j2^2*k2^3*l2 + 6*i1*i2*j2^2*k2^3*l2 + 6*i1*j2^3*k2^3*l2 + 2*i1^5*k2^4*l2 + i1^3*i2*k2^4*l2 +
            2*i1*i2^2*k2^4*l2 + 6*i1*i2*j2*k2^4*l2 + i1*j2^2*k2^4*l2 + 3*i1*i2*k2^5*l2 + i1*j2*k2^5*l2 + 4*i1*k2^6*l2
            + 2*i1^11*l2^2 + 6*i1^9*i2*l2^2 + i1^7*i2^2*l2^2 + 4*i1^5*i2^3*l2^2 + 5*i1^3*i2^4*l2^2 + i1*i2^5*l2^2 +
            2*i1^9*j2*l2^2 + 2*i1^7*i2*j2*l2^2 + 6*i1^5*i2^2*j2*l2^2 + 2*i1^3*i2^3*j2*l2^2 + 3*i1*i2^4*j2*l2^2 +
            4*i1^7*j2^2*l2^2 + 6*i1^5*i2*j2^2*l2^2 + i1^3*i2^2*j2^2*l2^2 + 6*i1*i2^3*j2^2*l2^2 + i1^5*j2^3*l2^2 +
            3*i1*i2^2*j2^3*l2^2 + 6*i1^3*j2^4*l2^2 + i1*i2*j2^4*l2^2 + 4*i1^9*k2*l2^2 + 2*i1^7*i2*k2*l2^2 +
            6*i1^5*i2^2*k2*l2^2 + i1^3*i2^3*k2*l2^2 + i1*i2^4*k2*l2^2 + 6*i1^7*j2*k2*l2^2 + 2*i1^3*i2^2*j2*k2*l2^2 +
            6*i1*i2^3*j2*k2*l2^2 + 6*i1^5*j2^2*k2*l2^2 + 6*i1^3*i2*j2^2*k2*l2^2 + 3*i1*i2^2*j2^2*k2*l2^2 +
            6*i1^3*j2^3*k2*l2^2 + 4*i1*j2^4*k2*l2^2 + 5*i1^5*i2*k2^2*l2^2 + 4*i1*i2^3*k2^2*l2^2 +
            4*i1^3*i2*j2*k2^2*l2^2 + 5*i1*i2^2*j2*k2^2*l2^2 + i1^3*j2^2*k2^2*l2^2 + i1*i2*j2^2*k2^2*l2^2 +
            4*i1*j2^3*k2^2*l2^2 + 6*i1^5*k2^3*l2^2 + 4*i1^3*i2*k2^3*l2^2 + 4*i1*i2^2*k2^3*l2^2 + 3*i1^3*j2*k2^3*l2^2 +
            6*i1*i2*j2*k2^3*l2^2 + 2*i1*j2^2*k2^3*l2^2 + i1^3*k2^4*l2^2 + 2*i1*i2*k2^4*l2^2 + 4*i1*j2*k2^4*l2^2 +
            4*i1*k2^5*l2^2 + 6*i1^9*l2^3 + 4*i1^7*i2*l2^3 + i1^5*i2^2*l2^3 + 4*i1^3*i2^3*l2^3 + 2*i1^7*j2*l2^3 +
            i1^5*i2*j2*l2^3 + 4*i1^3*i2^2*j2*l2^3 + 6*i1*i2^3*j2*l2^3 + 6*i1^5*j2^2*l2^3 + i1^3*i2*j2^2*l2^3 +
            3*i1*i2^2*j2^2*l2^3 + 5*i1^3*j2^3*l2^3 + 4*i1*i2*j2^3*l2^3 + i1*j2^4*l2^3 + 3*i1^7*k2*l2^3 +
            3*i1^5*i2*k2*l2^3 + 2*i1^3*i2^2*k2*l2^3 + 6*i1^5*j2*k2*l2^3 + 5*i1^3*i2*j2*k2*l2^3 + 6*i1*i2^2*j2*k2*l2^3
            + 3*i1^3*j2^2*k2*l2^3 + 2*i1*j2^3*k2*l2^3 + 4*i1^3*i2*k2^2*l2^3 + 2*i1*i2^2*k2^2*l2^3 + i1^3*j2*k2^2*l2^3
            + 5*i1*i2*j2*k2^2*l2^3 + 4*i1*j2^2*k2^2*l2^3 + i1^3*k2^3*l2^3 + 3*i1*i2*k2^3*l2^3 + 6*i1*j2*k2^3*l2^3 +
            i1*k2^4*l2^3 + 3*i1^7*l2^4 + i1^3*i2^2*l2^4 + 3*i1*i2^3*l2^4 + 4*i1^5*j2*l2^4 + 3*i1^3*i2*j2*l2^4 +
            3*i1*j2^3*l2^4 + 6*i1^5*k2*l2^4 + 6*i1^3*i2*k2*l2^4 + i1^3*j2*k2*l2^4 + i1*i2*j2*k2*l2^4 +
            3*i1*j2^2*k2*l2^4 + 3*i1^3*k2^2*l2^4 + 6*i1*i2*k2^2*l2^4 + 2*i1*j2*k2^2*l2^4 + i1*k2^3*l2^4 + 2*i1^5*l2^5
            + 6*i1^3*i2*l2^5 + 5*i1*i2^2*l2^5 + 6*i1^3*j2*l2^5 + 4*i1*i2*j2*l2^5 + 6*i1*j2^2*l2^5 + 6*i1^3*k2*l2^5 +
            5*i1*i2*k2*l2^5 + i1*j2*k2*l2^5 + 6*i1*k2^2*l2^5 + 3*i1^3*l2^6 + 2*i1*i2*l2^6 + 4*i1*j2*l2^6 +
            3*i1*k2*l2^6 + 2*i1*l2^7 + 6*i1^12*i3 + 4*i1^10*i2*i3 + 6*i1^8*i2^2*i3 + 3*i1^6*i2^3*i3 + i1^4*i2^4*i3 +
            5*i1^10*j2*i3 + i1^8*i2*j2*i3 + 3*i1^6*i2^2*j2*i3 + 3*i1^4*i2^3*j2*i3 + i1^8*j2^2*i3 + 6*i1^6*i2*j2^2*i3 +
            6*i1^4*i2^2*j2^2*i3 + 2*i1^6*j2^3*i3 + 3*i1^4*i2*j2^3*i3 + i1^4*j2^4*i3 + 2*i1^10*k2*i3 + 2*i1^8*i2*k2*i3
            + 2*i1^6*i2^2*k2*i3 + 3*i1^4*i2^3*k2*i3 + 4*i1^2*i2^4*k2*i3 + i1^8*j2*k2*i3 + 3*i1^4*i2^2*j2*k2*i3 +
            5*i1^2*i2^3*j2*k2*i3 + 5*i1^6*j2^2*k2*i3 + 6*i1^4*i2*j2^2*k2*i3 + 3*i1^2*i2^2*j2^2*k2*i3 +
            2*i1^4*j2^3*k2*i3 + 5*i1^2*i2*j2^3*k2*i3 + 4*i1^2*j2^4*k2*i3 + i1^8*k2^2*i3 + i1^4*i2^2*k2^2*i3 +
            2*i1^2*i2^3*k2^2*i3 + 5*i2^4*k2^2*i3 + 4*i1^6*j2*k2^2*i3 + 5*i1^4*i2*j2*k2^2*i3 + 3*i1^2*i2^2*j2*k2^2*i3 +
            i2^3*j2*k2^2*i3 + i1^4*j2^2*k2^2*i3 + 2*i1^2*i2*j2^2*k2^2*i3 + 2*i2^2*j2^2*k2^2*i3 + i2*j2^3*k2^2*i3 +
            5*j2^4*k2^2*i3 + 3*i1^6*k2^3*i3 + 2*i1^4*i2*k2^3*i3 + i1^2*i2^2*k2^3*i3 + 5*i2^3*k2^3*i3 + i1^4*j2*k2^3*i3
            + i1^2*i2*j2*k2^3*i3 + 5*i2^2*j2*k2^3*i3 + 3*i1^2*j2^2*k2^3*i3 + 3*i2*j2^2*k2^3*i3 + j2^3*k2^3*i3 +
            5*i1^2*i2*k2^4*i3 + i2^2*k2^4*i3 + 2*i1^2*j2*k2^4*i3 + 4*i2*j2*k2^4*i3 + 4*j2*k2^5*i3 + 3*k2^6*i3 +
            i1^10*l2*i3 + 5*i1^8*i2*l2*i3 + 2*i1^6*i2^2*l2*i3 + 4*i1^4*i2^3*l2*i3 + 2*i1^2*i2^4*l2*i3 +
            3*i1^8*j2*l2*i3 + 5*i1^6*i2*j2*l2*i3 + 6*i1^2*i2^3*j2*l2*i3 + 3*i1^6*j2^2*l2*i3 + 2*i1^4*i2*j2^2*l2*i3 +
            5*i1^2*i2^2*j2^2*l2*i3 + i1^4*j2^3*l2*i3 + 6*i1^2*i2*j2^3*l2*i3 + 2*i1^2*j2^4*l2*i3 + 6*i1^8*k2*l2*i3 +
            i1^6*i2*k2*l2*i3 + 5*i1^4*i2^2*k2*l2*i3 + 3*i2^4*k2*l2*i3 + i1^6*j2*k2*l2*i3 + 4*i1^4*i2*j2*k2*l2*i3 +
            6*i1^2*i2^2*j2*k2*l2*i3 + 2*i2^3*j2*k2*l2*i3 + i1^4*j2^2*k2*l2*i3 + 2*i1^2*i2*j2^2*k2*l2*i3 +
            4*i2^2*j2^2*k2*l2*i3 + 6*i1^2*j2^3*k2*l2*i3 + 2*i2*j2^3*k2*l2*i3 + 3*j2^4*k2*l2*i3 + 6*i1^4*i2*k2^2*l2*i3
            + 4*i1^2*i2^2*k2^2*l2*i3 + 4*i2^3*k2^2*l2*i3 + 3*i1^4*j2*k2^2*l2*i3 + 4*i1^2*i2*j2*k2^2*l2*i3 +
            3*i2^2*j2*k2^2*l2*i3 + 2*i1^2*j2^2*k2^2*l2*i3 + 3*i2*j2^2*k2^2*l2*i3 + 4*j2^3*k2^2*l2*i3 +
            3*i1^4*k2^3*l2*i3 + 2*i1^2*j2*k2^3*l2*i3 + 6*i2*j2*k2^3*l2*i3 + 6*i1^2*k2^4*l2*i3 + 4*i2*k2^4*l2*i3 +
            2*j2*k2^4*l2*i3 + 6*i1^8*l2^2*i3 + 4*i1^6*i2*l2^2*i3 + 5*i1^4*i2^2*l2^2*i3 + i1^2*i2^3*l2^2*i3 +
            6*i2^4*l2^2*i3 + 2*i1^6*j2*l2^2*i3 + 6*i1^4*i2*j2*l2^2*i3 + 2*i1^2*i2^2*j2*l2^2*i3 + 4*i2^3*j2*l2^2*i3 +
            5*i1^4*j2^2*l2^2*i3 + i2^2*j2^2*l2^2*i3 + 4*i1^2*j2^3*l2^2*i3 + 4*i2*j2^3*l2^2*i3 + 6*j2^4*l2^2*i3 +
            3*i1^6*k2*l2^2*i3 + 2*i1^2*i2^2*k2*l2^2*i3 + 4*i2^3*k2*l2^2*i3 + i1^4*j2*k2*l2^2*i3 +
            5*i1^2*i2*j2*k2*l2^2*i3 + 2*i2^2*j2*k2*l2^2*i3 + 5*i1^2*j2^2*k2*l2^2*i3 + 5*i2*j2^2*k2*l2^2*i3 +
            3*j2^3*k2*l2^2*i3 + 4*i1^4*k2^2*l2^2*i3 + 3*i1^2*i2*k2^2*l2^2*i3 + 3*i2^2*k2^2*l2^2*i3 +
            i2*j2*k2^2*l2^2*i3 + 5*j2^2*k2^2*l2^2*i3 + 6*i1^2*k2^3*l2^2*i3 + 6*j2*k2^3*l2^2*i3 + k2^4*l2^2*i3 +
            2*i1^6*l2^3*i3 + 5*i1^2*i2^2*l2^3*i3 + i1^2*i2*j2*l2^3*i3 + i1^2*j2^2*l2^3*i3 + 5*i1^4*k2*l2^3*i3 +
            5*i1^2*i2*k2*l2^3*i3 + 6*i2^2*k2*l2^3*i3 + 2*i1^2*j2*k2*l2^3*i3 + i2*j2*k2*l2^3*i3 + 3*j2^2*k2*l2^3*i3 +
            2*i1^2*k2^2*l2^3*i3 + 2*i2*k2^2*l2^3*i3 + k2^3*l2^3*i3 + i1^4*l2^4*i3 + 3*i1^2*i2*l2^4*i3 + j2^2*l2^4*i3 +
            5*i1^2*k2*l2^4*i3 + 4*i2*k2*l2^4*i3 + 4*j2*k2*l2^4*i3 + 3*k2^2*l2^4*i3 + 2*i1^2*l2^5*i3 + 2*i2*l2^5*i3 +
            3*k2*l2^5*i3 + 2*l2^6*i3 + 2*i1^9*i3^2 + 2*i1^7*i2*i3^2 + 3*i1^5*i2^2*i3^2 + 2*i1^7*j2*i3^2 +
            i1^5*i2*j2*i3^2 + 3*i1^5*j2^2*i3^2 + 5*i1^7*k2*i3^2 + 4*i1^5*i2*k2*i3^2 + 5*i1^5*j2*k2*i3^2 +
            4*i1^5*k2^2*i3^2 + 5*i1^3*i2*k2^2*i3^2 + i1*i2^2*k2^2*i3^2 + 5*i1*i2*j2*k2^2*i3^2 + i1*j2^2*k2^2*i3^2 +
            3*i1^3*k2^3*i3^2 + 3*i1*i2*k2^3*i3^2 + 4*i1*j2*k2^3*i3^2 + 4*i1*k2^4*i3^2 + 2*i1^7*l2*i3^2 +
            i1^5*i2*l2*i3^2 + 4*i1^3*i2^2*l2*i3^2 + 6*i1^3*i2*j2*l2*i3^2 + 4*i1^3*j2^2*l2*i3^2 + 3*i1^5*k2*l2*i3^2 +
            3*i1^3*i2*k2*l2*i3^2 + i1*i2^2*k2*l2*i3^2 + 5*i1*i2*j2*k2*l2*i3^2 + i1*j2^2*k2*l2*i3^2 +
            4*i1^3*k2^2*l2*i3^2 + 6*i1*i2*k2^2*l2*i3^2 + 4*i1*j2*k2^2*l2*i3^2 + 2*i1*k2^3*l2*i3^2 + 4*i1^5*l2^2*i3^2 +
            i1*i2^2*l2^2*i3^2 + 5*i1*i2*j2*l2^2*i3^2 + i1*j2^2*l2^2*i3^2 + 6*i1^3*k2*l2^2*i3^2 + 3*i1*i2*k2*l2^2*i3^2
            + 2*i1*j2*k2*l2^2*i3^2 + 6*i1*k2^2*l2^2*i3^2 + 2*i1^3*l2^3*i3^2 + 6*i1*i2*l2^3*i3^2 + 4*i1*k2*l2^3*i3^2 +
            6*i1*l2^4*i3^2 + i1^6*i3^3 + i1^4*k2*i3^3 + 2*i1^2*k2^2*i3^3 + 3*i1^4*l2*i3^3 + 3*k2^2*l2*i3^3 +
            i1^2*l2^2*i3^3 + 3*k2*l2^2*i3^3 + 6*l2^3*i3^3 + 2*i1^12*j3 + 3*i1^10*i2*j3 + 4*i1^8*i2^2*j3 +
            5*i1^6*i2^3*j3 + 6*i1^4*i2^4*j3 + i1^2*i2^5*j3 + 5*i1^10*j2*j3 + 6*i1^8*i2*j2*j3 + 4*i1^6*i2^2*j2*j3 +
            5*i1^4*i2^3*j2*j3 + 2*i1^2*i2^4*j2*j3 + 5*i1^8*j2^2*j3 + 5*i1^4*i2^2*j2^2*j3 + 3*i1^2*i2^3*j2^2*j3 +
            5*i1^6*j2^3*j3 + 4*i1^2*i2^2*j2^3*j3 + 5*i1^4*j2^4*j3 + 5*i1^2*i2*j2^4*j3 + 6*i1^2*j2^5*j3 + 3*i1^10*k2*j3
            + 5*i1^8*i2*k2*j3 + 5*i1^4*i2^3*k2*j3 + i2^5*k2*j3 + i1^8*j2*k2*j3 + 6*i1^6*i2*j2*k2*j3 +
            4*i1^4*i2^2*j2*k2*j3 + 3*i1^2*i2^3*j2*k2*j3 + 2*i2^4*j2*k2*j3 + i1^6*j2^2*k2*j3 + 3*i1^4*i2*j2^2*k2*j3 +
            5*i1^2*i2^2*j2^2*k2*j3 + 3*i2^3*j2^2*k2*j3 + 2*i1^4*j2^3*k2*j3 + 2*i1^2*i2*j2^3*k2*j3 + 4*i2^2*j2^3*k2*j3
            + 4*i1^2*j2^4*k2*j3 + 5*i2*j2^4*k2*j3 + 6*j2^5*k2*j3 + 4*i1^8*k2^2*j3 + 2*i1^6*i2*k2^2*j3 +
            2*i1^4*i2^2*k2^2*j3 + 6*i2^4*k2^2*j3 + 4*i1^6*j2*k2^2*j3 + 2*i1^4*i2*j2*k2^2*j3 + 5*i1^2*i2^2*j2*k2^2*j3 +
            3*i1^2*i2*j2^2*k2^2*j3 + 6*i2^2*j2^2*k2^2*j3 + 6*i1^2*j2^3*k2^2*j3 + 6*i2*j2^3*k2^2*j3 + 3*j2^4*k2^2*j3 +
            4*i1^6*k2^3*j3 + 2*i1^4*i2*k2^3*j3 + i1^2*i2^2*k2^3*j3 + 2*i2^3*k2^3*j3 + 2*i1^4*j2*k2^3*j3 +
            6*i1^2*i2*j2*k2^3*j3 + 4*i2^2*j2*k2^3*j3 + 2*i1^2*j2^2*k2^3*j3 + 4*i2*j2^2*k2^3*j3 + 4*j2^3*k2^3*j3 +
            6*i1^4*k2^4*j3 + 4*i2^2*k2^4*j3 + 6*i1^2*j2*k2^4*j3 + 5*i2*j2*k2^4*j3 + 4*j2^2*k2^4*j3 + 5*i1^2*k2^5*j3 +
            3*i2*k2^5*j3 + 5*j2*k2^5*j3 + k2^6*j3 + 4*i1^10*l2*j3 + 5*i1^8*i2*l2*j3 + 4*i1^6*i2^2*l2*j3 +
            3*i1^4*i2^3*l2*j3 + 5*i1^2*i2^4*l2*j3 + 3*i2^5*l2*j3 + i1^6*i2*j2*l2*j3 + i1^4*i2^2*j2*l2*j3 +
            4*i1^2*i2^3*j2*l2*j3 + 6*i2^4*j2*l2*j3 + 2*i2^3*j2^2*l2*j3 + 3*i1^4*j2^3*l2*j3 + 3*i1^2*i2*j2^3*l2*j3 +
            5*i2^2*j2^3*l2*j3 + 2*i1^2*j2^4*l2*j3 + i2*j2^4*l2*j3 + 4*j2^5*l2*j3 + 3*i1^8*k2*l2*j3 +
            6*i1^6*i2*k2*l2*j3 + 6*i1^4*i2^2*k2*l2*j3 + 2*i2^4*k2*l2*j3 + 4*i1^4*i2*j2*k2*l2*j3 +
            5*i1^2*i2^2*j2*k2*l2*j3 + 6*i2^3*j2*k2*l2*j3 + 2*i1^4*j2^2*k2*l2*j3 + i1^2*i2*j2^2*k2*l2*j3 +
            5*i2^2*j2^2*k2*l2*j3 + i1^2*j2^3*k2*l2*j3 + 6*i2*j2^3*k2*l2*j3 + 2*j2^4*k2*l2*j3 + 3*i1^6*k2^2*l2*j3 +
            2*i1^4*i2*k2^2*l2*j3 + 4*i2^3*k2^2*l2*j3 + 6*i1^4*j2*k2^2*l2*j3 + 6*i1^2*i2*j2*k2^2*l2*j3 +
            2*i2^2*j2*k2^2*l2*j3 + i1^2*j2^2*k2^2*l2*j3 + 3*i2*j2^2*k2^2*l2*j3 + 5*j2^3*k2^2*l2*j3 + 4*i1^4*k2^3*l2*j3
            + 3*i2^2*k2^3*l2*j3 + 4*i1^2*j2*k2^3*l2*j3 + 4*i2*j2*k2^3*l2*j3 + 3*j2^2*k2^3*l2*j3 + 3*i1^2*k2^4*l2*j3 +
            i2*k2^4*l2*j3 + 5*j2*k2^4*l2*j3 + 4*k2^5*l2*j3 + 4*i1^8*l2^2*j3 + 5*i1^6*i2*l2^2*j3 + 3*i1^4*i2^2*l2^2*j3
            + i1^2*i2^3*l2^2*j3 + 3*i2^4*l2^2*j3 + 4*i1^6*j2*l2^2*j3 + 3*i1^4*i2*j2*l2^2*j3 + 3*i1^2*i2^2*j2*l2^2*j3 +
            4*i2^3*j2*l2^2*j3 + 3*i1^4*j2^2*l2^2*j3 + 6*i1^2*i2*j2^2*l2^2*j3 + 5*i2^2*j2^2*l2^2*j3 +
            4*i1^2*j2^3*l2^2*j3 + i2*j2^3*l2^2*j3 + j2^4*l2^2*j3 + 3*i1^6*k2*l2^2*j3 + i2^3*k2*l2^2*j3 +
            6*i1^4*j2*k2*l2^2*j3 + 4*i1^2*i2*j2*k2*l2^2*j3 + i1^2*j2^2*k2*l2^2*j3 + 6*j2^3*k2*l2^2*j3 +
            2*i1^4*k2^2*l2^2*j3 + 3*i1^2*i2*k2^2*l2^2*j3 + 4*i2^2*k2^2*l2^2*j3 + 6*i2*j2*k2^2*l2^2*j3 +
            3*j2^2*k2^2*l2^2*j3 + i1^2*k2^3*l2^2*j3 + 6*i2*k2^3*l2^2*j3 + 5*k2^4*l2^2*j3 + 6*i1^6*l2^3*j3 +
            3*i1^4*i2*l2^3*j3 + 6*i2^3*l2^3*j3 + 3*i1^4*j2*l2^3*j3 + 6*i2^2*j2*l2^3*j3 + 2*i1^2*j2^2*l2^3*j3 +
            5*i2*j2^2*l2^3*j3 + 4*j2^3*l2^3*j3 + 6*i1^4*k2*l2^3*j3 + 6*i1^2*i2*k2*l2^3*j3 + 6*i2^2*k2*l2^3*j3 +
            2*i1^2*j2*k2*l2^3*j3 + 3*i2*j2*k2*l2^3*j3 + 3*j2^2*k2*l2^3*j3 + i1^2*k2^2*l2^3*j3 + 2*j2*k2^2*l2^3*j3 +
            3*k2^3*l2^3*j3 + i1^4*l2^4*j3 + 5*i1^2*i2*l2^4*j3 + 6*i2^2*l2^4*j3 + i1^2*j2*l2^4*j3 + 5*i2*j2*l2^4*j3 +
            4*j2^2*l2^4*j3 + 2*i1^2*k2*l2^4*j3 + j2*k2*l2^4*j3 + 2*k2^2*l2^4*j3 + 6*i1^2*l2^5*j3 + 5*i2*l2^5*j3 +
            6*k2*l2^5*j3 + 4*l2^6*j3 + 4*i1^9*i3*j3 + 5*i1^7*i2*i3*j3 + 2*i1^5*i2^2*i3*j3 + 5*i1^3*i2^3*i3*j3 +
            2*i1^7*j2*i3*j3 + 4*i1^5*i2*j2*i3*j3 + 6*i1^3*i2^2*j2*i3*j3 + i1^5*j2^2*i3*j3 + i1^3*i2*j2^2*i3*j3 +
            2*i1^3*j2^3*i3*j3 + i1^7*k2*i3*j3 + i1^3*i2^2*k2*i3*j3 + 5*i1^5*j2*k2*i3*j3 + 5*i1^3*i2*j2*k2*i3*j3 +
            i1^3*j2^2*k2*i3*j3 + 3*i1^3*i2*k2^2*i3*j3 + 2*i1*i2^2*k2^2*i3*j3 + 6*i1^3*j2*k2^2*i3*j3 +
            3*i1*i2*j2*k2^2*i3*j3 + 2*i1*j2^2*k2^2*i3*j3 + i1^3*k2^3*i3*j3 + 4*i1*i2*k2^3*i3*j3 + 4*i1*j2*k2^3*i3*j3 +
            3*i1*k2^4*i3*j3 + 4*i1^7*l2*i3*j3 + 4*i1^5*i2*l2*i3*j3 + 4*i1^3*i2^2*l2*i3*j3 + 3*i1*i2^3*l2*i3*j3 +
            i1^5*j2*l2*i3*j3 + 3*i1^3*i2*j2*l2*i3*j3 + 5*i1*i2^2*j2*l2*i3*j3 + 2*i1*i2*j2^2*l2*i3*j3 +
            4*i1*j2^3*l2*i3*j3 + 5*i1^5*k2*l2*i3*j3 + i1^3*i2*k2*l2*i3*j3 + i1*i2^2*k2*l2*i3*j3 +
            6*i1^3*j2*k2*l2*i3*j3 + 6*i1*i2*j2*k2*l2*i3*j3 + 6*i1^3*k2^2*l2*i3*j3 + 4*i1*j2*k2^2*l2*i3*j3 +
            5*i1^5*l2^2*i3*j3 + 6*i1^3*i2*l2^2*i3*j3 + 6*i1*i2*j2*l2^2*i3*j3 + i1*j2^2*l2^2*i3*j3 +
            4*i1^3*k2*l2^2*i3*j3 + 2*i1*i2*k2*l2^2*i3*j3 + 2*i1*j2*k2*l2^2*i3*j3 + 2*i1*k2^2*l2^2*i3*j3 +
            5*i1^3*l2^3*i3*j3 + 5*i1*i2*l2^3*i3*j3 + 6*i1*j2*l2^3*i3*j3 + 4*i1*k2*l2^3*i3*j3 + i1*l2^4*i3*j3 +
            6*i1^6*i3^2*j3 + 6*i1^4*i2*i3^2*j3 + i1^4*j2*i3^2*j3 + 2*i1^4*k2*i3^2*j3 + 5*i1^2*i2*k2*i3^2*j3 +
            2*i1^2*j2*k2*i3^2*j3 + 4*i1^2*k2^2*i3^2*j3 + 2*i2*k2^2*i3^2*j3 + 5*j2*k2^2*i3^2*j3 + 6*k2^3*i3^2*j3 +
            6*i1^4*l2*i3^2*j3 + 4*i1^2*i2*l2*i3^2*j3 + 3*i1^2*j2*l2*i3^2*j3 + 2*k2^2*l2*i3^2*j3 + 5*i1^2*l2^2*i3^2*j3
            + 3*i2*l2^2*i3^2*j3 + 4*j2*l2^2*i3^2*j3 + 5*k2*l2^2*i3^2*j3 + 6*i1^9*j3^2 + 5*i1^7*i2*j3^2 +
            5*i1^5*i2^2*j3^2 + 3*i1^3*i2^3*j3^2 + 3*i1*i2^4*j3^2 + 3*i1^7*j2*j3^2 + 2*i1^5*i2*j2*j3^2 +
            6*i1^3*i2^2*j2*j3^2 + 2*i1*i2^3*j2*j3^2 + 3*i1^5*j2^2*j3^2 + 4*i1*i2^2*j2^2*j3^2 + 5*i1^3*j2^3*j3^2 +
            2*i1*i2*j2^3*j3^2 + 3*i1*j2^4*j3^2 + 5*i1^7*k2*j3^2 + 4*i1^5*i2*k2*j3^2 + 4*i1^3*i2^2*k2*j3^2 +
            5*i1*i2^3*k2*j3^2 + 4*i1^5*j2*k2*j3^2 + 5*i1^3*i2*j2*k2*j3^2 + 5*i1*i2^2*j2*k2*j3^2 + 5*i1^3*j2^2*k2*j3^2
            + 3*i1*i2*j2^2*k2*j3^2 + i1*j2^3*k2*j3^2 + 6*i1^5*k2^2*j3^2 + 2*i1^3*i2*k2^2*j3^2 + i1*i2^2*k2^2*j3^2 +
            i1^3*j2*k2^2*j3^2 + 6*i1^3*k2^3*j3^2 + 6*i1*i2*k2^3*j3^2 + i1*j2*k2^3*j3^2 + 5*i1*k2^4*j3^2 +
            3*i1^7*l2*j3^2 + 4*i1^5*i2*l2*j3^2 + i1^3*i2^2*l2*j3^2 + 3*i1*i2^3*l2*j3^2 + 3*i1^5*j2*l2*j3^2 +
            6*i1*i2^2*j2*l2*j3^2 + 5*i1^3*j2^2*l2*j3^2 + 5*i1*j2^3*l2*j3^2 + 6*i1^5*k2*l2*j3^2 + 5*i1^3*i2*k2*l2*j3^2
            + 4*i1*i2^2*k2*l2*j3^2 + 4*i1^3*j2*k2*l2*j3^2 + 6*i1*i2*j2*k2*l2*j3^2 + 4*i1*j2^2*k2*l2*j3^2 +
            6*i1^3*k2^2*l2*j3^2 + 2*i1*i2*k2^2*l2*j3^2 + 4*i1*j2*k2^2*l2*j3^2 + 4*i1*k2^3*l2*j3^2 + 4*i1^5*l2^2*j3^2 +
            i1^3*i2*l2^2*j3^2 + 6*i1*i2^2*l2^2*j3^2 + 5*i1^3*j2*l2^2*j3^2 + 5*i1*i2*j2*l2^2*j3^2 + 4*i1*j2^2*l2^2*j3^2
            + 4*i1*i2*k2*l2^2*j3^2 + 3*i1*j2*k2*l2^2*j3^2 + 6*i1*k2^2*l2^2*j3^2 + 3*i1^3*l2^3*j3^2 + 2*i1*i2*l2^3*j3^2
            + 6*i1*j2*l2^3*j3^2 + 2*i1*k2*l2^3*j3^2 + 5*i1*l2^4*j3^2 + 6*i1^6*i3*j3^2 + 4*i1^4*i2*i3*j3^2 +
            5*i1^2*i2^2*i3*j3^2 + 5*i1^4*j2*i3*j3^2 + 4*i1^2*i2*j2*i3*j3^2 + 5*i1^2*j2^2*i3*j3^2 + 3*i1^4*k2*i3*j3^2 +
            i2^2*k2*i3*j3^2 + 3*i1^2*j2*k2*i3*j3^2 + 5*i2*j2*k2*i3*j3^2 + j2^2*k2*i3*j3^2 + i2*k2^2*i3*j3^2 +
            3*j2*k2^2*i3*j3^2 + 6*k2^3*i3*j3^2 + 6*i1^4*l2*i3*j3^2 + 4*i1^2*i2*l2*i3*j3^2 + 2*i2^2*l2*i3*j3^2 +
            2*i1^2*j2*l2*i3*j3^2 + 3*i2*j2*l2*i3*j3^2 + 2*j2^2*l2*i3*j3^2 + 5*i1^2*k2*l2*i3*j3^2 + 3*j2*k2*l2*i3*j3^2
            + 2*i1^2*l2^2*i3*j3^2 + 2*i2*l2^2*i3*j3^2 + 5*j2*l2^2*i3*j3^2 + 5*k2*l2^2*i3*j3^2 + 2*l2^3*i3*j3^2 +
            2*i1^3*i3^2*j3^2 + 2*i1*k2*i3^2*j3^2 + 4*i1*l2*i3^2*j3^2 + 5*i1^6*j3^3 + 4*i1^4*i2*j3^3 + 2*i1^2*i2^2*j3^3
            + i2^3*j3^3 + 2*i1^4*j2*j3^3 + 5*i1^2*i2*j2*j3^3 + 4*i2^2*j2*j3^3 + 3*i2*j2^2*j3^3 + 6*j2^3*j3^3 +
            4*i1^4*k2*j3^3 + 2*i2^2*k2*j3^3 + 3*i1^2*j2*k2*j3^3 + 2*i2*j2*k2*j3^3 + 3*j2^2*k2*j3^3 + i1^2*k2^2*j3^3 +
            3*i2*k2^2*j3^3 + 6*j2*k2^2*j3^3 + 4*k2^3*j3^3 + 4*i1^4*l2*j3^3 + 5*i1^2*i2*l2*j3^3 + 5*i2^2*l2*j3^3 +
            2*i1^2*j2*l2*j3^3 + 3*i2*j2*l2*j3^3 + 6*j2^2*l2*j3^3 + 4*i1^2*k2*l2*j3^3 + 4*i2*k2*l2*j3^3 +
            5*j2*k2*l2*j3^3 + 6*k2^2*l2*j3^3 + 3*i2*l2^2*j3^3 + 3*j2*l2^2*j3^3 + 5*k2*l2^2*j3^3 + 6*l2^3*j3^3 +
            6*i1^3*i3*j3^3 + 3*i1*i2*i3*j3^3 + 4*i1*j2*i3*j3^3 + i1*l2*i3*j3^3 + 5*i1^3*j3^4 + 3*i1*i2*j3^4 +
            4*i1*j2*j3^4 + i1*k2*j3^4 + 3*i1*l2*j3^4 + 6*i1^12*k3 + 6*i1^10*i2*k3 + 3*i1^8*i2^2*k3 + i1^6*i2^3*k3 +
            5*i1^4*i2^4*k3 + 4*i1^2*i2^5*k3 + 5*i1^8*i2*j2*k3 + 5*i1^6*i2^2*j2*k3 + 2*i1^4*i2^3*j2*k3 +
            i1^2*i2^4*j2*k3 + 5*i1^8*j2^2*k3 + 3*i1^6*i2*j2^2*k3 + 6*i1^4*i2^2*j2^2*k3 + 5*i1^2*i2^3*j2^2*k3 +
            5*i1^6*j2^3*k3 + 4*i1^4*i2*j2^3*k3 + 2*i1^2*i2^2*j2^3*k3 + 4*i1^4*j2^4*k3 + 6*i1^2*i2*j2^4*k3 +
            3*i1^2*j2^5*k3 + i1^10*k2*k3 + 4*i1^8*i2*k2*k3 + 3*i1^6*i2^2*k2*k3 + 6*i1^4*i2^3*k2*k3 + 2*i1^2*i2^4*k2*k3
            + i2^5*k2*k3 + i1^8*j2*k2*k3 + i1^6*i2*j2*k2*k3 + 3*i1^4*i2^2*j2*k2*k3 + 3*i1^2*i2^3*j2*k2*k3 +
            2*i2^4*j2*k2*k3 + 3*i1^6*j2^2*k2*k3 + 5*i1^4*i2*j2^2*k2*k3 + 3*i2^3*j2^2*k2*k3 + 4*i1^2*i2*j2^3*k2*k3 +
            4*i2^2*j2^3*k2*k3 + 5*i1^2*j2^4*k2*k3 + 5*i2*j2^4*k2*k3 + 6*j2^5*k2*k3 + 6*i1^6*i2*k2^2*k3 +
            2*i1^4*i2^2*k2^2*k3 + i1^2*i2^3*k2^2*k3 + 3*i2^4*k2^2*k3 + 2*i1^6*j2*k2^2*k3 + 6*i1^4*i2*j2*k2^2*k3 +
            5*i2^3*j2*k2^2*k3 + 6*i1^4*j2^2*k2^2*k3 + i1^2*i2*j2^2*k2^2*k3 + 2*i2^2*j2^2*k2^2*k3 + 5*i1^2*j2^3*k2^2*k3
            + 4*i2*j2^3*k2^2*k3 + 3*i1^6*k2^3*k3 + 3*i1^4*i2*k2^3*k3 + i1^2*i2^2*k2^3*k3 + 6*i2^3*k2^3*k3 +
            i1^4*j2*k2^3*k3 + 3*i1^2*i2*j2*k2^3*k3 + 4*i2^2*j2*k2^3*k3 + i1^2*j2^2*k2^3*k3 + 3*i2*j2^2*k2^3*k3 +
            j2^3*k2^3*k3 + 2*i1^2*i2*k2^4*k3 + 2*i2^2*k2^4*k3 + 3*i1^2*j2*k2^4*k3 + 4*i2*j2*k2^4*k3 + 6*j2^2*k2^4*k3 +
            6*i1^2*k2^5*k3 + 3*i2*k2^5*k3 + 5*j2*k2^5*k3 + 2*k2^6*k3 + 3*i1^10*l2*k3 + 4*i1^8*i2*l2*k3 +
            4*i1^6*i2^2*l2*k3 + 6*i1^2*i2^4*l2*k3 + 4*i2^5*l2*k3 + i1^8*j2*l2*k3 + i1^4*i2^2*j2*l2*k3 +
            2*i1^2*i2^3*j2*l2*k3 + i2^4*j2*l2*k3 + 5*i1^6*j2^2*l2*k3 + 4*i1^4*i2*j2^2*l2*k3 + 5*i2^3*j2^2*l2*k3 +
            2*i1^4*j2^3*l2*k3 + 5*i1^2*i2*j2^3*l2*k3 + 2*i2^2*j2^3*l2*k3 + i1^2*j2^4*l2*k3 + 6*i2*j2^4*l2*k3 +
            3*j2^5*l2*k3 + i1^8*k2*l2*k3 + 6*i1^6*i2*k2*l2*k3 + 6*i1^4*i2^2*k2*l2*k3 + 3*i1^6*j2*k2*l2*k3 +
            5*i1^4*i2*j2*k2*l2*k3 + 2*i1^2*i2^2*j2*k2*l2*k3 + 6*i2^3*j2*k2*l2*k3 + i1^4*j2^2*k2*l2*k3 +
            i1^2*i2*j2^2*k2*l2*k3 + 3*i2^2*j2^2*k2*l2*k3 + 4*i1^2*j2^3*k2*l2*k3 + 4*i2*j2^3*k2*l2*k3 + j2^4*k2*l2*k3 +
            2*i1^6*k2^2*l2*k3 + 2*i1^4*i2*k2^2*l2*k3 + 3*i1^2*i2^2*k2^2*l2*k3 + 2*i2^3*k2^2*l2*k3 +
            2*i1^4*j2*k2^2*l2*k3 + 6*i2^2*j2*k2^2*l2*k3 + 4*i1^2*j2^2*k2^2*l2*k3 + 4*i2*j2^2*k2^2*l2*k3 +
            2*j2^3*k2^2*l2*k3 + 4*i1^2*i2*k2^3*l2*k3 + 6*i2^2*k2^3*l2*k3 + i2*j2*k2^3*l2*k3 + 2*j2^2*k2^3*l2*k3 +
            3*i1^2*k2^4*l2*k3 + 6*j2*k2^4*l2*k3 + 2*k2^5*l2*k3 + i1^8*l2^2*k3 + i1^6*i2*l2^2*k3 + 5*i1^4*i2^2*l2^2*k3
            + 4*i1^2*i2^3*l2^2*k3 + 2*i2^4*l2^2*k3 + 6*i1^2*i2^2*j2*l2^2*k3 + 6*i2^3*j2*l2^2*k3 + 3*i1^4*j2^2*l2^2*k3
            + 4*i1^2*i2*j2^2*l2^2*k3 + 5*i2^2*j2^2*l2^2*k3 + 6*i2*j2^3*l2^2*k3 + 2*j2^4*l2^2*k3 + 5*i1^6*k2*l2^2*k3 +
            2*i1^2*i2^2*k2*l2^2*k3 + 2*i2^3*k2*l2^2*k3 + i2^2*j2*k2*l2^2*k3 + 4*i1^2*j2^2*k2*l2^2*k3 +
            3*i2*j2^2*k2*l2^2*k3 + j2^3*k2*l2^2*k3 + i1^4*k2^2*l2^2*k3 + 3*i1^2*i2*k2^2*l2^2*k3 + 4*i2^2*k2^2*l2^2*k3
            + i1^2*j2*k2^2*l2^2*k3 + i2*j2*k2^2*l2^2*k3 + 2*j2^2*k2^2*l2^2*k3 + 2*i1^2*k2^3*l2^2*k3 +
            6*i2*k2^3*l2^2*k3 + 5*j2*k2^3*l2^2*k3 + i1^6*l2^3*k3 + i1^4*i2*l2^3*k3 + 6*i1^2*i2^2*l2^3*k3 +
            3*i2^3*l2^3*k3 + 3*i1^4*j2*l2^3*k3 + 6*i1^2*i2*j2*l2^3*k3 + 5*i2^2*j2*l2^3*k3 + 2*i1^2*j2^2*l2^3*k3 +
            2*i2*j2^2*l2^3*k3 + 4*j2^3*l2^3*k3 + 5*i1^4*k2*l2^3*k3 + 3*i1^2*i2*k2*l2^3*k3 + i2^2*k2*l2^3*k3 +
            2*i2*j2*k2*l2^3*k3 + 2*j2^2*k2*l2^3*k3 + 2*i1^2*k2^2*l2^3*k3 + 5*i2*k2^2*l2^3*k3 + 3*j2*k2^2*l2^3*k3 +
            3*k2^3*l2^3*k3 + 3*i1^4*l2^4*k3 + i1^2*i2*l2^4*k3 + 5*i2^2*l2^4*k3 + 5*i1^2*j2*l2^4*k3 + 6*i2*j2*l2^4*k3 +
            j2^2*l2^4*k3 + 6*i1^2*k2*l2^4*k3 + 2*i2*k2*l2^4*k3 + 4*j2*k2*l2^4*k3 + 5*j2*l2^5*k3 + 3*k2*l2^5*k3 +
            l2^6*k3 + 3*i1^9*j3*k3 + 3*i1^7*i2*j3*k3 + 6*i1^5*i2^2*j3*k3 + 2*i1*i2^4*j3*k3 + 2*i1^7*j2*j3*k3 +
            i1^5*i2*j2*j3*k3 + i1^3*i2^2*j2*j3*k3 + 6*i1*i2^3*j2*j3*k3 + 6*i1^5*j2^2*j3*k3 + 5*i1^3*i2*j2^2*j3*k3 +
            5*i1*i2^2*j2^2*j3*k3 + i1^3*j2^3*j3*k3 + 6*i1*i2*j2^3*j3*k3 + 2*i1*j2^4*j3*k3 + i1^7*k2*j3*k3 +
            4*i1^3*i2^2*k2*j3*k3 + 6*i1^5*j2*k2*j3*k3 + 4*i1*i2^2*j2*k2*j3*k3 + 6*i1*i2*j2^2*k2*j3*k3 +
            4*i1*j2^3*k2*j3*k3 + 6*i1^5*k2^2*j3*k3 + i1^3*i2*k2^2*j3*k3 + 6*i1*i2^2*k2^2*j3*k3 + 5*i1^3*j2*k2^2*j3*k3
            + 5*i1*i2*j2*k2^2*j3*k3 + i1*j2^2*k2^2*j3*k3 + i1^3*k2^3*j3*k3 + 5*i1*i2*k2^3*j3*k3 + 2*i1*j2*k2^3*j3*k3 +
            4*i1*k2^4*j3*k3 + 2*i1^7*l2*j3*k3 + 4*i1^5*i2*l2*j3*k3 + i1*i2^3*l2*j3*k3 + 4*i1^5*j2*l2*j3*k3 +
            3*i1^3*i2*j2*l2*j3*k3 + 5*i1*i2^2*j2*l2*j3*k3 + 2*i1^3*j2^2*l2*j3*k3 + i1*i2*j2^2*l2*j3*k3 +
            2*i1^5*k2*l2*j3*k3 + 5*i1^3*i2*k2*l2*j3*k3 + 4*i1*i2^2*k2*l2*j3*k3 + 6*i1^3*j2*k2*l2*j3*k3 +
            5*i1*j2^2*k2*l2*j3*k3 + 4*i1^3*k2^2*l2*j3*k3 + i1*i2*k2^2*l2*j3*k3 + 3*i1*j2*k2^2*l2*j3*k3 +
            5*i1*k2^3*l2*j3*k3 + 5*i1^5*l2^2*j3*k3 + 3*i1^3*i2*l2^2*j3*k3 + 3*i1*i2^2*l2^2*j3*k3 +
            5*i1^3*j2*l2^2*j3*k3 + 6*i1*i2*j2*l2^2*j3*k3 + 2*i1*j2^2*l2^2*j3*k3 + 2*i1*i2*k2*l2^2*j3*k3 +
            i1*j2*k2*l2^2*j3*k3 + 4*i1*k2^2*l2^2*j3*k3 + 6*i1^3*l2^3*j3*k3 + 4*i1*i2*l2^3*j3*k3 + 3*i1*j2*l2^3*j3*k3 +
            5*i1*l2^4*j3*k3 + 4*i1^6*j3^2*k3 + 6*i1^2*i2^2*j3^2*k3 + 3*i2^3*j3^2*k3 + 3*i1^4*j2*j3^2*k3 +
            2*i1^2*i2*j2*j3^2*k3 + 5*i2^2*j2*j3^2*k3 + 6*i1^2*j2^2*j3^2*k3 + 2*i2*j2^2*j3^2*k3 + 4*j2^3*j3^2*k3 +
            3*i1^4*k2*j3^2*k3 + 4*i1^2*i2*k2*j3^2*k3 + 4*i2^2*k2*j3^2*k3 + 4*i1^2*j2*k2*j3^2*k3 + 5*i2*j2*k2*j3^2*k3 +
            5*j2^2*k2*j3^2*k3 + 5*i1^2*k2^2*j3^2*k3 + 6*j2*k2^2*j3^2*k3 + 3*i1^4*l2*j3^2*k3 + 2*i1^2*i2*l2*j3^2*k3 +
            6*i2*j2*l2*j3^2*k3 + j2^2*l2*j3^2*k3 + 6*i1^2*k2*l2*j3^2*k3 + 2*k2^2*l2*j3^2*k3 + 3*i1^2*l2^2*j3^2*k3 +
            j2*l2^2*j3^2*k3 + k2*l2^2*j3^2*k3 + 4*l2^3*j3^2*k3 + 3*i1^3*j3^3*k3 + 2*i1*i2*j3^3*k3 + i1*j2*j3^3*k3 +
            5*i1*k2*j3^3*k3 + 3*i1*l2*j3^3*k3 + i1^9*k3^2 + 2*i1^7*i2*k3^2 + 2*i1^5*i2^2*k3^2 + 5*i1^3*i2^3*k3^2 +
            3*i1*i2^4*k3^2 + i1^7*j2*k3^2 + i1^5*i2*j2*k3^2 + 2*i1^3*i2^2*j2*k3^2 + 2*i1*i2^3*j2*k3^2 +
            3*i1^5*j2^2*k3^2 + 2*i1^3*i2*j2^2*k3^2 + 4*i1*i2^2*j2^2*k3^2 + 5*i1^3*j2^3*k3^2 + 2*i1*i2*j2^3*k3^2 +
            3*i1*j2^4*k3^2 + 4*i1^7*k2*k3^2 + 4*i1^5*i2*k2*k3^2 + 3*i1^3*i2^2*k2*k3^2 + 5*i1*i2^3*k2*k3^2 +
            6*i1^5*j2*k2*k3^2 + 2*i1^3*i2*j2*k2*k3^2 + 5*i1*i2^2*j2*k2*k3^2 + 3*i1*i2*j2^2*k2*k3^2 + i1*j2^3*k2*k3^2 +
            i1^5*k2^2*k3^2 + 5*i1^3*i2*k2^2*k3^2 + 3*i1*i2^2*k2^2*k3^2 + 6*i1^3*j2*k2^2*k3^2 + 4*i1*i2*j2*k2^2*k3^2 +
            3*i1*j2^2*k2^2*k3^2 + 5*i1^3*k2^3*k3^2 + 2*i1*i2*k2^3*k3^2 + 6*i1*j2*k2^3*k3^2 + 2*i1*k2^4*k3^2 +
            5*i1^7*l2*k3^2 + i1^5*i2*l2*k3^2 + i1^3*i2^2*l2*k3^2 + i1*i2^3*l2*k3^2 + 4*i1^5*j2*l2*k3^2 +
            5*i1*i2^2*j2*l2*k3^2 + 6*i1^3*j2^2*l2*k3^2 + i1*i2*j2^2*l2*k3^2 + 4*i1^5*k2*l2*k3^2 + 2*i1^3*i2*k2*l2*k3^2
            + 6*i1*i2^2*k2*l2*k3^2 + 3*i1^3*j2*k2*l2*k3^2 + 2*i1*i2*j2*k2*l2*k3^2 + 6*i1^3*k2^2*l2*k3^2 +
            6*i1*i2*k2^2*l2*k3^2 + 4*i1^5*l2^2*k3^2 + 6*i1^3*i2*l2^2*k3^2 + 3*i1*i2^2*l2^2*k3^2 + 3*i1^3*j2*l2^2*k3^2
            + 2*i1*i2*j2*l2^2*k3^2 + 5*i1*j2^2*l2^2*k3^2 + 3*i1^3*k2*l2^2*k3^2 + 2*i1*j2*k2*l2^2*k3^2 +
            2*i1*k2^2*l2^2*k3^2 + i1*i2*l2^3*k3^2 + 2*i1*j2*l2^3*k3^2 + i1*k2*l2^3*k3^2 + i1*l2^4*k3^2 +
            5*i1^6*j3*k3^2 + 6*i1^4*i2*j3*k3^2 + 6*i1^2*i2^2*j3*k3^2 + 4*i2^3*j3*k3^2 + i1^2*i2*j2*j3*k3^2 +
            2*i2^2*j2*j3*k3^2 + 5*i2*j2^2*j3*k3^2 + 3*j2^3*j3*k3^2 + 2*i1^4*k2*j3*k3^2 + i2^2*k2*j3*k3^2 +
            2*i1^2*j2*k2*j3*k3^2 + 3*i2*j2*k2*j3*k3^2 + 3*j2^2*k2*j3*k3^2 + 3*i1^2*k2^2*j3*k3^2 + 4*i2*k2^2*j3*k3^2 +
            5*j2*k2^2*j3*k3^2 + 4*k2^3*j3*k3^2 + 3*i1^2*i2*l2*j3*k3^2 + 4*i1^2*j2*l2*j3*k3^2 + 4*i2*j2*l2*j3*k3^2 +
            3*j2^2*l2*j3*k3^2 + i1^2*k2*l2*j3*k3^2 + 3*i2*k2*l2*j3*k3^2 + 2*j2*k2*l2*j3*k3^2 + 4*k2^2*l2*j3*k3^2 +
            6*i1^2*l2^2*j3*k3^2 + 2*i2*l2^2*j3*k3^2 + 6*j2*l2^2*j3*k3^2 + 5*k2*l2^2*j3*k3^2 + 2*l2^3*j3*k3^2 +
            6*i1^3*j3^2*k3^2 + 5*i1*i2*j3^2*k3^2 + 2*i1*j2*j3^2*k3^2 + 3*i1*k2*j3^2*k3^2 + i1*l2*j3^2*k3^2 +
            4*j3^3*k3^2 + 4*i1^6*k3^3 + 6*i1^4*i2*k3^3 + i1^2*i2^2*k3^3 + i1^4*j2*k3^3 + 4*i1^2*i2*j2*k3^3 +
            2*i1^2*j2^2*k3^3 + 5*i1^4*k2*k3^3 + i1^2*i2*k2*k3^3 + 3*i2^2*k2*k3^3 + 2*i1^2*j2*k2*k3^3 + 2*i2*j2*k2*k3^3
            + 2*j2^2*k2*k3^3 + 4*i1^2*k2^2*k3^3 + i2*k2^2*k3^3 + 2*j2*k2^2*k3^3 + 2*k2^3*k3^3 + 6*i1^4*l2*k3^3 +
            4*i1^2*i2*l2*k3^3 + 2*i2^2*l2*k3^3 + 4*i2*j2*l2*k3^3 + j2^2*l2*k3^3 + 3*i1^2*k2*l2*k3^3 + 2*i2*k2*l2*k3^3
            + 3*i1^2*l2^2*k3^3 + 5*i2*l2^2*k3^3 + j2*l2^2*k3^3 + 2*k2*l2^2*k3^3 + l2^3*k3^3 + 2*i1^3*j3*k3^3 +
            2*i1*j2*j3*k3^3 + 5*i1*k2*j3*k3^3 + 2*i1*l2*j3*k3^3 + 4*i1*i2*k3^4 + 4*i1*j2*k3^4 + 6*i1*l2*k3^4 + j3*k3^4
            + 6*k3^5 + 6*i1^12*l3 + 5*i1^10*i2*l3 + 5*i1^6*i2^3*l3 + 5*i1^4*i2^4*l3 + 3*i1^2*i2^5*l3 + 6*i2^6*l3 +
            4*i1^10*j2*l3 + 5*i1^8*i2*j2*l3 + 4*i1^6*i2^2*j2*l3 + 2*i1^4*i2^3*j2*l3 + 2*i1^2*i2^4*j2*l3 + 6*i2^5*j2*l3
            + i1^8*j2^2*l3 + 6*i1^6*i2*j2^2*l3 + 5*i1^4*i2^2*j2^2*l3 + 4*i1^2*i2^3*j2^2*l3 + 6*i2^4*j2^2*l3 +
            i1^6*j2^3*l3 + 6*i1^4*i2*j2^3*l3 + 2*i1^2*i2^2*j2^3*l3 + 6*i2^3*j2^3*l3 + 3*i1^4*j2^4*l3 +
            3*i1^2*i2*j2^4*l3 + 6*i2^2*j2^4*l3 + 6*i2*j2^5*l3 + 6*j2^6*l3 + 3*i1^10*k2*l3 + i1^6*i2^2*k2*l3 +
            3*i1^4*i2^3*k2*l3 + 4*i2^5*k2*l3 + 2*i1^8*j2*k2*l3 + 2*i1^6*i2*j2*k2*l3 + 2*i1^4*i2^2*j2*k2*l3 +
            5*i1^2*i2^3*j2*k2*l3 + 5*i2^4*j2*k2*l3 + 2*i1^6*j2^2*k2*l3 + 5*i1^4*i2*j2^2*k2*l3 + 4*i1^2*i2^2*j2^2*k2*l3
            + 3*i2^3*j2^2*k2*l3 + i1^4*j2^3*k2*l3 + 5*i1^2*i2*j2^3*k2*l3 + 5*i2^2*j2^3*k2*l3 + 4*i2*j2^4*k2*l3 +
            6*i1^6*i2*k2^2*l3 + 4*i1^4*i2^2*k2^2*l3 + 6*i1^2*i2^3*k2^2*l3 + i2^4*k2^2*l3 + 4*i1^6*j2*k2^2*l3 +
            6*i1^2*i2^2*j2*k2^2*l3 + 6*i2^3*j2*k2^2*l3 + 6*i1^4*j2^2*k2^2*l3 + 3*i1^2*i2*j2^2*k2^2*l3 +
            4*i2^2*j2^2*k2^2*l3 + i1^2*j2^3*k2^2*l3 + 5*i2*j2^3*k2^2*l3 + 5*j2^4*k2^2*l3 + 4*i1^4*i2*k2^3*l3 +
            i2^3*k2^3*l3 + 5*i1^4*j2*k2^3*l3 + 5*i1^2*i2*j2*k2^3*l3 + 5*i2^2*j2*k2^3*l3 + 5*i1^2*j2^2*k2^3*l3 +
            i2*j2^2*k2^3*l3 + 3*j2^3*k2^3*l3 + 6*i1^4*k2^4*l3 + 4*i1^2*i2*k2^4*l3 + i2^2*k2^4*l3 + 2*i2*j2*k2^4*l3 +
            2*i1^2*k2^5*l3 + 4*i2*k2^5*l3 + 6*j2*k2^5*l3 + k2^6*l3 + 4*i1^8*i2*l2*l3 + 4*i1^6*i2^2*l2*l3 +
            5*i1^4*i2^3*l2*l3 + 5*i1^2*i2^4*l2*l3 + i2^5*l2*l3 + 2*i1^8*j2*l2*l3 + 4*i1^6*i2*j2*l2*l3 +
            i1^4*i2^2*j2*l2*l3 + 4*i1^2*i2^3*j2*l2*l3 + i2^4*j2*l2*l3 + i1^6*j2^2*l2*l3 + 5*i1^4*j2^3*l2*l3 +
            3*i1^2*i2*j2^3*l2*l3 + 5*i2^2*j2^3*l2*l3 + 2*i1^2*j2^4*l2*l3 + 2*i2*j2^4*l2*l3 + 5*j2^5*l2*l3 +
            6*i1^8*k2*l2*l3 + 4*i1^6*i2*k2*l2*l3 + 5*i1^2*i2^3*k2*l2*l3 + 3*i1^6*j2*k2*l2*l3 + i1^4*i2*j2*k2*l2*l3 +
            6*i1^2*i2^2*j2*k2*l2*l3 + 5*i2^3*j2*k2*l2*l3 + 2*i1^4*j2^2*k2*l2*l3 + 5*i1^2*i2*j2^2*k2*l2*l3 +
            4*i2^2*j2^2*k2*l2*l3 + 5*i2*j2^3*k2*l2*l3 + 5*i1^6*k2^2*l2*l3 + 2*i1^4*i2*k2^2*l2*l3 +
            6*i1^2*i2^2*k2^2*l2*l3 + 3*i2^3*k2^2*l2*l3 + 6*i1^4*j2*k2^2*l2*l3 + 3*i1^2*i2*j2*k2^2*l2*l3 +
            4*i2^2*j2*k2^2*l2*l3 + 5*i2*j2^2*k2^2*l2*l3 + 4*j2^3*k2^2*l2*l3 + 3*i1^4*k2^3*l2*l3 + 5*i1^2*i2*k2^3*l2*l3
            + 4*i2^2*k2^3*l2*l3 + 6*i1^2*j2*k2^3*l2*l3 + i2*j2*k2^3*l2*l3 + j2^2*k2^3*l2*l3 + i1^2*k2^4*l2*l3 +
            6*i2*k2^4*l2*l3 + 3*j2*k2^4*l2*l3 + 3*k2^5*l2*l3 + 6*i1^8*l2^2*l3 + 4*i1^6*i2*l2^2*l3 +
            5*i1^4*i2^2*l2^2*l3 + 2*i1^2*i2^3*l2^2*l3 + 6*i1^6*j2*l2^2*l3 + 2*i1^2*i2^2*j2*l2^2*l3 + i2^3*j2*l2^2*l3 +
            6*i1^4*j2^2*l2^2*l3 + 5*i1^2*i2*j2^2*l2^2*l3 + i2^2*j2^2*l2^2*l3 + 3*i1^2*j2^3*l2^2*l3 + 2*i2*j2^3*l2^2*l3
            + 3*j2^4*l2^2*l3 + 2*i1^6*k2*l2^2*l3 + 4*i1^4*i2*k2*l2^2*l3 + 2*i1^2*i2^2*k2*l2^2*l3 + 2*i2^3*k2*l2^2*l3 +
            i1^4*j2*k2*l2^2*l3 + 5*i2^2*j2*k2*l2^2*l3 + 3*i1^2*j2^2*k2*l2^2*l3 + 4*i2*j2^2*k2*l2^2*l3 +
            2*j2^3*k2*l2^2*l3 + 6*i1^4*k2^2*l2^2*l3 + 3*i1^2*i2*k2^2*l2^2*l3 + 5*i2^2*k2^2*l2^2*l3 +
            3*i1^2*j2*k2^2*l2^2*l3 + 4*i2*j2*k2^2*l2^2*l3 + 5*j2^2*k2^2*l2^2*l3 + 3*i1^2*k2^3*l2^2*l3 +
            i2*k2^3*l2^2*l3 + 6*j2*k2^3*l2^2*l3 + 4*k2^4*l2^2*l3 + 3*i1^6*l2^3*l3 + 4*i1^4*i2*l2^3*l3 +
            5*i1^2*i2^2*l2^3*l3 + 2*i2^3*l2^3*l3 + i1^4*j2*l2^3*l3 + 4*i1^2*i2*j2*l2^3*l3 + 5*i2^2*j2*l2^3*l3 +
            2*i1^2*j2^2*l2^3*l3 + 4*i2*j2^2*l2^3*l3 + 4*j2^3*l2^3*l3 + 3*i1^4*k2*l2^3*l3 + 6*i1^2*i2*k2*l2^3*l3 +
            5*i2^2*k2*l2^3*l3 + 2*i2*j2*k2*l2^3*l3 + 4*j2^2*k2*l2^3*l3 + i1^2*k2^2*l2^3*l3 + 5*j2*k2^2*l2^3*l3 +
            k2^3*l2^3*l3 + 6*i1^4*l2^4*l3 + 5*i2^2*l2^4*l3 + 5*i1^2*j2*l2^4*l3 + 2*i2*j2*l2^4*l3 + 4*j2^2*l2^4*l3 +
            6*i1^2*k2*l2^4*l3 + 3*i2*k2*l2^4*l3 + j2*k2*l2^4*l3 + 4*k2^2*l2^4*l3 + 6*i1^2*l2^5*l3 + 3*i2*l2^5*l3 +
            2*j2*l2^5*l3 + 3*k2*l2^5*l3 + 6*l2^6*l3 + 3*i1^7*i2*i3*l3 + i1^5*i2^2*i3*l3 + 3*i1^3*i2^3*i3*l3 +
            3*i1*i2^4*i3*l3 + 3*i1^7*j2*i3*l3 + 5*i1^5*i2*j2*i3*l3 + 3*i1^3*i2^2*j2*i3*l3 + 2*i1*i2^3*j2*i3*l3 +
            3*i1^5*j2^2*i3*l3 + 6*i1^3*i2*j2^2*i3*l3 + 4*i1*i2^2*j2^2*i3*l3 + 2*i1^3*j2^3*i3*l3 + 2*i1*i2*j2^3*i3*l3 +
            3*i1*j2^4*i3*l3 + 2*i1^7*k2*i3*l3 + i1^5*i2*k2*i3*l3 + 5*i1^3*i2^2*k2*i3*l3 + 6*i1*i2^3*k2*i3*l3 +
            2*i1^5*j2*k2*i3*l3 + 3*i1^3*i2*j2*k2*i3*l3 + 6*i1*i2^2*j2*k2*i3*l3 + 3*i1^3*j2^2*k2*i3*l3 +
            5*i1*i2*j2^2*k2*i3*l3 + 4*i1*j2^3*k2*i3*l3 + 4*i1^5*k2^2*i3*l3 + 3*i1^3*i2*k2^2*i3*l3 + i1*i2^2*k2^2*i3*l3
            + 6*i1^3*j2*k2^2*i3*l3 + 4*i1*i2*j2*k2^2*i3*l3 + 4*i1*j2^2*k2^2*i3*l3 + i1^3*k2^3*i3*l3 +
            2*i1*i2*k2^3*i3*l3 + 6*i1*j2*k2^3*i3*l3 + 6*i1*k2^4*i3*l3 + 2*i1^7*l2*i3*l3 + 3*i1^5*i2*l2*i3*l3 +
            6*i1^3*i2^2*l2*i3*l3 + 6*i1*i2^3*l2*i3*l3 + i1^5*j2*l2*i3*l3 + i1^3*i2*j2*l2*i3*l3 + 2*i1*i2^2*j2*l2*i3*l3
            + 4*i1^3*j2^2*l2*i3*l3 + 6*i1*i2*j2^2*l2*i3*l3 + 3*i1^5*k2*l2*i3*l3 + i1^3*i2*k2*l2*i3*l3 +
            6*i1*i2^2*k2*l2*i3*l3 + 3*i1^3*j2*k2*l2*i3*l3 + 5*i1*i2*j2*k2*l2*i3*l3 + 6*i1*j2^2*k2*l2*i3*l3 +
            5*i1^3*k2^2*l2*i3*l3 + 3*i1*i2*k2^2*l2*i3*l3 + 4*i1*j2*k2^2*l2*i3*l3 + i1*k2^3*l2*i3*l3 +
            6*i1^3*i2*l2^2*i3*l3 + 6*i1*i2^2*l2^2*i3*l3 + 5*i1^3*j2*l2^2*i3*l3 + i1*i2*j2*l2^2*i3*l3 +
            3*i1*j2^2*l2^2*i3*l3 + 6*i1^3*k2*l2^2*i3*l3 + i1*i2*k2*l2^2*i3*l3 + i1*k2^2*l2^2*i3*l3 + i1^3*l2^3*i3*l3 +
            3*i1*i2*l2^3*i3*l3 + 6*i1*j2*l2^3*i3*l3 + 2*i1*k2*l2^3*i3*l3 + 4*i1*l2^4*i3*l3 + i1^6*i3^2*l3 +
            5*i1^4*i2*i3^2*l3 + 3*i1^2*i2^2*i3^2*l3 + 5*i1^4*j2*i3^2*l3 + i1^2*i2*j2*i3^2*l3 + 3*i1^2*j2^2*i3^2*l3 +
            4*i1^4*k2*i3^2*l3 + 5*i1^2*i2*k2*i3^2*l3 + 5*i2^2*k2*i3^2*l3 + 4*i2*j2*k2*i3^2*l3 + 5*j2^2*k2*i3^2*l3 +
            3*i1^2*k2^2*i3^2*l3 + 2*i2*k2^2*i3^2*l3 + j2*k2^2*i3^2*l3 + 2*k2^3*i3^2*l3 + 3*i1^4*l2*i3^2*l3 +
            3*i1^2*i2*l2*i3^2*l3 + 3*i2^2*l2*i3^2*l3 + 5*i1^2*j2*l2*i3^2*l3 + i2*j2*l2*i3^2*l3 + 3*j2^2*l2*i3^2*l3 +
            2*i1^2*k2*l2*i3^2*l3 + 6*i2*k2*l2*i3^2*l3 + 5*j2*k2*l2*i3^2*l3 + 4*k2^2*l2*i3^2*l3 + 2*i2*l2^2*i3^2*l3 +
            4*j2*l2^2*i3^2*l3 + 5*k2*l2^2*i3^2*l3 + 4*l2^3*i3^2*l3 + 3*i1^3*i3^3*l3 + 4*i1*k2*i3^3*l3 +
            3*i1*l2*i3^3*l3 + 2*i1^9*j3*l3 + 3*i1^7*i2*j3*l3 + 4*i1^5*i2^2*j3*l3 + 5*i1^3*i2^3*j3*l3 + 2*i1*i2^4*j3*l3
            + 4*i1^7*j2*j3*l3 + 3*i1^5*i2*j2*j3*l3 + 2*i1^3*i2^2*j2*j3*l3 + 2*i1*i2^3*j2*j3*l3 + 2*i1^5*j2^2*j3*l3 +
            2*i1^3*i2*j2^2*j3*l3 + 3*i1*i2^2*j2^2*j3*l3 + 5*i1^3*j2^3*j3*l3 + i1*i2*j2^3*j3*l3 + 6*i1*j2^4*j3*l3 +
            i1^7*k2*j3*l3 + 3*i1^5*i2*k2*j3*l3 + 3*i1^3*i2^2*k2*j3*l3 + 2*i1^5*j2*k2*j3*l3 + 3*i1^3*i2*j2*k2*j3*l3 +
            6*i1*i2^2*j2*k2*j3*l3 + i1^3*j2^2*k2*j3*l3 + i1*i2*j2^2*k2*j3*l3 + 2*i1^5*k2^2*j3*l3 +
            2*i1^3*i2*k2^2*j3*l3 + 5*i1*i2*j2*k2^2*j3*l3 + 6*i1*j2^2*k2^2*j3*l3 + 5*i1^3*k2^3*j3*l3 +
            3*i1*i2*k2^3*j3*l3 + i1*j2*k2^3*j3*l3 + 3*i1*k2^4*j3*l3 + 4*i1^5*i2*l2*j3*l3 + i1^3*i2^2*l2*j3*l3 +
            3*i1*i2^3*l2*j3*l3 + 6*i1^5*j2*l2*j3*l3 + 2*i1^3*i2*j2*l2*j3*l3 + 3*i1*i2^2*j2*l2*j3*l3 + i1*j2^3*l2*j3*l3
            + 3*i1^5*k2*l2*j3*l3 + i1^3*i2*k2*l2*j3*l3 + 3*i1^3*j2*k2*l2*j3*l3 + 6*i1*i2*j2*k2*l2*j3*l3 +
            3*i1^3*k2^2*l2*j3*l3 + 5*i1*i2*k2^2*l2*j3*l3 + 5*i1*j2*k2^2*l2*j3*l3 + 3*i1*k2^3*l2*j3*l3 +
            2*i1^5*l2^2*j3*l3 + 2*i1^3*i2*l2^2*j3*l3 + 2*i1*i2^2*l2^2*j3*l3 + 3*i1*i2*j2*l2^2*j3*l3 +
            i1*i2*k2*l2^2*j3*l3 + 4*i1*j2*k2*l2^2*j3*l3 + i1^3*l2^3*j3*l3 + 5*i1*i2*l2^3*j3*l3 + 2*i1*j2*l2^3*j3*l3 +
            6*i1*k2*l2^3*j3*l3 + 2*i1^6*j3^2*l3 + i1^4*i2*j3^2*l3 + 6*i2^3*j3^2*l3 + 4*i1^4*j2*j3^2*l3 +
            5*i1^2*i2*j2*j3^2*l3 + 2*i2^2*j2*j3^2*l3 + 4*i1^2*j2^2*j3^2*l3 + 6*i2*j2^2*j3^2*l3 + 5*i1^2*i2*k2*j3^2*l3
            + i2^2*k2*j3^2*l3 + 5*i1^2*j2*k2*j3^2*l3 + 2*i2*j2*k2*j3^2*l3 + 5*i2*k2^2*j3^2*l3 + 4*j2*k2^2*j3^2*l3 +
            2*k2^3*j3^2*l3 + i1^2*i2*l2*j3^2*l3 + 6*i2^2*l2*j3^2*l3 + 3*i1^2*j2*l2*j3^2*l3 + 6*i2*j2*l2*j3^2*l3 +
            6*j2^2*l2*j3^2*l3 + j2*k2*l2*j3^2*l3 + 2*k2^2*l2*j3^2*l3 + i1^2*l2^2*j3^2*l3 + 6*i2*l2^2*j3^2*l3 +
            j2*l2^2*j3^2*l3 + 6*k2*l2^2*j3^2*l3 + 5*l2^3*j3^2*l3 + 6*i1^3*j3^3*l3 + i1*j2*j3^3*l3 + 3*i1*k2*j3^3*l3 +
            2*i1*l2*j3^3*l3 + 6*i1^9*l3^2 + 5*i1^5*i2^2*l3^2 + 4*i1^3*i2^3*l3^2 + 6*i1*i2^4*l3^2 + 4*i1^7*j2*l3^2 +
            6*i1^5*i2*j2*l3^2 + 6*i1^3*i2^2*j2*l3^2 + 4*i1*i2^3*j2*l3^2 + 4*i1^5*j2^2*l3^2 + 3*i1*i2^2*j2^2*l3^2 +
            6*i1^3*j2^3*l3^2 + i1*j2^4*l3^2 + 2*i1^7*k2*l3^2 + 4*i1^5*i2*k2*l3^2 + 2*i1^3*i2^2*k2*l3^2 +
            5*i1*i2^3*k2*l3^2 + 3*i1^5*j2*k2*l3^2 + 3*i1^3*i2*j2*k2*l3^2 + 6*i1*i2^2*j2*k2*l3^2 + 2*i1*i2*j2^2*k2*l3^2
            + 2*i1*j2^3*k2*l3^2 + 2*i1^5*k2^2*l3^2 + 5*i1^3*i2*k2^2*l3^2 + 3*i1*i2^2*k2^2*l3^2 + 2*i1*i2*j2*k2^2*l3^2
            + 6*i1^3*k2^3*l3^2 + 2*i1*i2*k2^3*l3^2 + 2*i1*j2*k2^3*l3^2 + 5*i1*k2^4*l3^2 + 3*i1^7*l2*l3^2 +
            4*i1^3*i2^2*l2*l3^2 + 5*i1*i2^3*l2*l3^2 + 4*i1^5*j2*l2*l3^2 + i1^3*i2*j2*l2*l3^2 + 3*i1*i2^2*j2*l2*l3^2 +
            i1^3*j2^2*l2*l3^2 + 5*i1*i2*j2^2*l2*l3^2 + 3*i1*j2^3*l2*l3^2 + i1^5*k2*l2*l3^2 + 4*i1*i2^2*k2*l2*l3^2 +
            4*i1*i2*j2*k2*l2*l3^2 + 6*i1*j2^2*k2*l2*l3^2 + 4*i1^3*k2^2*l2*l3^2 + 4*i1*j2*k2^2*l2*l3^2 +
            6*i1*k2^3*l2*l3^2 + 4*i1^5*l2^2*l3^2 + 3*i1^3*i2*l2^2*l3^2 + 3*i1*i2*j2*l2^2*l3^2 + 4*i1^3*k2*l2^2*l3^2 +
            2*i1*i2*k2*l2^2*l3^2 + 2*i1*j2*k2*l2^2*l3^2 + 3*i1*k2^2*l2^2*l3^2 + i1^3*l2^3*l3^2 + i1*i2*l2^3*l3^2 +
            i1*j2*l2^3*l3^2 + 3*i1*k2*l2^3*l3^2 + 2*i1*l2^4*l3^2 + 4*i1^6*i3*l3^2 + 3*i1^2*i2^2*i3*l3^2 +
            i1^4*j2*i3*l3^2 + 3*i1^2*i2*j2*i3*l3^2 + 6*i2^2*j2*i3*l3^2 + 5*i1^2*j2^2*i3*l3^2 + 2*i2*j2^2*i3*l3^2 +
            6*j2^3*i3*l3^2 + 2*i1^4*k2*i3*l3^2 + 3*i1^2*i2*k2*i3*l3^2 + 4*i2^2*k2*i3*l3^2 + 3*i1^2*j2*k2*i3*l3^2 +
            3*i2*j2*k2*i3*l3^2 + 5*i1^2*k2^2*i3*l3^2 + 4*i2*k2^2*i3*l3^2 + j2*k2^2*i3*l3^2 + 4*i1^4*l2*i3*l3^2 +
            3*i1^2*i2*l2*i3*l3^2 + 5*i2^2*l2*i3*l3^2 + 4*i1^2*j2*l2*i3*l3^2 + 4*i2*j2*l2*i3*l3^2 +
            5*i1^2*k2*l2*i3*l3^2 + 3*i2*k2*l2*i3*l3^2 + 2*j2*k2*l2*i3*l3^2 + k2^2*l2*i3*l3^2 + 4*i1^2*l2^2*i3*l3^2 +
            j2*l2^2*i3*l3^2 + 4*k2*l2^2*i3*l3^2 + 5*l2^3*i3*l3^2 + 2*i1*i2*i3^2*l3^2 + i1*j2*i3^2*l3^2 +
            2*i1*l2*i3^2*l3^2 + 2*i1^6*j3*l3^2 + i1^4*i2*j3*l3^2 + 4*i1^2*i2^2*j3*l3^2 + 6*i2^3*j3*l3^2 +
            i1^4*j2*j3*l3^2 + 4*i2^2*j2*j3*l3^2 + 2*i1^2*j2^2*j3*l3^2 + 2*i2*j2^2*j3*l3^2 + 2*j2^3*j3*l3^2 +
            i1^4*k2*j3*l3^2 + i1^2*i2*k2*j3*l3^2 + 2*i2^2*k2*j3*l3^2 + 5*i1^2*j2*k2*j3*l3^2 + 6*i2*j2*k2*j3*l3^2 +
            6*i1^2*k2^2*j3*l3^2 + 3*i2*k2^2*j3*l3^2 + 2*i1^4*l2*j3*l3^2 + 4*i1^2*i2*l2*j3*l3^2 + i2^2*l2*j3*l3^2 +
            3*i1^2*j2*l2*j3*l3^2 + 2*j2^2*l2*j3*l3^2 + 5*i1^2*k2*l2*j3*l3^2 + 3*i2*k2*l2*j3*l3^2 + 4*j2*k2*l2*j3*l3^2
            + 2*k2^2*l2*j3*l3^2 + 5*i1^2*l2^2*j3*l3^2 + 2*j2*l2^2*j3*l3^2 + 5*k2*l2^2*j3*l3^2 + 4*l2^3*j3*l3^2 +
            3*i1*i2*j3^2*l3^2 + i1*j2*j3^2*l3^2 + 4*i1*k2*j3^2*l3^2 + 5*i1*l2*j3^2*l3^2 + 6*i1^6*l3^3 + 5*i1^4*i2*l3^3
            + 3*i1^2*i2^2*l3^3 + 2*i2^3*l3^3 + 6*i1^2*i2*j2*l3^3 + 6*i2^2*j2*l3^3 + 5*i1^2*j2^2*l3^3 + 2*i2*j2^2*l3^3
            + 4*j2^3*l3^3 + 5*i1^2*i2*k2*l3^3 + 4*i2^2*k2*l3^3 + 5*i1^2*j2*k2*l3^3 + 6*i2*j2*k2*l3^3 + 6*j2^2*k2*l3^3
            + 2*i1^2*k2^2*l3^3 + i2*k2^2*l3^3 + 6*j2*k2^2*l3^3 + 6*k2^3*l3^3 + 3*i1^4*l2*l3^3 + 3*i1^2*i2*l2*l3^3 +
            2*i2^2*l2*l3^3 + 2*j2^2*l2*l3^3 + 6*i1^2*k2*l2*l3^3 + 5*i2*k2*l2*l3^3 + 2*k2^2*l2*l3^3 + 5*i1^2*l2^2*l3^3
            + 3*i2*l2^2*l3^3 + 5*j2*l2^2*l3^3 + 5*k2*l2^2*l3^3 + 5*l2^3*l3^3 + 6*i1^3*i3*l3^3 + 5*i1*i2*i3*l3^3 +
            i1*j2*i3*l3^3 + 4*i1*k2*i3*l3^3 + 3*i1*l2*i3*l3^3 + 6*i3^2*l3^3 + 2*i1^3*j3*l3^3 + 2*i1*i2*j3*l3^3 +
            i1*j2*j3*l3^3 + 4*i1*k2*j3*l3^3 + 5*i1*l2*j3*l3^3 + 6*i1^3*l3^4 + 2*i1*i2*l3^4 + 3*i1*j2*l3^4 + 5*i1^11*i4
            + 3*i1^9*i2*i4 + 3*i1^7*i2^2*i4 + 5*i1^5*i2^3*i4 + 6*i1^3*i2^4*i4 + 4*i1*i2^5*i4 + 5*i1^7*i2*j2*i4 +
            6*i1^5*i2^2*j2*i4 + 4*i1^3*i2^3*j2*i4 + i1*i2^4*j2*i4 + 6*i1^7*j2^2*i4 + 5*i1^5*i2*j2^2*i4 +
            i1^3*i2^2*j2^2*i4 + 5*i1*i2^3*j2^2*i4 + 5*i1^5*j2^3*i4 + 4*i1^3*i2*j2^3*i4 + 2*i1*i2^2*j2^3*i4 +
            6*i1^3*j2^4*i4 + 6*i1*i2*j2^4*i4 + 3*i1*j2^5*i4 + 5*i1^9*k2*i4 + 4*i1^7*i2*k2*i4 + 4*i1*i2^4*k2*i4 +
            5*i1^7*j2*k2*i4 + 3*i1^5*i2*j2*k2*i4 + i1^3*i2^2*j2*k2*i4 + 6*i1^5*j2^2*k2*i4 + 6*i1^3*i2*j2^2*k2*i4 +
            4*i1*i2^2*j2^2*k2*i4 + 4*i1*i2*j2^3*k2*i4 + 2*i1*j2^4*k2*i4 + 2*i1^7*k2^2*i4 + 3*i1^5*i2*k2^2*i4 +
            5*i1^3*i2^2*k2^2*i4 + 2*i1*i2^3*k2^2*i4 + 3*i1^5*j2*k2^2*i4 + 6*i1^3*j2^2*k2^2*i4 + 4*i1*i2*j2^2*k2^2*i4 +
            i1*j2^3*k2^2*i4 + 2*i1^5*k2^3*i4 + 2*i1^3*i2*k2^3*i4 + 5*i1*i2^2*k2^3*i4 + 2*i1^3*j2*k2^3*i4 +
            2*i1*i2*j2*k2^3*i4 + 6*i1^3*k2^4*i4 + 2*i1*i2*k2^4*i4 + 6*i1*j2*k2^4*i4 + 6*i1*k2^5*i4 + 5*i1^9*l2*i4 +
            2*i1^7*i2*l2*i4 + 2*i1^5*i2^2*l2*i4 + i1^3*i2^3*l2*i4 + 6*i1*i2^4*l2*i4 + i1^7*j2*l2*i4 + i1^5*i2*j2*l2*i4
            + 2*i1^3*i2^2*j2*l2*i4 + 2*i1*i2^3*j2*l2*i4 + 3*i1^5*j2^2*l2*i4 + 5*i1^3*i2*j2^2*l2*i4 + 6*i1^3*j2^3*l2*i4
            + 5*i1*i2*j2^3*l2*i4 + i1*j2^4*l2*i4 + 3*i1^7*k2*l2*i4 + 3*i1^5*i2*k2*l2*i4 + 5*i1^3*i2^2*k2*l2*i4 +
            2*i1^5*j2*k2*l2*i4 + 4*i1^3*i2*j2*k2*l2*i4 + 2*i1^3*j2^2*k2*l2*i4 + i1*i2*j2^2*k2*l2*i4 +
            6*i1*j2^3*k2*l2*i4 + 3*i1^5*k2^2*l2*i4 + i1^3*i2*k2^2*l2*i4 + 2*i1*i2^2*k2^2*l2*i4 + 4*i1^3*j2*k2^2*l2*i4
            + i1*i2*j2*k2^2*l2*i4 + 2*i1*j2^2*k2^2*l2*i4 + 3*i1^3*k2^3*l2*i4 + i1*i2*k2^3*l2*i4 + 4*i1*k2^4*l2*i4 +
            5*i1^7*l2^2*i4 + 5*i1^5*i2*l2^2*i4 + 4*i1^3*i2^2*l2^2*i4 + 3*i1*i2^3*l2^2*i4 + 5*i1^5*j2*l2^2*i4 +
            i1*i2^2*j2*l2^2*i4 + 3*i1^3*j2^2*l2^2*i4 + 3*i1*j2^3*l2^2*i4 + 2*i1^5*k2*l2^2*i4 + 5*i1^3*i2*k2*l2^2*i4 +
            i1*i2^2*k2*l2^2*i4 + 4*i1^3*j2*k2*l2^2*i4 + 2*i1*i2*j2*k2*l2^2*i4 + 4*i1*j2^2*k2*l2^2*i4 +
            3*i1^3*k2^2*l2^2*i4 + 2*i1*j2*k2^2*l2^2*i4 + 5*i1^5*l2^3*i4 + 5*i1^3*i2*l2^3*i4 + 6*i1*i2*j2*l2^3*i4 +
            6*i1*j2^2*l2^3*i4 + 4*i1^3*k2*l2^3*i4 + i1*i2*k2*l2^3*i4 + 4*i1*j2*k2*l2^3*i4 + 4*i1*k2^2*l2^3*i4 +
            2*i1^3*l2^4*i4 + 3*i1*j2*l2^4*i4 + 3*i1*k2*l2^4*i4 + 6*i1*l2^5*i4 + 5*i1^4*i2^2*i3*i4 + i1^2*i2^3*i3*i4 +
            3*i1^6*j2*i3*i4 + 5*i1^4*i2*j2*i3*i4 + 4*i1^2*i2^2*j2*i3*i4 + 4*i1^4*j2^2*i3*i4 + 3*i1^2*i2*j2^2*i3*i4 +
            6*i1^2*j2^3*i3*i4 + 3*i1^6*k2*i3*i4 + 6*i1^4*i2*k2*i3*i4 + i1^2*i2^2*k2*i3*i4 + i2^3*k2*i3*i4 +
            6*i1^4*j2*k2*i3*i4 + i1^2*i2*j2*k2*i3*i4 + 4*i2^2*j2*k2*i3*i4 + 5*i1^2*j2^2*k2*i3*i4 + 3*i2*j2^2*k2*i3*i4
            + 6*j2^3*k2*i3*i4 + 3*i1^4*k2^2*i3*i4 + 4*i1^2*i2*k2^2*i3*i4 + 4*i2^2*k2^2*i3*i4 + 4*i1^2*j2*k2^2*i3*i4 +
            3*j2^2*k2^2*i3*i4 + i2*k2^3*i3*i4 + 4*k2^4*i3*i4 + 6*i1^4*i2*l2*i3*i4 + 2*i1^2*i2^2*l2*i3*i4 +
            i2^3*l2*i3*i4 + 4*i1^4*j2*l2*i3*i4 + 2*i1^2*i2*j2*l2*i3*i4 + 4*i2^2*j2*l2*i3*i4 + 3*i1^2*j2^2*l2*i3*i4 +
            3*i2*j2^2*l2*i3*i4 + 6*j2^3*l2*i3*i4 + i1^4*k2*l2*i3*i4 + 2*i1^2*i2*k2*l2*i3*i4 + 2*i2^2*k2*l2*i3*i4 +
            3*i1^2*j2*k2*l2*i3*i4 + 5*i2*j2*k2*l2*i3*i4 + 4*i1^2*k2^2*l2*i3*i4 + 4*i2*k2^2*l2*i3*i4 +
            5*j2*k2^2*l2*i3*i4 + 2*k2^3*l2*i3*i4 + i1^2*i2*l2^2*i3*i4 + i2^2*l2^2*i3*i4 + 3*i1^2*j2*l2^2*i3*i4 +
            i2*j2*l2^2*i3*i4 + 5*j2^2*l2^2*i3*i4 + 6*i1^2*k2*l2^2*i3*i4 + 5*i2*k2*l2^2*i3*i4 + 3*j2*k2*l2^2*i3*i4 +
            5*k2^2*l2^2*i3*i4 + i1^2*l2^3*i3*i4 + 2*i2*l2^3*i3*i4 + 5*j2*l2^3*i3*i4 + 6*k2*l2^3*i3*i4 + 6*l2^4*i3*i4 +
            2*i1^5*i3^2*i4 + 6*i1^3*i2*i3^2*i4 + i1^3*j2*i3^2*i4 + 6*i1^3*k2*i3^2*i4 + 6*i1*i2*k2*i3^2*i4 +
            i1*j2*k2*i3^2*i4 + 2*i1^3*l2*i3^2*i4 + i1*i2*l2*i3^2*i4 + 6*i1*j2*l2*i3^2*i4 + 3*i1*k2*l2*i3^2*i4 +
            4*i1*l2^2*i3^2*i4 + i1^8*j3*i4 + 4*i1^6*i2*j3*i4 + 3*i1^2*i2^3*j3*i4 + 4*i2^4*j3*i4 + i1^6*j2*j3*i4 +
            4*i1^4*i2*j2*j3*i4 + 5*i2^3*j2*j3*i4 + 5*i1^2*i2*j2^2*j3*i4 + 3*i2^2*j2^2*j3*i4 + 6*i1^2*j2^3*j3*i4 +
            5*i2*j2^3*j3*i4 + 4*j2^4*j3*i4 + i1^4*i2*k2*j3*i4 + 2*i1^2*i2^2*k2*j3*i4 + 3*i2^3*k2*j3*i4 +
            i1^4*j2*k2*j3*i4 + 6*i1^2*i2*j2*k2*j3*i4 + 4*i1^2*j2^2*k2*j3*i4 + 5*i2*j2^2*k2*j3*i4 + 6*j2^3*k2*j3*i4 +
            6*i1^4*k2^2*j3*i4 + 4*i1^2*i2*k2^2*j3*i4 + 6*i2^2*k2^2*j3*i4 + 4*i1^2*j2*k2^2*j3*i4 + 4*i2*j2*k2^2*j3*i4 +
            6*j2^2*k2^2*j3*i4 + 2*i1^2*k2^3*j3*i4 + 6*i2*k2^3*j3*i4 + 6*j2*k2^3*j3*i4 + i1^6*l2*j3*i4 +
            5*i1^4*i2*l2*j3*i4 + i1^2*i2^2*l2*j3*i4 + 3*i2^3*l2*j3*i4 + 5*i1^4*j2*l2*j3*i4 + 5*i1^2*i2*j2*l2*j3*i4 +
            i2^2*j2*l2*j3*i4 + 6*i1^2*j2^2*l2*j3*i4 + 3*i2*j2^2*l2*j3*i4 + i1^4*k2*l2*j3*i4 + 6*i1^2*i2*k2*l2*j3*i4 +
            2*i2^2*k2*l2*j3*i4 + i1^2*j2*k2*l2*j3*i4 + i2*j2*k2*l2*j3*i4 + 2*j2^2*k2*l2*j3*i4 + 5*i2*k2^2*l2*j3*i4 +
            6*j2*k2^2*l2*j3*i4 + 4*k2^3*l2*j3*i4 + 5*i1^4*l2^2*j3*i4 + 5*i1^2*i2*l2^2*j3*i4 + 2*i1^2*j2*l2^2*j3*i4 +
            6*i2*j2*l2^2*j3*i4 + 5*j2^2*l2^2*j3*i4 + i1^2*k2*l2^2*j3*i4 + 3*k2^2*l2^2*j3*i4 + i1^2*l2^3*j3*i4 +
            j2*l2^3*j3*i4 + 4*k2*l2^3*j3*i4 + 4*l2^4*j3*i4 + 5*i1^5*i3*j3*i4 + 3*i1^3*i2*i3*j3*i4 + 3*i1*i2^2*i3*j3*i4
            + 6*i1^3*j2*i3*j3*i4 + i1*i2*j2*i3*j3*i4 + 3*i1*j2^2*i3*j3*i4 + 6*i1^3*k2*i3*j3*i4 + 5*i1*i2*k2*i3*j3*i4 +
            4*i1*j2*k2*i3*j3*i4 + 4*i1*k2^2*i3*j3*i4 + 4*i1*i2*l2*i3*j3*i4 + 5*i1*j2*l2*i3*j3*i4 + 2*i1*l2^2*i3*j3*i4
            + i1^2*i3^2*j3*i4 + 4*l2*i3^2*j3*i4 + 6*i1^5*j3^2*i4 + 4*i1^3*i2*j3^2*i4 + 4*i1*i2^2*j3^2*i4 +
            i1^3*j2*j3^2*i4 + 4*i1*i2*j2*j3^2*i4 + 6*i1*j2^2*j3^2*i4 + 3*i1^3*k2*j3^2*i4 + 3*i1*i2*k2*j3^2*i4 +
            5*i1*k2^2*j3^2*i4 + 3*i1^3*l2*j3^2*i4 + 3*i1*j2*l2*j3^2*i4 + 3*i1*k2*l2*j3^2*i4 + 6*i1*l2^2*j3^2*i4 +
            i1^2*i3*j3^2*i4 + 5*i2*i3*j3^2*i4 + 2*j2*i3*j3^2*i4 + 2*k2*i3*j3^2*i4 + l2*i3*j3^2*i4 + 5*i1^2*j3^3*i4 +
            4*i2*j3^3*i4 + 3*j2*j3^3*i4 + 6*k2*j3^3*i4 + 4*l2*j3^3*i4 + 4*i1^8*k3*i4 + 6*i1^6*i2*k3*i4 +
            2*i1^4*i2^2*k3*i4 + i1^2*i2^3*k3*i4 + i2^4*k3*i4 + 6*i1^6*j2*k3*i4 + 3*i1^4*i2*j2*k3*i4 +
            6*i1^2*i2^2*j2*k3*i4 + 3*i2^3*j2*k3*i4 + 2*i1^4*j2^2*k3*i4 + 6*i1^2*i2*j2^2*k3*i4 + 6*i2^2*j2^2*k3*i4 +
            i1^2*j2^3*k3*i4 + 3*i2*j2^3*k3*i4 + j2^4*k3*i4 + 5*i1^6*k2*k3*i4 + 4*i1^4*i2*k2*k3*i4 + 2*i2^3*k2*k3*i4 +
            6*i1^2*i2*j2*k2*k3*i4 + 6*i2^2*j2*k2*k3*i4 + 2*i1^2*j2^2*k2*k3*i4 + 3*i2*j2^2*k2*k3*i4 + 3*j2^3*k2*k3*i4 +
            6*i1^4*k2^2*k3*i4 + 6*i1^2*i2*k2^2*k3*i4 + 5*i2^2*k2^2*k3*i4 + 3*i1^2*j2*k2^2*k3*i4 + 2*i2*j2*k2^2*k3*i4 +
            4*j2^2*k2^2*k3*i4 + 3*i1^2*k2^3*k3*i4 + 4*i2*k2^3*k3*i4 + 5*j2*k2^3*k3*i4 + 2*k2^4*k3*i4 +
            3*i1^4*i2*l2*k3*i4 + 2*i1^2*i2^2*l2*k3*i4 + 6*i2^3*l2*k3*i4 + 5*i1^4*j2*l2*k3*i4 + 3*i1^2*i2*j2*l2*k3*i4 +
            6*i2^2*j2*l2*k3*i4 + 2*i1^2*j2^2*l2*k3*i4 + 5*i2*j2^2*l2*k3*i4 + 4*j2^3*l2*k3*i4 + 6*i1^4*k2*l2*k3*i4 +
            4*i2^2*k2*l2*k3*i4 + 3*i1^2*j2*k2*l2*k3*i4 + 3*i2*j2*k2*l2*k3*i4 + 2*j2^2*k2*l2*k3*i4 +
            5*i1^2*k2^2*l2*k3*i4 + i2*k2^2*l2*k3*i4 + 4*j2*k2^2*l2*k3*i4 + 4*k2^3*l2*k3*i4 + 4*i1^4*l2^2*k3*i4 +
            6*i1^2*i2*l2^2*k3*i4 + i2^2*l2^2*k3*i4 + i1^2*j2*l2^2*k3*i4 + 6*i2*j2*l2^2*k3*i4 + 3*i1^2*k2*l2^2*k3*i4 +
            4*i2*k2*l2^2*k3*i4 + 5*j2*k2*l2^2*k3*i4 + 3*k2^2*l2^2*k3*i4 + 5*i1^2*l2^3*k3*i4 + 5*i2*l2^3*k3*i4 +
            3*j2*l2^3*k3*i4 + 5*k2*l2^3*k3*i4 + 2*i1^5*j3*k3*i4 + 4*i1^3*i2*j3*k3*i4 + 4*i1*i2^2*j3*k3*i4 +
            2*i1^3*j2*j3*k3*i4 + 5*i1*i2*j2*j3*k3*i4 + 5*i1*j2^2*j3*k3*i4 + 4*i1^3*k2*j3*k3*i4 + 2*i1*i2*k2*j3*k3*i4 +
            2*i1*k2^2*j3*k3*i4 + 6*i1^3*l2*j3*k3*i4 + i1*i2*l2*j3*k3*i4 + 5*i1*j2*l2*j3*k3*i4 + 5*i1*k2*l2*j3*k3*i4 +
            i1*l2^2*j3*k3*i4 + 5*i1^2*j3^2*k3*i4 + 6*i2*j3^2*k3*i4 + 2*j2*j3^2*k3*i4 + 3*k2*j3^2*k3*i4 +
            3*l2*j3^2*k3*i4 + 5*i1^11*j4 + 4*i1^9*i2*j4 + 4*i1^7*i2^2*j4 + 5*i1^5*i2^3*j4 + 2*i1^3*i2^4*j4 +
            6*i1^9*j2*j4 + 3*i1^7*i2*j2*j4 + 6*i1^3*i2^3*j2*j4 + 4*i1^7*j2^2*j4 + 6*i1^5*i2*j2^2*j4 +
            5*i1^3*i2^2*j2^2*j4 + 3*i1^5*j2^3*j4 + 6*i1^3*i2*j2^3*j4 + 2*i1^3*j2^4*j4 + 2*i1^9*k2*j4 +
            2*i1^5*i2^2*k2*j4 + 5*i1^3*i2^3*k2*j4 + 5*i1*i2^4*k2*j4 + i1^7*j2*k2*j4 + 3*i1^5*i2*j2*k2*j4 +
            i1^3*i2^2*j2*k2*j4 + i1*i2^3*j2*k2*j4 + 3*i1^5*j2^2*k2*j4 + 4*i1^3*i2*j2^2*k2*j4 + 2*i1*i2^2*j2^2*k2*j4 +
            4*i1^3*j2^3*k2*j4 + i1*i2*j2^3*k2*j4 + 5*i1*j2^4*k2*j4 + 3*i1^7*k2^2*j4 + i1^5*i2*k2^2*j4 +
            5*i1^3*i2^2*k2^2*j4 + 6*i1*i2^3*k2^2*j4 + 3*i1^5*j2*k2^2*j4 + 6*i1^3*i2*j2*k2^2*j4 + 4*i1*i2^2*j2*k2^2*j4
            + 4*i1^3*j2^2*k2^2*j4 + 2*i1*i2*j2^2*k2^2*j4 + 2*i1*j2^3*k2^2*j4 + 5*i1^5*k2^3*j4 + 6*i1^3*i2*k2^3*j4 +
            4*i1*i2^2*k2^3*j4 + 4*i1^3*j2*k2^3*j4 + i1*i2*j2*k2^3*j4 + 4*i1*j2^2*k2^3*j4 + 2*i1^3*k2^4*j4 +
            i1*i2*k2^4*j4 + 4*i1*j2*k2^4*j4 + 2*i1*k2^5*j4 + i1^9*l2*j4 + 5*i1^7*i2*l2*j4 + 6*i1^5*i2^2*l2*j4 +
            3*i1^3*i2^3*l2*j4 + i1*i2^4*l2*j4 + 2*i1^7*j2*l2*j4 + i1^5*i2*j2*l2*j4 + 2*i1^3*i2^2*j2*l2*j4 +
            3*i1*i2^3*j2*l2*j4 + 5*i1^5*j2^2*l2*j4 + i1^3*i2*j2^2*l2*j4 + 6*i1*i2^2*j2^2*l2*j4 + i1^3*j2^3*l2*j4 +
            3*i1*i2*j2^3*l2*j4 + i1*j2^4*l2*j4 + i1^7*k2*l2*j4 + 4*i1^5*i2*k2*l2*j4 + 2*i1^3*i2^2*k2*l2*j4 +
            3*i1*i2^3*k2*l2*j4 + i1^5*j2*k2*l2*j4 + 2*i1^3*i2*j2*k2*l2*j4 + 2*i1*i2^2*j2*k2*l2*j4 +
            4*i1^3*j2^2*k2*l2*j4 + i1*i2*j2^2*k2*l2*j4 + i1*j2^3*k2*l2*j4 + i1^5*k2^2*l2*j4 + 5*i1*i2^2*k2^2*l2*j4 +
            5*i1^3*j2*k2^2*l2*j4 + 6*i1*i2*j2*k2^2*l2*j4 + 4*i1*i2*k2^3*l2*j4 + i1*j2*k2^3*l2*j4 + 2*i1^7*l2^2*j4 +
            i1^5*i2*l2^2*j4 + 6*i1*i2^3*l2^2*j4 + 2*i1^5*j2*l2^2*j4 + 4*i1*i2^2*j2*l2^2*j4 + 4*i1^3*j2^2*l2^2*j4 +
            2*i1*i2*j2^2*l2^2*j4 + 2*i1*j2^3*l2^2*j4 + 2*i1^5*k2*l2^2*j4 + 5*i1^3*i2*k2*l2^2*j4 + 3*i1*i2^2*k2*l2^2*j4
            + 2*i1^3*j2*k2*l2^2*j4 + 4*i1*i2*j2*k2*l2^2*j4 + 5*i1*j2^2*k2*l2^2*j4 + 2*i1^3*k2^2*l2^2*j4 +
            2*i1*j2*k2^2*l2^2*j4 + 2*i1*k2^3*l2^2*j4 + 5*i1^5*l2^3*j4 + 6*i1^3*i2*l2^3*j4 + 2*i1*i2^2*l2^3*j4 +
            5*i1^3*j2*l2^3*j4 + i1*i2*j2*l2^3*j4 + 5*i1*j2^2*l2^3*j4 + 6*i1*i2*k2*l2^3*j4 + 2*i1*k2^2*l2^3*j4 +
            5*i1^3*l2^4*j4 + 2*i1*i2*l2^4*j4 + 5*i1*j2*l2^4*j4 + 4*i1*k2*l2^4*j4 + 5*i1*l2^5*j4 + i1^8*j3*j4 +
            4*i1^4*i2^2*j3*j4 + 3*i1^2*i2^3*j3*j4 + 3*i1^6*j2*j3*j4 + i1^4*i2*j2*j3*j4 + 5*i1^2*i2^2*j2*j3*j4 +
            2*i1^4*j2^2*j3*j4 + 2*i1^2*i2*j2^2*j3*j4 + 4*i1^2*j2^3*j3*j4 + 3*i1^2*i2^2*k2*j3*j4 + 3*i2^3*k2*j3*j4 +
            2*i1^4*j2*k2*j3*j4 + 6*i1^2*i2*j2*k2*j3*j4 + 5*i2^2*j2*k2*j3*j4 + 5*i1^2*j2^2*k2*j3*j4 +
            2*i2*j2^2*k2*j3*j4 + 4*j2^3*k2*j3*j4 + 4*i1^2*i2*k2^2*j3*j4 + 3*i2^2*k2^2*j3*j4 + 3*i1^2*j2*k2^2*j3*j4 +
            3*i2*j2*k2^2*j3*j4 + j2^2*k2^2*j3*j4 + i1^2*k2^3*j3*j4 + 5*i2*k2^3*j3*j4 + 6*j2*k2^3*j3*j4 +
            2*i1^4*i2*l2*j3*j4 + 3*i1^2*i2^2*l2*j3*j4 + 4*i2^3*l2*j3*j4 + 4*i1^4*j2*l2*j3*j4 + 4*i1^2*i2*j2*l2*j3*j4 +
            2*i2^2*j2*l2*j3*j4 + 5*i2*j2^2*l2*j3*j4 + 3*j2^3*l2*j3*j4 + 3*i1^4*k2*l2*j3*j4 + 5*i1^2*i2*k2*l2*j3*j4 +
            i2^2*k2*l2*j3*j4 + 6*i1^2*j2*k2*l2*j3*j4 + i2*j2*k2*l2*j3*j4 + 5*j2^2*k2*l2*j3*j4 + 6*i1^2*k2^2*l2*j3*j4 +
            6*i2*k2^2*l2*j3*j4 + 4*j2*k2^2*l2*j3*j4 + 3*i1^4*l2^2*j3*j4 + i1^2*i2*l2^2*j3*j4 + 6*i2^2*l2^2*j3*j4 +
            6*i1^2*j2*l2^2*j3*j4 + j2^2*l2^2*j3*j4 + 3*i1^2*k2*l2^2*j3*j4 + 6*i2*k2*l2^2*j3*j4 + 6*j2*k2*l2^2*j3*j4 +
            3*k2^2*l2^2*j3*j4 + 2*i1^2*l2^3*j3*j4 + 5*i2*l2^3*j3*j4 + 6*j2*l2^3*j3*j4 + 3*k2*l2^3*j3*j4 + l2^4*j3*j4 +
            5*i1^5*j3^2*j4 + 3*i1^3*i2*j3^2*j4 + 3*i1*i2^2*j3^2*j4 + i1*i2*j2*j3^2*j4 + 3*i1*j2^2*j3^2*j4 +
            2*i1^3*k2*j3^2*j4 + 4*i1*i2*k2*j3^2*j4 + i1*j2*k2*j3^2*j4 + 6*i1*k2^2*j3^2*j4 + 4*i1^3*l2*j3^2*j4 +
            5*i1*i2*l2*j3^2*j4 + 2*i1*j2*l2*j3^2*j4 + 2*i1*k2*l2*j3^2*j4 + 5*i1^2*j3^3*j4 + 6*i2*j3^3*j4 + j2*j3^3*j4
            + 5*k2*j3^3*j4 + l2*j3^3*j4 + 3*i1^11*k4 + i1^7*i2^2*k4 + 2*i1^5*i2^3*k4 + 3*i1^9*j2*k4 + 3*i1^7*i2*j2*k4
            + i1^5*i2^2*j2*k4 + 3*i1^7*j2^2*k4 + 6*i1^5*i2*j2^2*k4 + 5*i1^5*j2^3*k4 + 6*i1^9*k2*k4 + 3*i1^7*i2*k2*k4 +
            2*i1^5*i2^2*k2*k4 + 4*i1^3*i2^3*k2*k4 + 5*i1^7*j2*k2*k4 + 4*i1^5*i2*j2*k2*k4 + 2*i1^3*i2^2*j2*k2*k4 +
            i1^5*j2^2*k2*k4 + 5*i1^3*i2*j2^2*k2*k4 + 3*i1^3*j2^3*k2*k4 + 6*i1^7*k2^2*k4 + 3*i1^5*i2*k2^2*k4 +
            i1*i2^3*k2^2*k4 + 3*i1^5*j2*k2^2*k4 + 3*i1^3*i2*j2*k2^2*k4 + 4*i1*i2^2*j2*k2^2*k4 + 4*i1^3*j2^2*k2^2*k4 +
            3*i1*i2*j2^2*k2^2*k4 + 6*i1*j2^3*k2^2*k4 + 4*i1^5*k2^3*k4 + i1^3*i2*k2^3*k4 + 2*i1*i2^2*k2^3*k4 +
            4*i1^3*j2*k2^3*k4 + 5*i1*i2*j2*k2^3*k4 + 4*i1^3*k2^4*k4 + i1*i2*k2^4*k4 + 3*i1*j2*k2^4*k4 + 3*i1*k2^5*k4 +
            4*i1^7*i2*l2*k4 + 6*i1^5*i2^2*l2*k4 + 5*i1^7*j2*l2*k4 + 2*i1^5*i2*j2*l2*k4 + 6*i1^5*j2^2*l2*k4 +
            2*i1^7*k2*l2*k4 + 5*i1^5*i2*k2*l2*k4 + 5*i1^3*i2^2*k2*l2*k4 + 2*i1*i2^3*k2*l2*k4 + 2*i1^5*j2*k2*l2*k4 +
            i1*i2^2*j2*k2*l2*k4 + 2*i1^3*j2^2*k2*l2*k4 + 6*i1*i2*j2^2*k2*l2*k4 + 5*i1*j2^3*k2*l2*k4 + i1^5*k2^2*l2*k4
            + i1*i2^2*k2^2*l2*k4 + 5*i1^3*j2*k2^2*l2*k4 + i1*i2*j2*k2^2*l2*k4 + 5*i1*j2^2*k2^2*l2*k4 +
            4*i1^3*k2^3*l2*k4 + 5*i1*i2*k2^3*l2*k4 + 2*i1*j2*k2^3*l2*k4 + 3*i1*k2^4*l2*k4 + 5*i1^7*l2^2*k4 +
            5*i1^5*i2*l2^2*k4 + 3*i1*i2^3*l2^2*k4 + 3*i1^5*j2*l2^2*k4 + 3*i1^3*i2*j2*l2^2*k4 + 5*i1*i2^2*j2*l2^2*k4 +
            4*i1^3*j2^2*l2^2*k4 + 2*i1*i2*j2^2*l2^2*k4 + 4*i1*j2^3*l2^2*k4 + 4*i1^5*k2*l2^2*k4 + i1^3*i2*k2*l2^2*k4 +
            5*i1*i2^2*k2*l2^2*k4 + 4*i1^3*j2*k2*l2^2*k4 + i1*i2*j2*k2*l2^2*k4 + i1*j2^2*k2*l2^2*k4 +
            3*i1*i2*k2^2*l2^2*k4 + 3*i1*j2*k2^2*l2^2*k4 + i1*k2^3*l2^2*k4 + i1^3*i2*l2^3*k4 + 4*i1*i2^2*l2^3*k4 +
            2*i1*i2*j2*l2^3*k4 + i1*j2^2*l2^3*k4 + 4*i1^3*k2*l2^3*k4 + 6*i1*i2*k2*l2^3*k4 + i1*j2*k2*l2^3*k4 +
            6*i1*k2^2*l2^3*k4 + i1^3*l2^4*k4 + 5*i1*i2*l2^4*k4 + 5*i1*j2*l2^4*k4 + i1*k2*l2^4*k4 + 3*i1*l2^5*k4 +
            i1^8*i3*k4 + 4*i1^6*i2*i3*k4 + 3*i1^6*j2*i3*k4 + 4*i1^4*i2*k2*i3*k4 + 3*i1^4*j2*k2*i3*k4 +
            3*i1^4*k2^2*i3*k4 + 4*i1^2*i2*k2^2*i3*k4 + 3*i1^2*j2*k2^2*i3*k4 + 2*i2*k2^3*i3*k4 + 5*j2*k2^3*i3*k4 +
            5*k2^4*i3*k4 + 6*i1^6*l2*i3*k4 + i1^4*i2*l2*i3*k4 + 6*i1^4*j2*l2*i3*k4 + i1^4*k2*l2*i3*k4 +
            6*i1^2*i2*k2*l2*i3*k4 + i1^2*j2*k2*l2*i3*k4 + 5*i1^2*k2^2*l2*i3*k4 + i2*k2^2*l2*i3*k4 + 6*j2*k2^2*l2*i3*k4
            + 6*k2^3*l2*i3*k4 + 4*i1^4*l2^2*i3*k4 + 6*i1^2*i2*l2^2*i3*k4 + i1^2*j2*l2^2*i3*k4 + i1^2*k2*l2^2*i3*k4 +
            i2*k2*l2^2*i3*k4 + 6*j2*k2*l2^2*i3*k4 + 6*k2^2*l2^2*i3*k4 + 3*i1^2*l2^3*i3*k4 + i2*l2^3*i3*k4 +
            6*j2*l2^3*i3*k4 + k2*l2^3*i3*k4 + 6*l2^4*i3*k4 + 2*i1^8*k3*k4 + 4*i1^6*i2*k3*k4 + 3*i1^4*i2^2*k3*k4 +
            2*i1^6*j2*k3*k4 + i1^4*i2*j2*k3*k4 + 3*i1^4*j2^2*k3*k4 + 5*i1^6*k2*k3*k4 + 4*i1^4*i2*k2*k3*k4 +
            6*i1^2*i2^2*k2*k3*k4 + 2*i1^4*j2*k2*k3*k4 + 2*i1^2*i2*j2*k2*k3*k4 + 6*i1^2*j2^2*k2*k3*k4 +
            4*i1^4*k2^2*k3*k4 + 2*i2^2*k2^2*k3*k4 + 2*i1^2*j2*k2^2*k3*k4 + 3*i2*j2*k2^2*k3*k4 + 2*j2^2*k2^2*k3*k4 +
            4*i1^2*k2^3*k3*k4 + 2*i2*k2^3*k3*k4 + 3*j2*k2^3*k3*k4 + k2^4*k3*k4 + 4*i1^6*l2*k3*k4 + i1^4*i2*l2*k3*k4 +
            3*i1^2*i2^2*l2*k3*k4 + 5*i1^4*j2*l2*k3*k4 + i1^2*i2*j2*l2*k3*k4 + 3*i1^2*j2^2*l2*k3*k4 +
            2*i1^4*k2*l2*k3*k4 + 2*i1^2*i2*k2*l2*k3*k4 + i2^2*k2*l2*k3*k4 + 5*i1^2*j2*k2*l2*k3*k4 +
            5*i2*j2*k2*l2*k3*k4 + j2^2*k2*l2*k3*k4 + i1^2*k2^2*l2*k3*k4 + i2*k2^2*l2*k3*k4 + 6*j2*k2^2*l2*k3*k4 +
            k2^3*l2*k3*k4 + 4*i1^4*l2^2*k3*k4 + 5*i1^2*i2*l2^2*k3*k4 + 3*i2^2*l2^2*k3*k4 + 2*i1^2*j2*l2^2*k3*k4 +
            i2*j2*l2^2*k3*k4 + 3*j2^2*l2^2*k3*k4 + 4*i1^2*k2*l2^2*k3*k4 + 2*i2*k2*l2^2*k3*k4 + 5*j2*k2*l2^2*k3*k4 +
            3*k2^2*l2^2*k3*k4 + i2*l2^3*k3*k4 + 6*j2*l2^3*k3*k4 + 2*k2*l2^3*k3*k4 + 4*l2^4*k3*k4 + i1^5*k3^2*k4 +
            6*i1^3*i2*k3^2*k4 + i1^3*j2*k3^2*k4 + 3*i1^3*k2*k3^2*k4 + 2*i1*i2*k2*k3^2*k4 + 5*i1*j2*k2*k3^2*k4 +
            i1*k2^2*k3^2*k4 + 6*i1^3*l2*k3^2*k4 + 3*i1*i2*l2*k3^2*k4 + 4*i1*j2*l2*k3^2*k4 + 3*i1*k2*l2*k3^2*k4 +
            2*i1*l2^2*k3^2*k4 + 2*i1^2*k3^3*k4 + 3*k2*k3^3*k4 + l2*k3^3*k4 + 4*i1^7*k4^2 + 4*i1^5*k2*k4^2 +
            i1^3*k2^2*k4^2 + 5*i1*k2^3*k4^2 + 2*i1^5*l2*k4^2 + i1^3*k2*l2*k4^2 + 6*i1*k2^2*l2*k4^2 + 3*i1^3*l2^2*k4^2
            + 4*i1*k2*l2^2*k4^2 + 5*i1*l2^3*k4^2 + 2*i1^11*l4 + i1^7*i2^2*l4 + 2*i1^5*i2^3*l4 + 5*i1^3*i2^4*l4 +
            6*i1*i2^5*l4 + 6*i1^9*j2*l4 + 6*i1^7*i2*j2*l4 + 6*i1^5*i2^2*j2*l4 + 4*i1^3*i2^3*j2*l4 + 5*i1*i2^4*j2*l4 +
            3*i1^7*j2^2*l4 + 3*i1^5*i2*j2^2*l4 + 4*i1*i2^3*j2^2*l4 + 3*i1^5*j2^3*l4 + 3*i1^3*i2*j2^3*l4 +
            3*i1*i2^2*j2^3*l4 + 2*i1^3*j2^4*l4 + 2*i1*i2*j2^4*l4 + i1*j2^5*l4 + 2*i1^9*k2*l4 + 5*i1^7*i2*k2*l4 +
            5*i1^5*i2^2*k2*l4 + 4*i1^3*i2^3*k2*l4 + 2*i1*i2^4*k2*l4 + 5*i1^7*j2*k2*l4 + 3*i1^5*i2*j2*k2*l4 +
            3*i1^3*i2^2*j2*k2*l4 + 6*i1^5*j2^2*k2*l4 + 2*i1^3*i2*j2^2*k2*l4 + 2*i1*i2^2*j2^2*k2*l4 + 5*i1^3*j2^3*k2*l4
            + 2*i1*i2*j2^3*k2*l4 + i1*j2^4*k2*l4 + i1^7*k2^2*l4 + 3*i1^5*i2*k2^2*l4 + 3*i1^3*i2^2*k2^2*l4 +
            5*i1*i2^3*k2^2*l4 + 2*i1^5*j2*k2^2*l4 + 6*i1^3*i2*j2*k2^2*l4 + 6*i1*i2^2*j2*k2^2*l4 + 3*i1^3*j2^2*k2^2*l4
            + 3*i1*i2*j2^2*k2^2*l4 + 5*i1^5*k2^3*l4 + 5*i1^3*i2*k2^3*l4 + 4*i1*i2^2*k2^3*l4 + 6*i1^3*j2*k2^3*l4 +
            i1*i2*j2*k2^3*l4 + 2*i1*j2^2*k2^3*l4 + 6*i1^3*k2^4*l4 + 4*i1*i2*k2^4*l4 + i1*j2*k2^4*l4 + i1*k2^5*l4 +
            2*i1^7*i2*l2*l4 + 3*i1^5*i2^2*l2*l4 + i1^3*i2^3*l2*l4 + 5*i1*i2^4*l2*l4 + 4*i1^7*j2*l2*l4 +
            5*i1^5*i2*j2*l2*l4 + 5*i1^5*j2^2*l2*l4 + 5*i1*i2^2*j2^2*l2*l4 + 6*i1^3*j2^3*l2*l4 + 5*i1*i2*j2^3*l2*l4 +
            6*i1*j2^4*l2*l4 + 5*i1^7*k2*l2*l4 + 6*i1^5*i2*k2*l2*l4 + 6*i1^3*i2^2*k2*l2*l4 + i1*i2^3*k2*l2*l4 +
            i1^5*j2*k2*l2*l4 + 6*i1*i2^2*j2*k2*l2*l4 + 3*i1^3*j2^2*k2*l2*l4 + i1*i2*j2^2*k2*l2*l4 + 6*i1*j2^3*k2*l2*l4
            + 2*i1^5*k2^2*l2*l4 + 2*i1^3*i2*k2^2*l2*l4 + 5*i1*i2^2*k2^2*l2*l4 + 3*i1^3*j2*k2^2*l2*l4 +
            6*i1*i2*j2*k2^2*l2*l4 + 5*i1*j2^2*k2^2*l2*l4 + 4*i1^3*k2^3*l2*l4 + 2*i1*i2*k2^3*l2*l4 + 6*i1*j2*k2^3*l2*l4
            + 4*i1*k2^4*l2*l4 + 5*i1^5*i2*l2^2*l4 + 3*i1^3*i2^2*l2^2*l4 + i1*i2^3*l2^2*l4 + 4*i1^5*j2*l2^2*l4 +
            2*i1^3*i2*j2*l2^2*l4 + 5*i1*i2^2*j2*l2^2*l4 + 5*i1^3*j2^2*l2^2*l4 + i1*j2^3*l2^2*l4 + 6*i1^5*k2*l2^2*l4 +
            4*i1^3*i2*k2*l2^2*l4 + 4*i1*i2^2*k2*l2^2*l4 + 5*i1^3*j2*k2*l2^2*l4 + 6*i1*i2*j2*k2*l2^2*l4 +
            6*i1*j2^2*k2*l2^2*l4 + i1*i2*k2^2*l2^2*l4 + 4*i1*j2*k2^2*l2^2*l4 + 2*i1*k2^3*l2^2*l4 + 6*i1^3*i2*l2^3*l4 +
            2*i1*i2^2*l2^3*l4 + 2*i1^3*j2*l2^3*l4 + i1*i2*j2*l2^3*l4 + 3*i1^3*k2*l2^3*l4 + 4*i1*i2*k2*l2^3*l4 +
            4*i1*j2*k2*l2^3*l4 + 2*i1*k2^2*l2^3*l4 + 3*i1^3*l2^4*l4 + 5*i1*i2*l2^4*l4 + 5*i1*j2*l2^4*l4 +
            5*i1*k2*l2^4*l4 + 4*i1*l2^5*l4 + i1^6*i2*j3*l4 + 2*i1^4*i2^2*j3*l4 + 5*i2^4*j3*l4 + 4*i1^6*j2*j3*l4 +
            2*i1^4*i2*j2*j3*l4 + 2*i1^2*i2^2*j2*j3*l4 + i2^3*j2*j3*l4 + 4*i1^4*j2^2*j3*l4 + 3*i1^2*i2*j2^2*j3*l4 +
            2*i2^2*j2^2*j3*l4 + 2*i1^2*j2^3*j3*l4 + i2*j2^3*j3*l4 + 5*j2^4*j3*l4 + 3*i1^6*k2*j3*l4 +
            4*i1^2*i2^2*k2*j3*l4 + 3*i2^3*k2*j3*l4 + 4*i1^4*j2*k2*j3*l4 + 4*i1^2*i2*j2*k2*j3*l4 + 3*i2^2*j2*k2*j3*l4 +
            6*i2*j2^2*k2*j3*l4 + 2*j2^3*k2*j3*l4 + 6*i1^2*i2*k2^2*j3*l4 + 6*i2^2*k2^2*j3*l4 + 2*i1^2*j2*k2^2*j3*l4 +
            6*i2*j2*k2^2*j3*l4 + 3*j2^2*k2^2*j3*l4 + i1^2*k2^3*j3*l4 + 6*i2*k2^3*j3*l4 + 2*j2*k2^3*j3*l4 +
            5*k2^4*j3*l4 + 3*i1^6*l2*j3*l4 + i1^4*i2*l2*j3*l4 + 5*i1^2*i2^2*l2*j3*l4 + 2*i2^3*l2*j3*l4 +
            5*i1^4*j2*l2*j3*l4 + 2*i1^2*i2*j2*l2*j3*l4 + 6*i2^2*j2*l2*j3*l4 + 6*i1^2*j2^2*l2*j3*l4 +
            3*i2*j2^2*l2*j3*l4 + 3*j2^3*l2*j3*l4 + 6*i1^2*i2*k2*l2*j3*l4 + 3*i2^2*k2*l2*j3*l4 + 3*i1^2*j2*k2*l2*j3*l4
            + 2*i2*j2*k2*l2*j3*l4 + 2*j2^2*k2*l2*j3*l4 + i1^2*k2^2*l2*j3*l4 + 4*i2*k2^2*l2*j3*l4 + 3*j2*k2^2*l2*j3*l4
            + 3*k2^3*l2*j3*l4 + 3*i1^4*l2^2*j3*l4 + 5*i1^2*i2*l2^2*j3*l4 + 5*i2^2*l2^2*j3*l4 + 6*i1^2*j2*l2^2*j3*l4 +
            3*i2*j2*l2^2*j3*l4 + 2*j2^2*l2^2*j3*l4 + 3*i1^2*k2*l2^2*j3*l4 + j2*k2*l2^2*j3*l4 + 6*k2^2*l2^2*j3*l4 +
            6*i1^2*l2^3*j3*l4 + 2*j2*l2^3*j3*l4 + 4*k2*l2^3*j3*l4 + 6*i1^5*j3^2*l4 + i1^3*i2*j3^2*l4 +
            5*i1*i2^2*j3^2*l4 + 5*i1^3*j2*j3^2*l4 + 3*i1*i2*j2*j3^2*l4 + 6*i1*j2^2*j3^2*l4 + 3*i1^3*k2*j3^2*l4 +
            2*i1*i2*k2*j3^2*l4 + i1*j2*k2*j3^2*l4 + 3*i1^3*l2*j3^2*l4 + i1*i2*l2*j3^2*l4 + i1*j2*l2*j3^2*l4 +
            2*i1*k2*l2*j3^2*l4 + 5*i1*l2^2*j3^2*l4 + 2*i1^2*j3^3*l4 + 5*i2*j3^3*l4 + 6*j2*j3^3*l4 + 5*i1^8*k3*l4 +
            5*i1^6*i2*k3*l4 + 4*i1^4*i2^2*k3*l4 + 3*i1^2*i2^3*k3*l4 + 3*i1^6*j2*k3*l4 + 4*i1^4*i2*j2*k3*l4 +
            3*i1^2*i2^2*j2*k3*l4 + 6*i1^2*i2*j2^2*k3*l4 + 2*i1^2*j2^3*k3*l4 + 5*i1^4*i2*k2*k3*l4 +
            4*i1^2*i2^2*k2*k3*l4 + 5*i2^3*k2*k3*l4 + 2*i1^4*j2*k2*k3*l4 + 5*i2^2*j2*k2*k3*l4 + 6*i1^2*j2^2*k2*k3*l4 +
            3*i2*j2^2*k2*k3*l4 + j2^3*k2*k3*l4 + 3*i1^4*k2^2*k3*l4 + 2*i1^2*i2*k2^2*k3*l4 + 5*i2^2*k2^2*k3*l4 +
            i1^2*j2*k2^2*k3*l4 + 6*i2*j2*k2^2*k3*l4 + j2^2*k2^2*k3*l4 + i1^2*k2^3*k3*l4 + 2*i2*k2^3*k3*l4 +
            5*j2*k2^3*k3*l4 + 5*k2^4*k3*l4 + 4*i1^4*i2*l2*k3*l4 + i1^2*i2^2*l2*k3*l4 + i2^3*l2*k3*l4 +
            i1^2*i2*j2*l2*k3*l4 + 6*i1^2*j2^2*l2*k3*l4 + 4*i2*j2^2*l2*k3*l4 + 2*j2^3*l2*k3*l4 + 6*i1^4*k2*l2*k3*l4 +
            6*i2^2*k2*l2*k3*l4 + i2*j2*k2*l2*k3*l4 + 3*i1^2*k2^2*l2*k3*l4 + 4*i2*k2^2*l2*k3*l4 + 2*j2*k2^2*l2*k3*l4 +
            5*k2^3*l2*k3*l4 + 4*i1^2*i2*l2^2*k3*l4 + 4*i2^2*l2^2*k3*l4 + i1^2*j2*l2^2*k3*l4 + i2*j2*l2^2*k3*l4 +
            5*j2^2*l2^2*k3*l4 + 6*i1^2*k2*l2^2*k3*l4 + 3*i2*k2*l2^2*k3*l4 + k2^2*l2^2*k3*l4 + 3*i1^2*l2^3*k3*l4 +
            5*i2*l2^3*k3*l4 + 4*j2*l2^3*k3*l4 + 2*l2^4*k3*l4 + 3*i1^5*j3*k3*l4 + i1^3*i2*j3*k3*l4 + 4*i1^3*j2*j3*k3*l4
            + 2*i1*i2*j2*j3*k3*l4 + 5*i1*j2^2*j3*k3*l4 + 3*i1^3*k2*j3*k3*l4 + 4*i1*i2*k2*j3*k3*l4 +
            6*i1*j2*k2*j3*k3*l4 + i1*k2^2*j3*k3*l4 + i1^3*l2*j3*k3*l4 + 5*i1*i2*l2*j3*k3*l4 + 4*i1*j2*l2*j3*k3*l4 +
            2*i1*k2*l2*j3*k3*l4 + 6*i1^2*j3^2*k3*l4 + 6*i2*j3^2*k3*l4 + 5*j2*j3^2*k3*l4 + 6*k2*j3^2*k3*l4 +
            5*l2*j3^2*k3*l4 + 2*i1^3*i2*k3^2*l4 + 4*i1*i2^2*k3^2*l4 + 2*i1^3*j2*k3^2*l4 + 3*i1*j2^2*k3^2*l4 +
            2*i1^3*k2*k3^2*l4 + 4*i1*i2*k2*k3^2*l4 + 5*i1*j2*k2*k3^2*l4 + 4*i1*k2^2*k3^2*l4 + 3*i1^3*l2*k3^2*l4 +
            i1*i2*l2*k3^2*l4 + i1*j2*l2*k3^2*l4 + 6*i1*k2*l2*k3^2*l4 + 5*i1*l2^2*k3^2*l4 + 3*i1^2*j3*k3^2*l4 +
            4*i2*j3*k3^2*l4 + 4*j2*j3*k3^2*l4 + 2*k2*j3*k3^2*l4 + 2*l2*j3*k3^2*l4 + i1^2*k3^3*l4 + 6*i2*k3^3*l4 +
            4*j2*k3^3*l4 + k2*k3^3*l4 + 6*l2*k3^3*l4 + 6*i1^8*l3*l4 + 2*i1^6*i2*l3*l4 + 2*i1^4*i2^2*l3*l4 +
            6*i1^2*i2^3*l3*l4 + 2*i1^6*j2*l3*l4 + 4*i1^4*i2*j2*l3*l4 + i1^2*i2^2*j2*l3*l4 + 6*i2^3*j2*l3*l4 +
            2*i1^2*i2*j2^2*l3*l4 + 3*i2^2*j2^2*l3*l4 + 5*i1^2*j2^3*l3*l4 + 4*i2*j2^3*l3*l4 + j2^4*l3*l4 +
            2*i1^6*k2*l3*l4 + 6*i1^4*i2*k2*l3*l4 + 5*i1^2*i2^2*k2*l3*l4 + 3*i2^3*k2*l3*l4 + 4*i1^2*i2*j2*k2*l3*l4 +
            2*i2^2*j2*k2*l3*l4 + 4*i1^2*j2^2*k2*l3*l4 + 2*i2*j2^2*k2*l3*l4 + 2*i1^4*k2^2*l3*l4 + i1^2*i2*k2^2*l3*l4 +
            4*i2^2*k2^2*l3*l4 + 5*i1^2*j2*k2^2*l3*l4 + 5*i2*j2*k2^2*l3*l4 + 2*j2^2*k2^2*l3*l4 + 5*i1^2*k2^3*l3*l4 +
            4*i2*k2^3*l3*l4 + 2*j2*k2^3*l3*l4 + 2*k2^4*l3*l4 + 6*i1^6*l2*l3*l4 + i1^4*i2*l2*l3*l4 + i2^3*l2*l3*l4 +
            3*i1^4*j2*l2*l3*l4 + i1^2*i2*j2*l2*l3*l4 + 2*i2^2*j2*l2*l3*l4 + 2*i1^2*j2^2*l2*l3*l4 + 3*i2*j2^2*l2*l3*l4
            + j2^3*l2*l3*l4 + 5*i1^4*k2*l2*l3*l4 + 5*i1^2*i2*k2*l2*l3*l4 + 3*i2^2*k2*l2*l3*l4 + 4*i1^2*j2*k2*l2*l3*l4
            + 3*i2*j2*k2*l2*l3*l4 + 2*j2^2*k2*l2*l3*l4 + 2*i1^2*k2^2*l2*l3*l4 + 6*i2*k2^2*l2*l3*l4 +
            4*j2*k2^2*l2*l3*l4 + 6*k2^3*l2*l3*l4 + 4*i1^4*l2^2*l3*l4 + 4*i1^2*i2*l2^2*l3*l4 + 6*i2^2*l2^2*l3*l4 +
            6*i1^2*j2*l2^2*l3*l4 + i2*j2*l2^2*l3*l4 + 3*j2^2*l2^2*l3*l4 + 3*i1^2*k2*l2^2*l3*l4 + 4*i2*k2*l2^2*l3*l4 +
            4*j2*k2*l2^2*l3*l4 + 2*k2^2*l2^2*l3*l4 + 3*i1^2*l2^3*l3*l4 + 2*i2*l2^3*l3*l4 + 3*j2*l2^3*l3*l4 +
            6*k2*l2^3*l3*l4 + 4*l2^4*l3*l4 + 2*i1^5*j3*l3*l4 + 3*i1^3*i2*j3*l3*l4 + 2*i1*i2^2*j3*l3*l4 +
            6*i1*i2*j2*j3*l3*l4 + 6*i1*j2^2*j3*l3*l4 + 3*i1*j2*k2*j3*l3*l4 + 3*i1*k2^2*j3*l3*l4 + i1^3*l2*j3*l3*l4 +
            4*i1*i2*l2*j3*l3*l4 + 2*i1*k2*l2*j3*l3*l4 + i1*l2^2*j3*l3*l4 + 6*i1^2*j3^2*l3*l4 + 4*i2*j3^2*l3*l4 +
            4*j2*j3^2*l3*l4 + 2*i1^5*l3^2*l4 + 3*i1^3*i2*l3^2*l4 + i1*i2^2*l3^2*l4 + 2*i1^3*j2*l3^2*l4 +
            2*i1*i2*j2*l3^2*l4 + 3*i1*j2^2*l3^2*l4 + 5*i1^3*k2*l3^2*l4 + 2*i1*i2*k2*l3^2*l4 + 3*i1*j2*k2*l3^2*l4 +
            5*i1^3*l2*l3^2*l4 + i1*i2*l2*l3^2*l4 + 3*i1*j2*l2*l3^2*l4 + 6*i1*k2*l2*l3^2*l4 + i1*l2^2*l3^2*l4 +
            i1^2*j3*l3^2*l4 + 6*i2*j3*l3^2*l4 + 5*j2*j3*l3^2*l4 + 3*k2*j3*l3^2*l4 + 2*l2*j3*l3^2*l4 + 4*i1^2*l3^3*l4 +
            2*i2*l3^3*l4 + 5*j2*l3^3*l4 + 2*i1^7*l4^2 + 4*i1^5*i2*l4^2 + 6*i1^3*i2^2*l4^2 + 4*i1^3*i2*j2*l4^2 +
            6*i1*i2^2*j2*l4^2 + 5*i1^3*j2^2*l4^2 + 2*i1*i2*j2^2*l4^2 + 6*i1*j2^3*l4^2 + i1^5*k2*l4^2 +
            4*i1^3*i2*k2*l4^2 + i1*i2^2*k2*l4^2 + 3*i1^3*j2*k2*l4^2 + 5*i1*j2^2*k2*l4^2 + 5*i1^3*k2^2*l4^2 +
            4*i1*k2^3*l4^2 + 4*i1^3*i2*l2*l4^2 + 2*i1*i2^2*l2*l4^2 + i1*i2*j2*l2*l4^2 + 4*i1*j2^2*l2*l4^2 +
            2*i1^3*k2*l2*l4^2 + 5*i1*i2*k2*l2*l4^2 + 3*i1*j2*k2*l2*l4^2 + i1^3*l2^2*l4^2 + 6*i1*i2*l2^2*l4^2 +
            i1*j2*l2^2*l4^2 + i1*k2*l2^2*l4^2 + 6*i1*l2^3*l4^2 + 2*i1^4*j3*l4^2 + 3*i1^2*i2*j3*l4^2 + 5*i2^2*j3*l4^2 +
            3*i1^2*j2*j3*l4^2 + 4*i2*j2*j3*l4^2 + 5*j2^2*j3*l4^2 + 2*i1^2*k2*j3*l4^2 + 6*i2*k2*j3*l4^2 +
            5*j2*k2*j3*l4^2 + 5*k2^2*j3*l4^2 + 4*i1^2*l2*j3*l4^2 + 3*j2*l2*j3*l4^2 + 5*k2*l2*j3*l4^2 + l2^2*j3*l4^2 +
            6*i1^4*k3*l4^2 + 2*i1^2*i2*k3*l4^2 + 5*i2^2*k3*l4^2 + 6*i1^2*j2*k3*l4^2 + 5*i2*j2*k3*l4^2 + 4*j2^2*k3*l4^2
            + 6*i1^2*k2*k3*l4^2 + i2*k2*k3*l4^2 + 4*j2*k2*k3*l4^2 + 4*i1^2*l2*k3*l4^2 + 6*i2*l2*k3*l4^2 +
            j2*l2*k3*l4^2 + 4*k2*l2*k3*l4^2 + 3*l2^2*k3*l4^2 + i1*j3*k3*l4^2 + i1^4*l3*l4^2 + 4*i1^2*i2*l3*l4^2 +
            3*i2^2*l3*l4^2 + 5*i1^2*j2*l3*l4^2 + 6*i2*j2*l3*l4^2 + 5*j2^2*l3*l4^2 + 6*i1^2*k2*l3*l4^2 + i2*k2*l3*l4^2
            + 2*j2*k2*l3*l4^2 + 4*k2^2*l3*l4^2 + 4*i1^2*l2*l3*l4^2 + 2*i2*l2*l3*l4^2 + 6*k2*l2*l3*l4^2 +
            2*l2^2*l3*l4^2 + 3*i1*j3*l3*l4^2 + 3*i1^3*l4^3 + 3*i1*i2*l4^3 + 4*i1*l2*l4^3 + j3*l4^3 + 2*l3*l4^3 +
            i1^11*m4 + i1^5*i2^3*m4 + 5*i1^3*i2^4*m4 + 6*i1*i2^5*m4 + 4*i1^9*j2*m4 + i1^7*i2*j2*m4 + 4*i1^5*i2^2*j2*m4
            + i1^3*i2^3*j2*m4 + 5*i1*i2^4*j2*m4 + 2*i1^5*i2*j2^2*m4 + 2*i1^3*i2^2*j2^2*m4 + 4*i1*i2^3*j2^2*m4 +
            i1^3*i2*j2^3*m4 + 3*i1*i2^2*j2^3*m4 + 5*i1^3*j2^4*m4 + 2*i1*i2*j2^4*m4 + i1*j2^5*m4 + 6*i1^7*i2*k2*m4 +
            3*i1^5*i2^2*k2*m4 + 4*i1^3*i2^3*k2*m4 + 3*i1*i2^4*k2*m4 + i1^7*j2*k2*m4 + i1^3*i2^2*j2*k2*m4 +
            4*i1*i2^3*j2*k2*m4 + 4*i1^5*j2^2*k2*m4 + 5*i1*i2^2*j2^2*k2*m4 + 2*i1^3*j2^3*k2*m4 + i1*i2*j2^3*k2*m4 +
            i1*j2^4*k2*m4 + 4*i1^7*k2^2*m4 + 6*i1^5*i2*k2^2*m4 + i1^3*i2^2*k2^2*m4 + 5*i1*i2^3*k2^2*m4 +
            6*i1^5*j2*k2^2*m4 + 4*i1^3*i2*j2*k2^2*m4 + i1*i2^2*j2*k2^2*m4 + 2*i1*i2*j2^2*k2^2*m4 + 6*i1*j2^3*k2^2*m4 +
            5*i1^3*i2*k2^3*m4 + 3*i1*i2^2*k2^3*m4 + 3*i1^3*j2*k2^3*m4 + 5*i1*i2*j2*k2^3*m4 + 2*i1*j2^2*k2^3*m4 +
            5*i1^3*k2^4*m4 + i1*i2*k2^4*m4 + i1*j2*k2^4*m4 + 5*i1*k2^5*m4 + 5*i1^5*i2^2*l2*m4 + 4*i1^3*i2^3*l2*m4 +
            2*i1*i2^4*l2*m4 + 6*i1^7*j2*l2*m4 + i1^5*i2*j2*l2*m4 + 4*i1^3*i2^2*j2*l2*m4 + 3*i1*i2^3*j2*l2*m4 +
            5*i1^3*i2*j2^2*l2*m4 + i1^3*j2^3*l2*m4 + 4*i1*i2*j2^3*l2*m4 + 5*i1*j2^4*l2*m4 + 5*i1^7*k2*l2*m4 +
            2*i1^5*i2*k2*l2*m4 + 4*i1^3*i2^2*k2*l2*m4 + i1*i2^3*k2*l2*m4 + 3*i1*i2^2*j2*k2*l2*m4 +
            6*i1^3*j2^2*k2*l2*m4 + 4*i1*i2*j2^2*k2*l2*m4 + 6*i1*j2^3*k2*l2*m4 + i1^5*k2^2*l2*m4 + 4*i1^3*i2*k2^2*l2*m4
            + 6*i1*i2^2*k2^2*l2*m4 + 5*i1^3*j2*k2^2*l2*m4 + 3*i1*i2*j2*k2^2*l2*m4 + 5*i1*j2^2*k2^2*l2*m4 +
            3*i1^3*k2^3*l2*m4 + 5*i1*i2*k2^3*l2*m4 + 6*i1*j2*k2^3*l2*m4 + 6*i1*k2^4*l2*m4 + 5*i1^7*l2^2*m4 +
            3*i1^5*i2*l2^2*m4 + 4*i1^3*i2^2*l2^2*m4 + 3*i1^5*j2*l2^2*m4 + 4*i1^3*i2*j2*l2^2*m4 + 2*i1*i2^2*j2*l2^2*m4
            + 5*i1^3*j2^2*l2^2*m4 + 3*i1*i2*j2^2*l2^2*m4 + 2*i1*j2^3*l2^2*m4 + 6*i1^3*i2*k2*l2^2*m4 +
            6*i1*i2^2*k2*l2^2*m4 + 4*i1^3*j2*k2*l2^2*m4 + 4*i1*i2*j2*k2*l2^2*m4 + 6*i1*j2^2*k2*l2^2*m4 +
            5*i1^3*k2^2*l2^2*m4 + 4*i1*i2*k2^2*l2^2*m4 + 5*i1*j2*k2^2*l2^2*m4 + 5*i1*k2^3*l2^2*m4 + 6*i1^5*l2^3*m4 +
            i1^3*i2*l2^3*m4 + i1^3*j2*l2^3*m4 + i1^3*k2*l2^3*m4 + 2*i1*i2*k2*l2^3*m4 + 6*i1*j2*k2*l2^3*m4 +
            i1*k2^2*l2^3*m4 + 6*i1^3*l2^4*m4 + 6*i1*j2*l2^4*m4 + 5*i1*l2^5*m4 + 2*i1^8*i3*m4 + 5*i1^4*i2^2*i3*m4 +
            4*i1^2*i2^3*i3*m4 + i1^6*j2*i3*m4 + 2*i1^4*i2*j2*i3*m4 + 2*i1^2*i2^2*j2*i3*m4 + 5*i1^2*i2*j2^2*i3*m4 +
            3*i1^2*j2^3*i3*m4 + i1^6*k2*i3*m4 + 2*i1^4*i2*k2*i3*m4 + 6*i1^2*i2^2*k2*i3*m4 + 5*i2^3*k2*i3*m4 +
            2*i1^4*j2*k2*i3*m4 + i1^2*i2*j2*k2*i3*m4 + 6*i2^2*j2*k2*i3*m4 + i2*j2^2*k2*i3*m4 + 2*j2^3*k2*i3*m4 +
            3*i1^4*k2^2*i3*m4 + 4*i1^2*i2*k2^2*i3*m4 + 6*i2^2*k2^2*i3*m4 + i1^2*j2*k2^2*i3*m4 + i2*j2*k2^2*i3*m4 +
            4*i1^2*k2^3*i3*m4 + i2*k2^3*i3*m4 + 3*j2*k2^3*i3*m4 + 2*k2^4*i3*m4 + 3*i1^6*l2*i3*m4 + i1^4*i2*l2*i3*m4 +
            3*i1^2*i2^2*l2*i3*m4 + 3*i2^3*l2*i3*m4 + 6*i1^4*j2*l2*i3*m4 + 5*i1^2*i2*j2*l2*i3*m4 + 5*i2^2*j2*l2*i3*m4 +
            6*i1^2*j2^2*l2*i3*m4 + 2*i2*j2^2*l2*i3*m4 + 4*j2^3*l2*i3*m4 + 6*i1^4*k2*l2*i3*m4 + 4*i1^2*i2*k2*l2*i3*m4 +
            i1^2*j2*k2*l2*i3*m4 + 4*i2*j2*k2*l2*i3*m4 + 3*j2^2*k2*l2*i3*m4 + 3*i1^2*k2^2*l2*i3*m4 + 5*i2*k2^2*l2*i3*m4
            + 4*j2*k2^2*l2*i3*m4 + 2*k2^3*l2*i3*m4 + 3*i1^4*l2^2*i3*m4 + 5*i1^2*i2*l2^2*i3*m4 + i2^2*l2^2*i3*m4 +
            5*i2*j2*l2^2*i3*m4 + j2^2*l2^2*i3*m4 + 3*i2*k2*l2^2*i3*m4 + 3*j2*k2*l2^2*i3*m4 + 6*k2^2*l2^2*i3*m4 +
            2*i1^2*l2^3*i3*m4 + i2*l2^3*i3*m4 + 3*j2*l2^3*i3*m4 + 2*k2*l2^3*i3*m4 + 2*l2^4*i3*m4 + 2*i1^5*i3^2*m4 +
            4*i1^3*i2*i3^2*m4 + 3*i1^3*j2*i3^2*m4 + 5*i1^3*k2*i3^2*m4 + 5*i1*i2*k2*i3^2*m4 + 2*i1*j2*k2*i3^2*m4 +
            2*i1*k2^2*i3^2*m4 + 6*i1^3*l2*i3^2*m4 + 2*i1*i2*l2*i3^2*m4 + 5*i1*j2*l2*i3^2*m4 + 4*i1*k2*l2*i3^2*m4 +
            6*i1*l2^2*i3^2*m4 + 3*i1^8*j3*m4 + 6*i1^6*i2*j3*m4 + 6*i1^4*i2^2*j3*m4 + 2*i1^2*i2^3*j3*m4 + 6*i2^4*j3*m4
            + 3*i1^6*j2*j3*m4 + 6*i1^4*i2*j2*j3*m4 + 6*i1^2*i2^2*j2*j3*m4 + 4*i2^3*j2*j3*m4 + 5*i1^4*j2^2*j3*m4 +
            3*i1^2*i2*j2^2*j3*m4 + i2^2*j2^2*j3*m4 + 3*i1^2*j2^3*j3*m4 + 4*i2*j2^3*j3*m4 + 6*j2^4*j3*m4 +
            5*i1^4*i2*k2*j3*m4 + 4*i1^2*i2^2*k2*j3*m4 + 2*i2^3*k2*j3*m4 + 3*i1^4*j2*k2*j3*m4 + 5*i1^2*i2*j2*k2*j3*m4 +
            3*i2^2*j2*k2*j3*m4 + 2*i2*j2^2*k2*j3*m4 + 2*i1^4*k2^2*j3*m4 + 5*i1^2*i2*k2^2*j3*m4 + 3*i2^2*k2^2*j3*m4 +
            3*i1^2*j2*k2^2*j3*m4 + 3*i2*j2*k2^2*j3*m4 + 6*j2^2*k2^2*j3*m4 + i1^2*k2^3*j3*m4 + 6*i2*k2^3*j3*m4 +
            5*j2*k2^3*j3*m4 + i1^6*l2*j3*m4 + 4*i1^2*i2^2*l2*j3*m4 + 4*i2^3*l2*j3*m4 + 4*i1^4*j2*l2*j3*m4 +
            4*i2^2*j2*l2*j3*m4 + 3*i1^2*j2^2*l2*j3*m4 + i2*j2^2*l2*j3*m4 + 5*j2^3*l2*j3*m4 + 6*i1^4*k2*l2*j3*m4 +
            6*i1^2*i2*k2*l2*j3*m4 + 2*i2^2*k2*l2*j3*m4 + 2*i1^2*j2*k2*l2*j3*m4 + 2*j2^2*k2*l2*j3*m4 +
            6*i1^2*k2^2*l2*j3*m4 + 2*i2*k2^2*l2*j3*m4 + 2*j2*k2^2*l2*j3*m4 + 6*k2^3*l2*j3*m4 + 5*i1^4*l2^2*j3*m4 +
            5*i1^2*i2*l2^2*j3*m4 + 6*i2^2*l2^2*j3*m4 + 3*i1^2*j2*l2^2*j3*m4 + 3*i2*j2*l2^2*j3*m4 + 4*j2^2*l2^2*j3*m4 +
            5*i1^2*k2*l2^2*j3*m4 + 5*i2*k2*l2^2*j3*m4 + 3*j2*k2*l2^2*j3*m4 + 6*k2^2*l2^2*j3*m4 + 2*i1^2*l2^3*j3*m4 +
            i2*l2^3*j3*m4 + 2*k2*l2^3*j3*m4 + l2^4*j3*m4 + 6*i1^5*j3^2*m4 + i1^3*j2*j3^2*m4 + 2*i1*i2*j2*j3^2*m4 +
            5*i1*j2^2*j3^2*m4 + 4*i1^3*k2*j3^2*m4 + 6*i1*j2*k2*j3^2*m4 + 6*i1*k2^2*j3^2*m4 + 3*i1^3*l2*j3^2*m4 +
            6*i1*i2*l2*j3^2*m4 + 2*i1*j2*l2*j3^2*m4 + 5*i1*k2*l2*j3^2*m4 + 5*i1*l2^2*j3^2*m4 + i1^2*j3^3*m4 +
            i2*j3^3*m4 + 5*j2*j3^3*m4 + 6*k2*j3^3*m4 + 4*l2*j3^3*m4 + i1^8*l3*m4 + 3*i1^6*i2*l3*m4 + 4*i1^4*i2^2*l3*m4
            + 4*i1^2*i2^3*l3*m4 + i2^4*l3*m4 + 3*i1^6*j2*l3*m4 + 3*i1^4*i2*j2*l3*m4 + i1^2*i2^2*j2*l3*m4 +
            5*i2^3*j2*l3*m4 + 3*i1^4*j2^2*l3*m4 + 3*i1^2*i2*j2^2*l3*m4 + 6*i1^2*j2^3*l3*m4 + 2*i2*j2^3*l3*m4 +
            6*j2^4*l3*m4 + 4*i1^6*k2*l3*m4 + 2*i1^4*i2*k2*l3*m4 + 5*i1^2*i2^2*k2*l3*m4 + 4*i2^3*k2*l3*m4 +
            5*i1^4*j2*k2*l3*m4 + 4*i1^2*i2*j2*k2*l3*m4 + 4*i1^2*j2^2*k2*l3*m4 + 2*i2*j2^2*k2*l3*m4 + j2^3*k2*l3*m4 +
            5*i1^4*k2^2*l3*m4 + 5*i1^2*i2*k2^2*l3*m4 + 3*i2^2*k2^2*l3*m4 + 4*i1^2*j2*k2^2*l3*m4 + 3*j2^2*k2^2*l3*m4 +
            5*i1^2*k2^3*l3*m4 + 6*i2*k2^3*l3*m4 + j2*k2^3*l3*m4 + 3*k2^4*l3*m4 + 2*i1^6*l2*l3*m4 + 4*i1^4*i2*l2*l3*m4
            + i1^2*i2^2*l2*l3*m4 + 4*i1^4*j2*l2*l3*m4 + 2*i2^2*j2*l2*l3*m4 + 6*i1^2*j2^2*l2*l3*m4 + 5*j2^3*l2*l3*m4 +
            i1^4*k2*l2*l3*m4 + 5*i2^2*k2*l2*l3*m4 + 4*i1^2*j2*k2*l2*l3*m4 + 2*i2*j2*k2*l2*l3*m4 + 3*j2^2*k2*l2*l3*m4 +
            3*i1^2*k2^2*l2*l3*m4 + 5*i2*k2^2*l2*l3*m4 + 2*j2*k2^2*l2*l3*m4 + k2^3*l2*l3*m4 + 5*i1^4*l2^2*l3*m4 +
            3*i1^2*i2*l2^2*l3*m4 + 4*i1^2*j2*l2^2*l3*m4 + 5*i2*j2*l2^2*l3*m4 + 4*j2^2*l2^2*l3*m4 + i1^2*k2*l2^2*l3*m4
            + 3*i2*k2*l2^2*l3*m4 + 4*j2*k2*l2^2*l3*m4 + k2^2*l2^2*l3*m4 + 2*i1^2*l2^3*l3*m4 + i2*l2^3*l3*m4 +
            6*j2*l2^3*l3*m4 + 2*k2*l2^3*l3*m4 + 5*l2^4*l3*m4 + 4*i1^5*i3*l3*m4 + 6*i1^3*i2*i3*l3*m4 +
            5*i1*i2^2*i3*l3*m4 + 6*i1^3*j2*i3*l3*m4 + 6*i1*i2*j2*i3*l3*m4 + 3*i1*j2^2*i3*l3*m4 + 4*i1*i2*k2*i3*l3*m4 +
            4*i1*j2*k2*i3*l3*m4 + i1*k2^2*i3*l3*m4 + 6*i1^3*l2*i3*l3*m4 + i1*i2*l2*i3*l3*m4 + 4*i1*j2*l2*i3*l3*m4 +
            4*i1*k2*l2*i3*l3*m4 + 5*i1*l2^2*i3*l3*m4 + i2*i3^2*l3*m4 + 6*j2*i3^2*l3*m4 + k2*i3^2*l3*m4 +
            3*l2*i3^2*l3*m4 + 5*i1^5*j3*l3*m4 + 2*i1^3*i2*j3*l3*m4 + i1*i2^2*j3*l3*m4 + i1^3*j2*j3*l3*m4 +
            2*i1*j2^2*j3*l3*m4 + i1^3*k2*j3*l3*m4 + 3*i1*j2*k2*j3*l3*m4 + i1*k2^2*j3*l3*m4 + i1*i2*l2*j3*l3*m4 +
            3*i1*j2*l2*j3*l3*m4 + i1*k2*l2*j3*l3*m4 + 5*i1*l2^2*j3*l3*m4 + 3*i1^2*j3^2*l3*m4 + 4*i2*j3^2*l3*m4 +
            j2*j3^2*l3*m4 + 6*k2*j3^2*l3*m4 + 4*l2*j3^2*l3*m4 + 5*i1^5*l3^2*m4 + 4*i1^3*i2*l3^2*m4 + 5*i1^3*j2*l3^2*m4
            + i1*j2^2*l3^2*m4 + 3*i1*i2*k2*l3^2*m4 + 6*i1*k2^2*l3^2*m4 + 2*i1^3*l2*l3^2*m4 + i1*j2*l2*l3^2*m4 +
            i1*k2*l2*l3^2*m4 + i1*l2^2*l3^2*m4 + 6*i1^2*i3*l3^2*m4 + i2*i3*l3^2*m4 + 6*j2*i3*l3^2*m4 + 4*k2*i3*l3^2*m4
            + 2*l2*i3*l3^2*m4 + 3*i1^2*j3*l3^2*m4 + 3*i2*j3*l3^2*m4 + 2*j2*j3*l3^2*m4 + k2*j3*l3^2*m4 +
            3*l2*j3*l3^2*m4 + 4*i1^2*l3^3*m4 + 5*i2*l3^3*m4 + 6*j2*l3^3*m4 + 2*k2*l3^3*m4 + 6*l2*l3^3*m4 + i1^10*i5 +
            2*i1^8*i2*i5 + 2*i1^6*i2^2*i5 + 3*i1^8*j2*i5 + 3*i1^6*i2*j2*i5 + 2*i1^6*j2^2*i5 + 3*i1^8*k2*i5 +
            2*i1^6*i2*k2*i5 + 4*i1^4*i2^2*k2*i5 + 2*i1^6*j2*k2*i5 + 6*i1^4*i2*j2*k2*i5 + 4*i1^4*j2^2*k2*i5 +
            2*i1^6*k2^2*i5 + 6*i1^4*i2*k2^2*i5 + i1^2*i2^2*k2^2*i5 + 5*i1^4*j2*k2^2*i5 + 5*i1^2*i2*j2*k2^2*i5 +
            i1^2*j2^2*k2^2*i5 + 5*i1^4*k2^3*i5 + 2*i1^2*i2*k2^3*i5 + 3*i1^2*j2*k2^3*i5 + i1^2*k2^4*i5 + 5*i2*k2^4*i5 +
            6*j2*k2^4*i5 + 5*k2^5*i5 + 6*i1^8*l2*i5 + 6*i1^6*i2*l2*i5 + 4*i1^4*i2^2*l2*i5 + 6*i1^6*j2*l2*i5 +
            6*i1^4*i2*j2*l2*i5 + 4*i1^4*j2^2*l2*i5 + 3*i1^6*k2*l2*i5 + 5*i1^2*i2^2*k2*l2*i5 + 5*i1^4*j2*k2*l2*i5 +
            4*i1^2*i2*j2*k2*l2*i5 + 5*i1^2*j2^2*k2*l2*i5 + 2*i1^4*k2^2*l2*i5 + 3*i1^2*i2*k2^2*l2*i5 +
            5*i2^2*k2^2*l2*i5 + i1^2*j2*k2^2*l2*i5 + 4*i2*j2*k2^2*l2*i5 + 5*j2^2*k2^2*l2*i5 + 2*i1^2*k2^3*l2*i5 +
            i2*k2^3*l2*i5 + 3*j2*k2^3*l2*i5 + k2^4*l2*i5 + 2*i1^6*l2^2*i5 + 5*i1^4*i2*l2^2*i5 + 6*i1^2*i2^2*l2^2*i5 +
            2*i1^2*i2*j2*l2^2*i5 + 6*i1^2*j2^2*l2^2*i5 + 3*i1^2*i2*k2*l2^2*i5 + 4*i2^2*k2*l2^2*i5 +
            2*i1^2*j2*k2*l2^2*i5 + 6*i2*j2*k2*l2^2*i5 + 4*j2^2*k2*l2^2*i5 + 5*i1^2*k2^2*l2^2*i5 + 6*i2*k2^2*l2^2*i5 +
            2*j2*k2^2*l2^2*i5 + 5*k2^3*l2^2*i5 + 2*i1^2*i2*l2^3*i5 + 3*i2^2*l2^3*i5 + 3*i1^2*j2*l2^3*i5 +
            i2*j2*l2^3*i5 + 3*j2^2*l2^3*i5 + i1^2*k2*l2^3*i5 + 6*i2*k2*l2^3*i5 + 4*j2*k2*l2^3*i5 + i1^2*l2^4*i5 +
            j2*l2^4*i5 + 6*k2*l2^4*i5 + l2^5*i5 + 4*i1^7*i3*i5 + i1^5*k2*i3*i5 + 5*i1^3*k2^2*i3*i5 + 5*i1*k2^3*i3*i5 +
            6*i1^5*l2*i3*i5 + 5*i1^3*k2*l2*i3*i5 + 5*i1*k2^2*l2*i3*i5 + 3*i1^3*l2^2*i3*i5 + 5*i1*k2*l2^2*i3*i5 +
            i1*l2^3*i3*i5 + 6*i1^7*j3*i5 + 5*i1^5*k2*j3*i5 + i1^3*i2*k2*j3*i5 + 6*i1^3*j2*k2*j3*i5 + 5*i1^3*k2^2*j3*i5
            + 3*i1*i2*k2^2*j3*i5 + 4*i1*j2*k2^2*j3*i5 + 6*i1*k2^3*j3*i5 + 3*i1^5*l2*j3*i5 + 3*i1^3*i2*l2*j3*i5 +
            4*i1^3*j2*l2*j3*i5 + 4*i1^3*k2*l2*j3*i5 + 2*i1*i2*k2*l2*j3*i5 + 5*i1*j2*k2*l2*j3*i5 + i1*k2^2*l2*j3*i5 +
            i1*i2*l2^2*j3*i5 + 6*i1*j2*l2^2*j3*i5 + 4*i1*l2^3*j3*i5 + 3*i1^4*j3^2*i5 + 4*i1^2*k2*j3^2*i5 +
            6*k2^2*j3^2*i5 + 3*k2*l2*j3^2*i5 + 6*l2^2*j3^2*i5 + 2*i1^10*j5 + 6*i1^8*i2*j5 + 2*i1^6*i2^2*j5 +
            2*i1^4*i2^3*j5 + 6*i1^2*i2^4*j5 + 5*i2^5*j5 + i1^8*j2*j5 + 3*i1^6*i2*j2*j5 + 3*i1^4*i2^2*j2*j5 +
            i1^2*i2^3*j2*j5 + 3*i2^4*j2*j5 + i1^6*j2^2*j5 + 3*i1^4*i2*j2^2*j5 + 3*i1^2*i2^2*j2^2*j5 + i2^3*j2^2*j5 +
            6*i1^4*j2^3*j5 + 2*i1^2*i2*j2^3*j5 + 6*i2^2*j2^3*j5 + 2*i1^2*j2^4*j5 + 4*i2*j2^4*j5 + 2*j2^5*j5 +
            5*i1^8*k2*j5 + 6*i1^6*i2*k2*j5 + 3*i1^4*i2^2*k2*j5 + 5*i1^2*i2^3*k2*j5 + 2*i2^4*k2*j5 + 2*i1^6*j2*k2*j5 +
            5*i1^4*i2*j2*k2*j5 + 4*i1^2*i2^2*j2*k2*j5 + i2^3*j2*k2*j5 + 4*i1^4*j2^2*k2*j5 + 4*i1^2*i2*j2^2*k2*j5 +
            6*i2^2*j2^2*k2*j5 + i1^2*j2^3*k2*j5 + 5*i2*j2^3*k2*j5 + i1^6*k2^2*j5 + i1^4*i2*k2^2*j5 +
            3*i1^2*i2^2*k2^2*j5 + i1^4*j2*k2^2*j5 + 3*i1^2*i2*j2*k2^2*j5 + 6*i2^2*j2*k2^2*j5 + i1^2*j2^2*k2^2*j5 +
            6*i2*j2^2*k2^2*j5 + 2*j2^3*k2^2*j5 + 5*i1^4*k2^3*j5 + 5*i1^2*i2*k2^3*j5 + 3*i2^2*k2^3*j5 +
            6*i1^2*j2*k2^3*j5 + 2*i2*j2*k2^3*j5 + 2*j2^2*k2^3*j5 + 4*i1^2*k2^4*j5 + 4*i2*k2^4*j5 + 6*j2*k2^4*j5 +
            6*k2^5*j5 + 3*i1^8*l2*j5 + 3*i1^6*i2*l2*j5 + 6*i1^4*i2^2*l2*j5 + 5*i1^2*i2^3*l2*j5 + 6*i2^4*l2*j5 +
            4*i1^6*j2*l2*j5 + 6*i1^4*i2*j2*l2*j5 + 6*i1^2*i2^2*j2*l2*j5 + i2^3*j2*l2*j5 + 3*i1^4*j2^2*l2*j5 +
            4*i1^2*i2*j2^2*l2*j5 + 3*i2^2*j2^2*l2*j5 + 6*i1^2*j2^3*l2*j5 + 2*i2*j2^3*l2*j5 + 2*j2^4*l2*j5 +
            i1^4*i2*k2*l2*j5 + 6*i1^2*i2^2*k2*l2*j5 + 4*i2^3*k2*l2*j5 + 2*i1^2*i2*j2*k2*l2*j5 + 6*i2^2*j2*k2*l2*j5 +
            i1^2*j2^2*k2*l2*j5 + 4*j2^3*k2*l2*j5 + 2*i1^4*k2^2*l2*j5 + 4*i1^2*i2*k2^2*l2*j5 + 2*i2^2*k2^2*l2*j5 +
            4*i1^2*j2*k2^2*l2*j5 + 2*i2*j2*k2^2*l2*j5 + 2*j2^2*k2^2*l2*j5 + i1^2*k2^3*l2*j5 + i2*k2^3*l2*j5 +
            k2^4*l2*j5 + 3*i1^6*l2^2*j5 + 4*i1^4*i2*l2^2*j5 + 3*i1^2*i2^2*l2^2*j5 + 2*i2^3*l2^2*j5 +
            i1^2*i2*j2*l2^2*j5 + 5*i2^2*j2*l2^2*j5 + 3*i1^2*j2^2*l2^2*j5 + i2*j2^2*l2^2*j5 + 6*j2^3*l2^2*j5 +
            2*i1^4*k2*l2^2*j5 + 2*i1^2*i2*k2*l2^2*j5 + 6*i2^2*k2*l2^2*j5 + 6*i1^2*j2*k2*l2^2*j5 + i2*j2*k2*l2^2*j5 +
            5*j2^2*k2*l2^2*j5 + i1^2*k2^2*l2^2*j5 + 3*j2*k2^2*l2^2*j5 + 3*i1^4*l2^3*j5 + 6*i2^2*l2^3*j5 +
            5*i1^2*j2*l2^3*j5 + 5*i2*j2*l2^3*j5 + 5*j2^2*l2^3*j5 + 2*i1^2*k2*l2^3*j5 + 6*i2*k2*l2^3*j5 + j2*k2*l2^3*j5
            + 2*k2^2*l2^3*j5 + 6*i1^2*l2^4*j5 + i2*l2^4*j5 + j2*l2^4*j5 + 3*k2*l2^4*j5 + 3*l2^5*j5 + 3*i1^3*i2^2*i3*j5
            + 6*i1^3*i2*j2*i3*j5 + 5*i1^3*j2^2*i3*j5 + 2*i1^5*k2*i3*j5 + 5*i1^3*i2*k2*i3*j5 + 5*i1^3*j2*k2*i3*j5 +
            5*i1*i2*j2*k2*i3*j5 + 2*i1*j2^2*k2*i3*j5 + 3*i1^3*k2^2*i3*j5 + 3*i1*i2*k2^2*i3*j5 + 5*i1*j2*k2^2*i3*j5 +
            2*i1*k2^3*i3*j5 + i1^3*i2*l2*i3*j5 + 6*i1*i2^2*l2*i3*j5 + 5*i1^3*j2*l2*i3*j5 + i1*i2*j2*l2*i3*j5 +
            5*i1^3*k2*l2*i3*j5 + 3*i1*i2*k2*l2*i3*j5 + 3*i1*j2*k2*l2*i3*j5 + 2*i1*k2^2*l2*i3*j5 + 4*i1^3*l2^2*i3*j5 +
            3*i1*i2*l2^2*i3*j5 + 5*i1*j2*l2^2*i3*j5 + 2*i1*l2^3*i3*j5 + 5*i1^4*i3^2*j5 + 4*i1^2*i2*i3^2*j5 +
            3*i1^2*j2*i3^2*j5 + 5*i1^2*k2*i3^2*j5 + 3*i2*k2*i3^2*j5 + 4*j2*k2*i3^2*j5 + 2*k2^2*i3^2*j5 +
            3*i1^2*l2*i3^2*j5 + 4*i2*l2*i3^2*j5 + 3*j2*l2*i3^2*j5 + 2*k2*l2*i3^2*j5 + 5*l2^2*i3^2*j5 + 4*i1^7*j3*j5 +
            3*i1^5*i2*j3*j5 + 6*i1*i2^3*j3*j5 + 4*i1^5*j2*j3*j5 + 5*i1^3*i2*j2*j3*j5 + 2*i1^3*j2^2*j3*j5 +
            3*i1*i2*j2^2*j3*j5 + 5*i1*j2^3*j3*j5 + 6*i1^5*k2*j3*j5 + 3*i1*i2^2*k2*j3*j5 + i1^3*j2*k2*j3*j5 +
            3*i1*i2*j2*k2*j3*j5 + 6*i1*j2^2*k2*j3*j5 + 2*i1^3*k2^2*j3*j5 + 3*i1*i2*k2^2*j3*j5 + i1*k2^3*j3*j5 +
            5*i1^5*l2*j3*j5 + 6*i1*i2^2*l2*j3*j5 + 3*i1*i2*j2*l2*j3*j5 + 4*i1*j2^2*l2*j3*j5 + 2*i1^3*k2*l2*j3*j5 +
            3*i1*i2*k2*l2*j3*j5 + 4*i1*j2*k2*l2*j3*j5 + 3*i1*k2^2*l2*j3*j5 + i1^3*l2^2*j3*j5 + 2*i1*i2*l2^2*j3*j5 +
            3*i1*j2*l2^2*j3*j5 + i1*k2*l2^2*j3*j5 + 2*i1*l2^3*j3*j5 + 4*i1^4*i3*j3*j5 + 3*i2^2*i3*j3*j5 +
            6*i1^2*j2*i3*j3*j5 + i2*j2*i3*j3*j5 + 3*j2^2*i3*j3*j5 + 6*i1^2*k2*i3*j3*j5 + 5*i2*k2*i3*j3*j5 +
            6*j2*k2*i3*j3*j5 + 5*k2^2*i3*j3*j5 + i1^2*l2*i3*j3*j5 + 6*i2*l2*i3*j3*j5 + 4*j2*l2*i3*j3*j5 +
            6*k2*l2*i3*j3*j5 + 3*l2^2*i3*j3*j5 + 2*i1*i3^2*j3*j5 + 2*i1^4*j3^2*j5 + 3*i1^2*i2*j3^2*j5 + i2^2*j3^2*j5 +
            6*i1^2*j2*j3^2*j5 + 6*j2^2*j3^2*j5 + 2*i1^2*k2*j3^2*j5 + i2*k2*j3^2*j5 + 6*j2*k2*j3^2*j5 + k2^2*j3^2*j5 +
            3*i1^2*l2*j3^2*j5 + 5*i2*l2*j3^2*j5 + 6*j2*l2*j3^2*j5 + 6*k2*l2*j3^2*j5 + l2^2*j3^2*j5 + 6*i1*j3^3*j5 +
            i1^7*k3*j5 + 2*i1^5*i2*k3*j5 + 3*i1^3*i2^2*k3*j5 + i1^5*j2*k3*j5 + 6*i1^3*i2*j2*k3*j5 + i1*i2^2*j2*k3*j5 +
            6*i1^3*j2^2*k3*j5 + 5*i1*i2*j2^2*k3*j5 + i1*j2^3*k3*j5 + i1^5*k2*k3*j5 + 4*i1^3*i2*k2*k3*j5 +
            6*i1*i2^2*k2*k3*j5 + 6*i1^3*j2*k2*k3*j5 + 3*i1*i2*j2*k2*k3*j5 + 5*i1*j2^2*k2*k3*j5 + 5*i1^3*k2^2*k3*j5 +
            i1*j2*k2^2*k3*j5 + 2*i1*k2^3*k3*j5 + 3*i1^5*l2*k3*j5 + 3*i1^3*i2*l2*k3*j5 + 5*i1*i2^2*l2*k3*j5 +
            6*i1^3*j2*l2*k3*j5 + 4*i1*i2*j2*l2*k3*j5 + 6*i1*j2^2*l2*k3*j5 + i1^3*k2*l2*k3*j5 + i1*i2*k2*l2*k3*j5 +
            4*i1*j2*k2*l2*k3*j5 + 3*i1*k2^2*l2*k3*j5 + i1^3*l2^2*k3*j5 + 4*i1*j2*l2^2*k3*j5 + 2*i1*k2*l2^2*k3*j5 +
            4*i1*l2^3*k3*j5 + 5*i1^4*j3*k3*j5 + 2*i1^2*i2*j3*k3*j5 + 3*i2^2*j3*k3*j5 + 2*i2*j2*j3*k3*j5 +
            2*j2^2*j3*k3*j5 + 4*i2*k2*j3*k3*j5 + j2*k2*j3*k3*j5 + 6*k2^2*j3*k3*j5 + 4*i2*l2*j3*k3*j5 +
            3*k2*l2*j3*k3*j5 + 6*l2^2*j3*k3*j5 + 3*i1*j3^2*k3*j5 + 5*i1^2*i2*k3^2*j5 + 2*i2^2*k3^2*j5 +
            5*i1^2*j2*k3^2*j5 + 4*i2*j2*k3^2*j5 + j2^2*k3^2*j5 + 6*i1^2*k2*k3^2*j5 + 6*i2*k2*k3^2*j5 + 2*k2^2*k3^2*j5
            + 2*i1^2*l2*k3^2*j5 + 4*i2*l2*k3^2*j5 + 2*j2*l2*k3^2*j5 + 5*k2*l2*k3^2*j5 + 2*i1*j3*k3^2*j5 + i1^7*l3*j5 +
            6*i1^5*i2*l3*j5 + 3*i1^3*i2^2*l3*j5 + 4*i1*i2^3*l3*j5 + 2*i1^5*j2*l3*j5 + i1^3*i2*j2*l3*j5 +
            i1*i2^2*j2*l3*j5 + 5*i1^3*j2^2*l3*j5 + 6*i1*i2*j2^2*l3*j5 + 3*i1*j2^3*l3*j5 + 6*i1^3*i2*k2*l3*j5 +
            3*i1*i2^2*k2*l3*j5 + 5*i1^3*j2*k2*l3*j5 + i1*j2^2*k2*l3*j5 + 5*i1^3*k2^2*l3*j5 + 3*i1*i2*k2^2*l3*j5 +
            5*i1^3*j2*l2*l3*j5 + 4*i1*i2*j2*l2*l3*j5 + 5*i1^3*k2*l2*l3*j5 + i1*i2*k2*l2*l3*j5 + 5*i1*j2*k2*l2*l3*j5 +
            3*i1*k2^2*l2*l3*j5 + 4*i1^3*l2^2*l3*j5 + i1*i2*l2^2*l3*j5 + 3*i1*j2*l2^2*l3*j5 + 3*i1*k2*l2^2*l3*j5 +
            2*i1*l2^3*l3*j5 + i1^4*i3*l3*j5 + 6*i1^2*i2*i3*l3*j5 + 6*i2^2*i3*l3*j5 + 3*i1^2*j2*i3*l3*j5 +
            4*i2*j2*i3*l3*j5 + 4*j2^2*i3*l3*j5 + 5*i1^2*k2*i3*l3*j5 + 2*i2*k2*i3*l3*j5 + 3*j2*k2*i3*l3*j5 +
            6*i1^2*l2*i3*l3*j5 + 3*i2*l2*i3*l3*j5 + 4*j2*l2*i3*l3*j5 + 3*k2*l2*i3*l3*j5 + 5*l2^2*i3*l3*j5 +
            3*i1^4*l3^2*j5 + 6*i1^2*i2*l3^2*j5 + 2*i2^2*l3^2*j5 + 5*i1^2*j2*l3^2*j5 + 3*i2*j2*l3^2*j5 + j2^2*l3^2*j5 +
            4*i1^2*k2*l3^2*j5 + 6*i2*k2*l3^2*j5 + 4*j2*k2*l3^2*j5 + 5*k2^2*l3^2*j5 + 5*i1^2*l2*l3^2*j5 +
            5*i2*l2*l3^2*j5 + 3*j2*l2*l3^2*j5 + 2*k2*l2*l3^2*j5 + l2^2*l3^2*j5 + 6*i1*i3*l3^2*j5 + 2*i1*l3^3*j5 +
            i1^6*i4*j5 + 5*i1^4*i2*i4*j5 + 4*i1^2*i2^2*i4*j5 + 2*i2^3*i4*j5 + 4*i1^2*i2*j2*i4*j5 + 4*i2^2*j2*i4*j5 +
            6*i1^2*j2^2*i4*j5 + j2^3*i4*j5 + 6*i1^2*i2*k2*i4*j5 + 2*i2^2*k2*i4*j5 + 5*i1^2*j2*k2*i4*j5 +
            i2*j2*k2*i4*j5 + j2^2*k2*i4*j5 + 3*i1^2*k2^2*i4*j5 + 3*i2*k2^2*i4*j5 + 2*j2*k2^2*i4*j5 + 5*k2^3*i4*j5 +
            5*i1^4*l2*i4*j5 + 6*i2^2*l2*i4*j5 + i1^2*j2*l2*i4*j5 + 3*i2*j2*l2*i4*j5 + 5*j2^2*l2*i4*j5 +
            5*i2*k2*l2*i4*j5 + 2*j2*k2*l2*i4*j5 + 3*i1^2*l2^2*i4*j5 + 2*i2*l2^2*i4*j5 + 6*j2*l2^2*i4*j5 +
            k2*l2^2*i4*j5 + 5*l2^3*i4*j5 + 5*i1*i2*i3*i4*j5 + 6*i1*j2*i3*i4*j5 + 4*i1*k2*i3*i4*j5 + 4*i1*l2*i3*i4*j5 +
            3*i3^2*i4*j5 + 4*i1^3*j3*i4*j5 + 6*i1*i2*j3*i4*j5 + 6*i1*k2*j3*i4*j5 + 5*i1*l2*j3*i4*j5 + 6*i3*j3*i4*j5 +
            5*j3^2*i4*j5 + 5*i1^4*i2*l4*j5 + 6*i1^2*i2^2*l4*j5 + i2^3*l4*j5 + 3*i1^2*i2*j2*l4*j5 + 3*i2^2*j2*l4*j5 +
            5*i2*j2^2*l4*j5 + 5*j2^3*l4*j5 + 2*i1^4*k2*l4*j5 + 6*i1^2*i2*k2*l4*j5 + 3*i2^2*k2*l4*j5 +
            5*i1^2*j2*k2*l4*j5 + i2*j2*k2*l4*j5 + 3*j2^2*k2*l4*j5 + 6*i1^2*k2^2*l4*j5 + 5*i2*k2^2*l4*j5 +
            3*i1^4*l2*l4*j5 + 3*i1^2*i2*l2*l4*j5 + 3*i2^2*l2*l4*j5 + 5*i1^2*j2*l2*l4*j5 + 6*i2*j2*l2*l4*j5 +
            i1^2*k2*l2*l4*j5 + 5*i2*k2*l2*l4*j5 + 2*j2*k2*l2*l4*j5 + 6*k2^2*l2*l4*j5 + 6*i1^2*l2^2*l4*j5 +
            i2*l2^2*l4*j5 + 2*j2*l2^2*l4*j5 + k2*l2^2*l4*j5 + 6*l2^3*l4*j5 + i1^3*j3*l4*j5 + 2*i1*i2*j3*l4*j5 +
            2*i1*j2*j3*l4*j5 + 4*i1*k2*j3*l4*j5 + 5*i1*l2*j3*l4*j5 + 4*j3^2*l4*j5 + 3*i1^3*k3*l4*j5 + 6*i1*i2*k3*l4*j5
            + 6*i1*j2*k3*l4*j5 + 5*i1*k2*k3*l4*j5 + 2*i1*l2*k3*l4*j5 + j3*k3*l4*j5 + 6*k3^2*l4*j5 + 2*i1^3*l3*l4*j5 +
            6*i1*i2*l3*l4*j5 + 6*i1*j2*l3*l4*j5 + 4*i1*l2*l3*l4*j5 + 2*l3^2*l4*j5 + 2*i1^2*l4^2*j5 + 3*i2*l4^2*j5 +
            k2*l4^2*j5 + 2*l2*l4^2*j5 + i1^5*j5^2 + 2*i1^3*i2*j5^2 + 6*i1*i2^2*j5^2 + 4*i1^3*j2*j5^2 + 3*i1*i2*j2*j5^2
            + 2*i1^3*k2*j5^2 + 4*i1*i2*k2*j5^2 + 5*i1*j2*k2*j5^2 + 3*i1*k2^2*j5^2 + 6*i1^3*l2*j5^2 + i1*i2*l2*j5^2 +
            5*i1*k2*l2*j5^2 + 5*i1*l2^2*j5^2 + i1^2*i3*j5^2 + 5*i2*i3*j5^2 + 6*k2*i3*j5^2 + l2*i3*j5^2 +
            3*i1^2*j3*j5^2 + i2*j3*j5^2 + 3*j2*j3*j5^2 + 3*k2*j3*j5^2 + 3*l2*j3*j5^2 + i1^2*k3*j5^2 + 5*i2*k3*j5^2 +
            6*j2*k3*j5^2 + 6*k2*k3*j5^2 + l2*k3*j5^2 + i1^2*l3*j5^2 + 5*i2*l3*j5^2 + 6*j2*l3*j5^2 + 2*k2*l3*j5^2 +
            l2*l3*j5^2 + 5*i1^10*k5 + 4*i1^6*i2^2*k5 + i1^4*i2^3*k5 + 5*i1^8*j2*k5 + 5*i1^6*i2*j2*k5 +
            4*i1^4*i2^2*j2*k5 + 5*i1^6*j2^2*k5 + 3*i1^4*i2*j2^2*k5 + 6*i1^4*j2^3*k5 + i1^8*k2*k5 + i1^6*i2*k2*k5 +
            6*i1^4*i2^2*k2*k5 + 6*i1^2*i2^3*k2*k5 + 6*i1^4*i2*j2*k2*k5 + 3*i1^2*i2^2*j2*k2*k5 + 2*i1^4*j2^2*k2*k5 +
            4*i1^2*i2*j2^2*k2*k5 + i1^2*j2^3*k2*k5 + i1^6*k2^2*k5 + 3*i1^2*i2^2*k2^2*k5 + 5*i2^3*k2^2*k5 +
            2*i1^4*j2*k2^2*k5 + i1^2*i2*j2*k2^2*k5 + 6*i2^2*j2*k2^2*k5 + 3*i1^2*j2^2*k2^2*k5 + i2*j2^2*k2^2*k5 +
            2*j2^3*k2^2*k5 + i1^4*k2^3*k5 + 3*i2^2*k2^3*k5 + 5*i1^2*j2*k2^3*k5 + 3*i2*j2*k2^3*k5 + j2^2*k2^3*k5 +
            3*i1^2*k2^4*k5 + 5*i2*k2^4*k5 + 5*j2*k2^4*k5 + 4*k2^5*k5 + 6*i1^8*l2*k5 + 5*i1^6*i2*l2*k5 +
            2*i1^4*i2^2*l2*k5 + 3*i1^2*i2^3*l2*k5 + 3*i1^6*j2*l2*k5 + 5*i1^2*i2^2*j2*l2*k5 + 5*i1^4*j2^2*l2*k5 +
            2*i1^2*i2*j2^2*l2*k5 + 4*i1^2*j2^3*l2*k5 + 4*i1^6*k2*l2*k5 + i1^4*j2*k2*l2*k5 + 2*i1^2*i2*j2*k2*l2*k5 +
            5*i1^2*j2^2*k2*l2*k5 + i1^4*k2^2*l2*k5 + 6*i1^2*i2*k2^2*l2*k5 + 2*i1^2*j2*k2^2*l2*k5 + 2*i2*j2*k2^2*l2*k5
            + 5*j2^2*k2^2*l2*k5 + 3*i1^2*k2^3*l2*k5 + 2*j2*k2^3*l2*k5 + 2*k2^4*l2*k5 + 5*i1^6*l2^2*k5 +
            3*i1^4*i2*l2^2*k5 + 4*i1^2*i2^2*l2^2*k5 + 4*i2^3*l2^2*k5 + 4*i1^4*j2*l2^2*k5 + 2*i1^2*i2*j2*l2^2*k5 +
            2*i2^2*j2*l2^2*k5 + i1^2*j2^2*l2^2*k5 + 5*i2*j2^2*l2^2*k5 + 3*j2^3*l2^2*k5 + 6*i1^2*i2*k2*l2^2*k5 +
            2*i2^2*k2*l2^2*k5 + 3*i1^2*j2*k2*l2^2*k5 + 6*i2*j2*k2*l2^2*k5 + 6*j2^2*k2*l2^2*k5 + 3*i1^2*k2^2*l2^2*k5 +
            4*i2*k2^2*l2^2*k5 + 3*j2*k2^2*l2^2*k5 + 5*k2^3*l2^2*k5 + 2*i1^4*l2^3*k5 + 3*i1^2*i2*l2^3*k5 +
            4*i2^2*l2^3*k5 + 4*i1^2*j2*l2^3*k5 + 2*i2*j2*l2^3*k5 + j2^2*l2^3*k5 + 6*i1^2*k2*l2^3*k5 + 4*j2*k2*l2^3*k5
            + 3*k2^2*l2^3*k5 + 2*j2*l2^4*k5 + 4*i1^7*j3*k5 + 4*i1^5*i2*j3*k5 + 6*i1^3*i2^2*j3*k5 + i1^5*j2*j3*k5 +
            2*i1^3*i2*j2*j3*k5 + 6*i1^3*j2^2*j3*k5 + 6*i1^5*k2*j3*k5 + 4*i1^3*i2*k2*j3*k5 + 3*i1*i2^2*k2*j3*k5 +
            3*i1^3*j2*k2*j3*k5 + i1*i2*j2*k2*j3*k5 + 3*i1*j2^2*k2*j3*k5 + 6*i1^3*k2^2*j3*k5 + 2*i1*i2*k2^2*j3*k5 +
            3*i1*j2*k2^2*j3*k5 + i1*k2^3*j3*k5 + 5*i1^5*l2*j3*k5 + 5*i1^3*i2*l2*j3*k5 + 5*i1*i2^2*l2*j3*k5 +
            4*i1*i2*j2*l2*j3*k5 + 5*i1*j2^2*l2*j3*k5 + 5*i1*i2*k2*l2*j3*k5 + i1*j2*k2*l2*j3*k5 + 2*i1*i2*l2^2*j3*k5 +
            6*i1*j2*l2^2*j3*k5 + i1*k2*l2^2*j3*k5 + 5*i1*l2^3*j3*k5 + 2*i1^4*j3^2*k5 + 2*i1^2*k2*j3^2*k5 +
            3*i2*k2*j3^2*k5 + 4*j2*k2*j3^2*k5 + 3*k2^2*j3^2*k5 + 5*i1^2*l2*j3^2*k5 + 2*i2*l2*j3^2*k5 + 5*j2*l2*j3^2*k5
            + 5*k2*l2*j3^2*k5 + l2^2*j3^2*k5 + 2*i1*j3^3*k5 + 4*i1^7*k3*k5 + 6*i1^3*i2^2*k3*k5 + 6*i1^5*j2*k3*k5 +
            2*i1^3*i2*j2*k3*k5 + 6*i1^3*j2^2*k3*k5 + 3*i1^5*k2*k3*k5 + 5*i1^3*i2*k2*k3*k5 + 4*i1*i2^2*k2*k3*k5 +
            6*i1*i2*j2*k2*k3*k5 + 4*i1*j2^2*k2*k3*k5 + 3*i1^3*k2^2*k3*k5 + 2*i1*i2*k2^2*k3*k5 + i1*k2^3*k3*k5 +
            6*i1^5*l2*k3*k5 + 2*i1^3*i2*l2*k3*k5 + 4*i1*i2^2*l2*k3*k5 + 4*i1^3*j2*l2*k3*k5 + 6*i1*i2*j2*l2*k3*k5 +
            4*i1*j2^2*l2*k3*k5 + 3*i1^3*k2*l2*k3*k5 + 5*i1*i2*k2*l2*k3*k5 + i1*k2^2*l2*k3*k5 + i1^3*l2^2*k3*k5 +
            i1*i2*l2^2*k3*k5 + 3*i1*j2*l2^2*k3*k5 + 6*i1*l2^3*k3*k5 + 6*i1^2*i2*k3^2*k5 + i1^2*j2*k3^2*k5 +
            6*i1^2*k2*k3^2*k5 + i2*k2*k3^2*k5 + 6*j2*k2*k3^2*k5 + 6*k2^2*k3^2*k5 + 6*i1^2*l2*k3^2*k5 + i2*l2*k3^2*k5 +
            6*j2*l2*k3^2*k5 + 3*k2*l2*k3^2*k5 + 2*l2^2*k3^2*k5 + 5*i1*k3^3*k5 + 5*i1^7*l3*k5 + 2*i1^5*i2*l3*k5 +
            6*i1^3*i2^2*l3*k5 + 3*i1*i2^3*l3*k5 + i1^5*j2*l3*k5 + 5*i1^3*i2*j2*l3*k5 + 5*i1*i2^2*j2*l3*k5 +
            3*i1^3*j2^2*l3*k5 + 2*i1*i2*j2^2*l3*k5 + 4*i1*j2^3*l3*k5 + 5*i1^5*k2*l3*k5 + 5*i1^3*i2*k2*l3*k5 +
            2*i1*i2^2*k2*l3*k5 + 2*i1^3*j2*k2*l3*k5 + 6*i1*i2*j2*k2*l3*k5 + 6*i1*j2^2*k2*l3*k5 + 2*i1^3*k2^2*l3*k5 +
            2*i1*i2*k2^2*l3*k5 + i1*j2*k2^2*l3*k5 + 3*i1*k2^3*l3*k5 + 2*i1^5*l2*l3*k5 + 4*i1^3*i2*l2*l3*k5 +
            4*i1^3*j2*l2*l3*k5 + 3*i1*i2*j2*l2*l3*k5 + 4*i1*j2^2*l2*l3*k5 + i1*i2*k2*l2*l3*k5 + 4*i1^3*l2^2*l3*k5 +
            6*i1*i2*l2^2*l3*k5 + 6*i1*j2*l2^2*l3*k5 + i1*k2*l2^2*l3*k5 + 4*i1^4*j3*l3*k5 + 6*i1^2*i2*j3*l3*k5 +
            5*i2^2*j3*l3*k5 + 5*i1^2*j2*j3*l3*k5 + 4*i2*j2*j3*l3*k5 + 5*j2^2*j3*l3*k5 + 2*i1^2*k2*j3*l3*k5 +
            6*i2*k2*j3*l3*k5 + 6*j2*k2*j3*l3*k5 + 3*i1^2*l2*j3*l3*k5 + 4*j2*l2*j3*l3*k5 + 4*k2*l2*j3*l3*k5 +
            4*l2^2*j3*l3*k5 + 3*i1*j3^2*l3*k5 + 5*i1^6*l4*k5 + 4*i1^4*i2*l4*k5 + 4*i1^2*i2^2*l4*k5 + 6*i1^4*j2*l4*k5 +
            6*i1^2*i2*j2*l4*k5 + 4*i1^2*j2^2*l4*k5 + 2*i1^2*i2*k2*l4*k5 + 2*i2^2*k2*l4*k5 + 2*i1^2*j2*k2*l4*k5 +
            3*i2*j2*k2*l4*k5 + 2*j2^2*k2*l4*k5 + 2*i1^2*k2^2*l4*k5 + 6*i2*k2^2*l4*k5 + j2*k2^2*l4*k5 + 5*k2^3*l4*k5 +
            4*i1^2*i2*l2*l4*k5 + 6*i2^2*l2*l4*k5 + 5*i1^2*j2*l2*l4*k5 + 2*i2*j2*l2*l4*k5 + 6*j2^2*l2*l4*k5 +
            i1^2*k2*l2*l4*k5 + 2*i2*k2*l2*l4*k5 + 4*j2*k2*l2*l4*k5 + k2^2*l2*l4*k5 + 6*i1^2*l2^2*l4*k5 +
            3*i2*l2^2*l4*k5 + 6*j2*l2^2*l4*k5 + k2*l2^2*l4*k5 + 4*l2^3*l4*k5 + 4*i1^3*j3*l4*k5 + 3*i1*k2*j3*l4*k5 +
            2*i1*l2*j3*l4*k5 + 2*j3^2*l4*k5 + 3*i1^3*k3*l4*k5 + 5*i1*i2*k3*l4*k5 + 2*i1*j2*k3*l4*k5 + 6*i1*k2*k3*l4*k5
            + 5*i1*l2*k3*l4*k5 + 4*k3^2*l4*k5 + 6*i1^3*l3*l4*k5 + 5*i1*i2*l3*l4*k5 + 4*i1*j2*l3*l4*k5 +
            3*i1*k2*l3*l4*k5 + 6*i1*l2*l3*l4*k5 + j3*l3*l4*k5 + 3*i1^2*l4^2*k5 + 5*i2*l4^2*k5 + 2*j2*l4^2*k5 +
            6*k2*l4^2*k5 + 4*l2*l4^2*k5 + 3*i1^10*l5 + 2*i1^8*i2*l5 + 5*i1^6*i2^2*l5 + 2*i1^4*i2^3*l5 + i1^2*i2^4*l5 +
            6*i1^8*j2*l5 + 6*i1^6*i2*j2*l5 + 4*i1^4*i2^2*j2*l5 + 3*i1^2*i2^3*j2*l5 + 6*i1^2*i2^2*j2^2*l5 +
            i1^4*j2^3*l5 + 3*i1^2*i2*j2^3*l5 + i1^2*j2^4*l5 + 2*i1^8*k2*l5 + 4*i1^6*i2*k2*l5 + 5*i1^4*i2^2*k2*l5 +
            5*i1^2*i2^3*k2*l5 + i2^4*k2*l5 + 4*i1^6*j2*k2*l5 + i1^4*i2*j2*k2*l5 + 3*i1^2*i2^2*j2*k2*l5 +
            3*i2^3*j2*k2*l5 + 5*i1^4*j2^2*k2*l5 + 6*i2^2*j2^2*k2*l5 + 6*i1^2*j2^3*k2*l5 + 3*i2*j2^3*k2*l5 + j2^4*k2*l5
            + 5*i1^6*k2^2*l5 + 2*i1^4*i2*k2^2*l5 + 6*i1^2*i2^2*k2^2*l5 + i2^3*k2^2*l5 + 2*i1^4*j2*k2^2*l5 +
            4*i1^2*i2*j2*k2^2*l5 + 6*i2^2*j2*k2^2*l5 + 3*i1^2*j2^2*k2^2*l5 + 6*i2*j2^2*k2^2*l5 + j2^3*k2^2*l5 +
            4*i1^4*k2^3*l5 + 4*i1^2*i2*k2^3*l5 + i2^2*k2^3*l5 + 6*i2*j2*k2^3*l5 + 4*j2^2*k2^3*l5 + 4*i1^2*k2^4*l5 +
            5*i2*k2^4*l5 + 2*k2^5*l5 + 6*i1^8*l2*l5 + 5*i1^6*i2*l2*l5 + 3*i1^4*i2^2*l2*l5 + 4*i1^2*i2^3*l2*l5 +
            2*i2^4*l2*l5 + 6*i1^6*j2*l2*l5 + i1^4*i2*j2*l2*l5 + 5*i1^2*i2^2*j2*l2*l5 + 6*i2^3*j2*l2*l5 +
            2*i1^4*j2^2*l2*l5 + 6*i1^2*i2*j2^2*l2*l5 + 5*i2^2*j2^2*l2*l5 + 6*i1^2*j2^3*l2*l5 + 6*i2*j2^3*l2*l5 +
            2*j2^4*l2*l5 + 4*i1^6*k2*l2*l5 + 4*i1^2*i2^2*k2*l2*l5 + 6*i2^3*k2*l2*l5 + 3*i1^4*j2*k2*l2*l5 +
            3*i1^2*i2*j2*k2*l2*l5 + 5*i2^2*j2*k2*l2*l5 + 3*j2^3*k2*l2*l5 + 6*i1^4*k2^2*l2*l5 + 3*i1^2*i2*k2^2*l2*l5 +
            3*i2^2*k2^2*l2*l5 + 6*i1^2*k2^3*l2*l5 + 2*i2*k2^3*l2*l5 + 6*j2*k2^3*l2*l5 + 2*k2^4*l2*l5 + i1^6*l2^2*l5 +
            2*i1^4*i2*l2^2*l5 + 5*i1^2*i2^2*l2^2*l5 + 5*i2^3*l2^2*l5 + 2*i1^4*j2*l2^2*l5 + 4*i1^2*i2*j2*l2^2*l5 +
            i2^2*j2*l2^2*l5 + 4*i1^2*j2^2*l2^2*l5 + 4*i2*j2^2*l2^2*l5 + 4*j2^3*l2^2*l5 + 2*i1^4*k2*l2^2*l5 +
            4*i1^2*i2*k2*l2^2*l5 + 2*i2^2*k2*l2^2*l5 + 6*i1^2*j2*k2*l2^2*l5 + 6*i2*j2*k2*l2^2*l5 + 6*j2^2*k2*l2^2*l5 +
            2*i1^2*k2^2*l2^2*l5 + j2*k2^2*l2^2*l5 + k2^3*l2^2*l5 + 2*i1^4*l2^3*l5 + 5*i1^2*i2*l2^3*l5 + 4*i2^2*l2^3*l5
            + 5*i1^2*j2*l2^3*l5 + 6*j2^2*l2^3*l5 + 3*i1^2*k2*l2^3*l5 + 3*i2*k2*l2^3*l5 + 5*j2*k2*l2^3*l5 +
            k2^2*l2^3*l5 + 4*i1^2*l2^4*l5 + 3*i2*l2^4*l5 + 4*k2*l2^4*l5 + 2*l2^5*l5 + 2*i1^7*i3*l5 + 6*i1^5*i2*i3*l5 +
            6*i1^5*j2*i3*l5 + 3*i1^5*k2*i3*l5 + 5*i1^3*i2*k2*i3*l5 + 4*i1*i2^2*k2*i3*l5 + 4*i1^3*j2*k2*i3*l5 +
            6*i1*i2*j2*k2*i3*l5 + 4*i1*j2^2*k2*i3*l5 + 5*i1^3*k2^2*i3*l5 + 2*i1*k2^3*i3*l5 + 3*i1^5*l2*i3*l5 +
            i1^3*i2*l2*i3*l5 + 6*i1*i2^2*l2*i3*l5 + 5*i1^3*j2*l2*i3*l5 + 2*i1*i2*j2*l2*i3*l5 + 6*i1*j2^2*l2*i3*l5 +
            4*i1^3*k2*l2*i3*l5 + 6*i1*i2*k2*l2*i3*l5 + 3*i1*j2*k2*l2*i3*l5 + 2*i1*k2^2*l2*i3*l5 + 2*i1^3*l2^2*i3*l5 +
            6*i1*i2*l2^2*i3*l5 + 3*i1*j2*l2^2*i3*l5 + 4*i1*l2^3*i3*l5 + 6*i1^4*i3^2*l5 + 5*i1^2*k2*i3^2*l5 +
            3*k2*l2*i3^2*l5 + 2*l2^2*i3^2*l5 + 3*i1^7*j3*l5 + 5*i1^5*i2*j3*l5 + 3*i1^3*i2^2*j3*l5 + 6*i1*i2^3*j3*l5 +
            2*i1^5*j2*j3*l5 + 3*i1*i2^2*j2*j3*l5 + 4*i1^3*j2^2*j3*l5 + 4*i1*i2*j2^2*j3*l5 + i1*j2^3*j3*l5 +
            i1^3*i2*k2*j3*l5 + i1*i2^2*k2*j3*l5 + 2*i1^3*j2*k2*j3*l5 + i1*i2*j2*k2*j3*l5 + 5*i1*j2^2*k2*j3*l5 +
            5*i1^3*k2^2*j3*l5 + i1*j2*k2^2*j3*l5 + 2*i1*k2^3*j3*l5 + 4*i1^5*l2*j3*l5 + i1*i2^2*l2*j3*l5 +
            2*i1^3*j2*l2*j3*l5 + i1*i2*j2*l2*j3*l5 + 5*i1*j2^2*l2*j3*l5 + 5*i1^3*k2*l2*j3*l5 + 5*i1*i2*k2*l2*j3*l5 +
            2*i1*j2*k2*l2*j3*l5 + 6*i1*k2^2*l2*j3*l5 + 3*i1^3*l2^2*j3*l5 + 2*i1*i2*l2^2*j3*l5 + 6*i1*j2*l2^2*j3*l5 +
            5*i1*k2*l2^2*j3*l5 + 5*i1*l2^3*j3*l5 + 6*i1^4*j3^2*l5 + 5*i1^2*i2*j3^2*l5 + 3*i2^2*j3^2*l5 +
            4*i1^2*j2*j3^2*l5 + i2*j2*j3^2*l5 + 3*j2^2*j3^2*l5 + 5*i1^2*k2*j3^2*l5 + 2*i2*k2*j3^2*l5 + 4*j2*k2*j3^2*l5
            + 4*k2^2*j3^2*l5 + 4*i1^2*l2*j3^2*l5 + i2*l2*j3^2*l5 + 6*j2*l2*j3^2*l5 + 6*k2*l2*j3^2*l5 + l2^2*j3^2*l5 +
            3*i1*j3^3*l5 + 4*i1^7*k3*l5 + 5*i1^5*i2*k3*l5 + i1^3*i2^2*k3*l5 + 2*i1*i2^3*k3*l5 + 2*i1^5*j2*k3*l5 +
            3*i1^3*i2*j2*k3*l5 + i1*i2^2*j2*k3*l5 + 3*i1^3*j2^2*k3*l5 + 6*i1*i2*j2^2*k3*l5 + 5*i1*j2^3*k3*l5 +
            2*i1^5*k2*k3*l5 + 2*i1^3*i2*k2*k3*l5 + 3*i1*i2^2*k2*k3*l5 + 6*i1*i2*j2*k2*k3*l5 + 5*i1*j2^2*k2*k3*l5 +
            2*i1^3*k2^2*k3*l5 + i1*i2*k2^2*k3*l5 + i1*j2*k2^2*k3*l5 + 4*i1*k2^3*k3*l5 + 6*i1^3*i2*l2*k3*l5 +
            i1*i2^2*l2*k3*l5 + 4*i1^3*j2*l2*k3*l5 + 6*i1*j2^2*l2*k3*l5 + 2*i1^3*k2*l2*k3*l5 + 3*i1*i2*k2*l2*k3*l5 +
            6*i1*j2*k2*l2*k3*l5 + 3*i1*k2^2*l2*k3*l5 + 2*i1^3*l2^2*k3*l5 + i1*i2*l2^2*k3*l5 + i1*k2*l2^2*k3*l5 +
            i1*l2^3*k3*l5 + 3*i1^4*j3*k3*l5 + i1^2*i2*j3*k3*l5 + 6*i2^2*j3*k3*l5 + 5*i1^2*j2*j3*k3*l5 +
            2*i2*j2*j3*k3*l5 + 6*j2^2*j3*k3*l5 + 2*i1^2*k2*j3*k3*l5 + 6*i2*k2*j3*k3*l5 + j2*k2*j3*k3*l5 +
            k2^2*j3*k3*l5 + i1^2*l2*j3*k3*l5 + 4*i2*l2*j3*k3*l5 + j2*l2*j3*k3*l5 + 6*k2*l2*j3*k3*l5 + 2*l2^2*j3*k3*l5
            + 6*i1*j3^2*k3*l5 + 5*i1^4*k3^2*l5 + 2*i1^2*i2*k3^2*l5 + 3*i1^2*j2*k3^2*l5 + 5*i1^2*k2*k3^2*l5 +
            5*i2*k2*k3^2*l5 + 4*j2*k2*k3^2*l5 + 5*k2^2*k3^2*l5 + 6*i1^2*l2*k3^2*l5 + 2*i2*l2*k3^2*l5 + 3*j2*l2*k3^2*l5
            + 4*k2*l2*k3^2*l5 + 3*l2^2*k3^2*l5 + 5*i1*j3*k3^2*l5 + 3*i1*k3^3*l5 + 4*i1^7*l3*l5 + i1^5*i2*l3*l5 +
            6*i1^3*i2^2*l3*l5 + 5*i1^5*j2*l3*l5 + 6*i1^3*i2*j2*l3*l5 + 6*i1*i2^2*j2*l3*l5 + 3*i1^3*j2^2*l3*l5 +
            2*i1*i2*j2^2*l3*l5 + 6*i1*j2^3*l3*l5 + 4*i1^5*k2*l3*l5 + 6*i1^3*i2*k2*l3*l5 + 6*i1*i2^2*k2*l3*l5 +
            2*i1^3*j2*k2*l3*l5 + 4*i1*i2*j2*k2*l3*l5 + 3*i1*j2^2*k2*l3*l5 + i1^3*k2^2*l3*l5 + 2*i1*i2*k2^2*l3*l5 +
            i1*j2*k2^2*l3*l5 + i1*k2^3*l3*l5 + 2*i1^5*l2*l3*l5 + 3*i1*i2^2*l2*l3*l5 + i1^3*j2*l2*l3*l5 +
            3*i1*i2*j2*l2*l3*l5 + 2*i1*j2^2*l2*l3*l5 + 3*i1*i2*k2*l2*l3*l5 + 2*i1*j2*k2*l2*l3*l5 + 3*i1*k2^2*l2*l3*l5
            + 2*i1^3*l2^2*l3*l5 + 4*i1*j2*l2^2*l3*l5 + 2*i1*k2*l2^2*l3*l5 + 5*i1*l2^3*l3*l5 + 6*i1^4*i3*l3*l5 +
            6*i1^2*i2*i3*l3*l5 + 6*i2^2*i3*l3*l5 + 2*i1^2*j2*i3*l3*l5 + 2*i2*j2*i3*l3*l5 + 6*j2^2*i3*l3*l5 +
            4*i1^2*k2*i3*l3*l5 + 3*j2*k2*i3*l3*l5 + 4*k2^2*i3*l3*l5 + 2*i1^2*l2*i3*l3*l5 + 4*i2*l2*i3*l3*l5 +
            3*j2*l2*i3*l3*l5 + 6*k2*l2*i3*l3*l5 + l2^2*i3*l3*l5 + 6*i1*i3^2*l3*l5 + 6*i1^4*j3*l3*l5 +
            6*i1^2*i2*j3*l3*l5 + i2^2*j3*l3*l5 + 6*i1^2*j2*j3*l3*l5 + 6*j2^2*j3*l3*l5 + 5*i1^2*k2*j3*l3*l5 +
            i2*k2*j3*l3*l5 + 4*j2*k2*j3*l3*l5 + 4*k2^2*j3*l3*l5 + i1^2*l2*j3*l3*l5 + 2*i2*l2*j3*l3*l5 +
            4*k2*l2*j3*l3*l5 + 6*l2^2*j3*l3*l5 + 3*i1^4*l3^2*l5 + 6*i1^2*i2*l3^2*l5 + 4*i2^2*l3^2*l5 +
            4*i1^2*j2*l3^2*l5 + 3*i2*j2*l3^2*l5 + 6*j2^2*l3^2*l5 + 6*i2*k2*l3^2*l5 + j2*k2*l3^2*l5 + 5*i2*l2*l3^2*l5 +
            5*j2*l2*l3^2*l5 + 2*k2*l2*l3^2*l5 + 2*l2^2*l3^2*l5 + 6*i1*i3*l3^2*l5 + 4*i1*j3*l3^2*l5 + 6*i1*l3^3*l5 +
            6*i1^6*k4*l5 + 4*i1^4*i2*k4*l5 + 3*i1^4*j2*k4*l5 + i1^4*k2*k4*l5 + 2*i1^2*i2*k2*k4*l5 + 5*i1^2*j2*k2*k4*l5
            + 6*i1^2*k2^2*k4*l5 + 5*i2*k2^2*k4*l5 + 2*j2*k2^2*k4*l5 + 3*k2^3*k4*l5 + i1^4*l2*k4*l5 + i1^2*i2*l2*k4*l5
            + 6*i1^2*j2*l2*k4*l5 + 5*i1^2*k2*l2*k4*l5 + 6*k2^2*l2*k4*l5 + 5*i2*l2^2*k4*l5 + 2*j2*l2^2*k4*l5 +
            k2*l2^2*k4*l5 + l2^3*k4*l5 + 2*i1^3*k3*k4*l5 + 2*i1*k2*k3*k4*l5 + 5*i1*l2*k3*k4*l5 + 2*i1^6*m4*l5 +
            5*i1^2*i2^2*m4*l5 + 6*i2^3*m4*l5 + 4*i1^2*i2*j2*m4*l5 + 3*i2^2*j2*m4*l5 + 5*i1^2*j2^2*m4*l5 +
            4*i2*j2^2*m4*l5 + j2^3*m4*l5 + 5*i1^4*k2*m4*l5 + 5*i1^2*i2*k2*m4*l5 + 6*i2^2*k2*m4*l5 + i1^2*j2*k2*m4*l5 +
            j2^2*k2*m4*l5 + 2*i1^2*k2^2*m4*l5 + 3*i2*k2^2*m4*l5 + 5*j2*k2^2*m4*l5 + 6*i1^2*i2*l2*m4*l5 +
            i1^2*j2*l2*m4*l5 + 2*i1^2*k2*l2*m4*l5 + 4*i2*k2*l2*m4*l5 + 5*k2^2*l2*m4*l5 + 5*i1^2*l2^2*m4*l5 +
            5*i2*l2^2*m4*l5 + j2*l2^2*m4*l5 + 5*k2*l2^2*m4*l5 + 4*l2^3*m4*l5 + 5*i1^3*i3*m4*l5 + 2*i1*i2*i3*m4*l5 +
            5*i1*j2*i3*m4*l5 + 3*i1*k2*i3*m4*l5 + 3*i1*l2*i3*m4*l5 + 3*i1^3*j3*m4*l5 + 6*i1*i2*j3*m4*l5 +
            2*i1*j2*j3*m4*l5 + 6*i1*l2*j3*m4*l5 + j3^2*m4*l5 + 3*i1*i2*l3*m4*l5 + 5*i1*j2*l3*m4*l5 + 5*i1*l2*l3*m4*l5
            + 4*j3*l3*m4*l5 + 6*l3^2*m4*l5 + 3*i1^5*l5^2 + 6*i1^3*i2*l5^2 + 3*i1*i2^2*l5^2 + i1^3*j2*l5^2 +
            i1*i2*j2*l5^2 + 3*i1*j2^2*l5^2 + 3*i1^3*k2*l5^2 + 2*i1*i2*k2*l5^2 + i1*j2*k2*l5^2 + 2*i1*k2^2*l5^2 +
            4*i1^3*l2*l5^2 + 2*i1*i2*l2*l5^2 + 6*i1*j2*l2*l5^2 + 2*i1*k2*l2*l5^2 + i1*l2^2*l5^2 + 4*i1^2*i3*l5^2 +
            5*k2*i3*l5^2 + 5*l2*i3*l5^2 + 6*i1^2*j3*l5^2 + 3*i2*j3*l5^2 + 4*j2*j3*l5^2 + 2*l2*j3*l5^2 + 3*i1^2*k3*l5^2
            + 3*i2*k3*l5^2 + 4*j2*k3*l5^2 + 5*k2*k3*l5^2 + 4*l2*k3*l5^2 + 6*i1^2*l3*l5^2 + 3*i2*l3*l5^2 + j2*l3*l5^2 +
            k2*l3*l5^2 + i1*m4*l5^2 + l5^3 + 4*i1^7*k2*i6 + 3*i1^5*i2*k2*i6 + 4*i1^5*j2*k2*i6 + 6*i1^5*k2^2*i6 +
            5*i1^3*i2*k2^2*i6 + 2*i1^3*j2*k2^2*i6 + 2*i1*i2*k2^3*i6 + 5*i1*j2*k2^3*i6 + i1*k2^4*i6 + 4*i1^7*l2*i6 +
            2*i1^5*i2*l2*i6 + 5*i1^5*j2*l2*i6 + 3*i1^3*i2*k2*l2*i6 + 4*i1^3*j2*k2*l2*i6 + i1*i2*k2^2*l2*i6 +
            6*i1*j2*k2^2*l2*i6 + 2*i1*k2^3*l2*i6 + i1^5*l2^2*i6 + 2*i1^3*i2*l2^2*i6 + 5*i1^3*j2*l2^2*i6 +
            5*i1^3*k2*l2^2*i6 + 6*i1*i2*k2*l2^2*i6 + i1*j2*k2*l2^2*i6 + 6*i1*k2^2*l2^2*i6 + 6*i1*i2*l2^3*i6 +
            i1*j2*l2^3*i6 + 2*i1*k2*l2^3*i6 + i1^6*j3*i6 + 2*i1^4*k2*j3*i6 + k2^3*j3*i6 + i1^2*k2*l2*j3*i6 +
            3*k2^2*l2*j3*i6 + 3*i1^2*l2^2*j3*i6 + 4*k2*l2^2*j3*i6 + 4*l2^3*j3*i6 + 4*i1^6*k3*i6 + 6*i1^2*k2^2*k3*i6 +
            5*k2^3*k3*i6 + i1^4*l2*k3*i6 + 2*i1^2*k2*l2*k3*i6 + 6*k2^2*l2*k3*i6 + i1^2*l2^2*k3*i6 + 2*k2*l2^2*k3*i6 +
            3*l2^3*k3*i6 + 4*i1^6*l3*i6 + 5*i1^4*k2*l3*i6 + 5*i1^2*i2*k2*l3*i6 + 2*i1^2*j2*k2*l3*i6 +
            3*i1^2*k2^2*l3*i6 + 6*i2*k2^2*l3*i6 + j2*k2^2*l3*i6 + 4*k2^3*l3*i6 + 4*i1^4*l2*l3*i6 + i1^2*i2*l2*l3*i6 +
            6*i1^2*j2*l2*l3*i6 + 4*i1^2*k2*l2*l3*i6 + 2*i2*k2*l2*l3*i6 + 5*j2*k2*l2*l3*i6 + 2*k2^2*l2*l3*i6 +
            4*i1^2*l2^2*l3*i6 + 3*i2*l2^2*l3*i6 + 4*j2*l2^2*l3*i6 + 6*k2*l2^2*l3*i6 + l2^3*l3*i6 + 3*i1^3*j3*l3*i6 +
            3*i1*k2*j3*l3*i6 + i1^5*l4*i6 + 4*i1^3*k2*l4*i6 + 6*i1*k2*l2*l4*i6 + 3*i1*l2^2*l4*i6 + 3*i1^2*l3*l4*i6 +
            2*k2*l3*l4*i6 + 6*l2*l3*l4*i6;
        Append(~JI, J15); Append(~Wght, 15);
    end if;

    return JI, Wght;

end function;

function ShiodaAlgebraicInvariantsChar7(FreeJI, ratsolve)

    FF := Universe(FreeJI);

    P7 := PolynomialRing(FF, 7);
    J2 := P7.1; J7 := P7.2; J8 := P7.3;
    J9 := P7.4; J11 := P7.5; J13 := P7.6; J15 := P7.7;

    if ratsolve eq false or not IsFinite(FF) then
	g := 1; LG := [ FF!1 ];
    else
	Support := [i : i in [1..#FreeJI] | FreeJI[i] ne 0];
	if #Support eq 0 then
	    g := 1;
	else
	    g := Gcd([[3, 4, 5, 6, 10, 14][i] : i in Support]);
	end if;
	if g ne 1 then
	    LG := PowerRepresentativesInFiniteFields(FF, g);
	else
	    LG := [ FF!1 ];
	end if;
    end if;

    JIs := [];
    for L in LG do

	J3, J4, J5, J6, J10, J14 := Explode([L^([ 3, 4, 5, 6, 10, 14][i] div g)*FreeJI[i] : i in [1..#FreeJI]]);

	RES := [];

	// deg. 11
	Append(~RES,
	    J2^4*J3 + 4*J2*J3^3 + 5*J2^2*J3*J4 + 3*J3*J4^2 + 2*J2^3*J5 + 3*J3^2*J5 + 6*J2*J4*J5 + 6*J2*J3*J6 + 5*J5*J6 +
	    4*J2^2*J7 + 3*J4*J7 + 2*J3*J8 + 6*J2*J9);

	// deg. 13
	Append(~RES,
	    J2^5*J3 + J2^2*J3^3 + 5*J2^3*J3*J4 + J3^3*J4 + 3*J2*J3*J4^2 + 2*J2*J3^2*J5 + 3*J2^2*J4*J5 + J3*J4*J6 +
	    2*J2*J5*J6 + 4*J2^3*J7 + 4*J3^2*J7 + 3*J2*J4*J7 + 6*J6*J7 + 4*J2^2*J9 + 5*J4*J9 + 4*J3*J10 + J2*J11);

	// deg. 14
	Append(~RES,
	    J2^7 + 2*J2*J3^4 + 3*J2^5*J4 + J2^2*J3^2*J4 + 6*J3^2*J4^2 + 5*J2*J4^3 + 5*J2^3*J3*J5 + 3*J3^3*J5 +
	    6*J2*J3*J4*J5 + J2^4*J6 + 2*J2*J3^2*J6 + J2^2*J4*J6 + J3*J5*J6 + 2*J2*J6^2 + 6*J3*J4*J7 + 3*J2*J5*J7 +
	    J7^2 + 4*J2^3*J8 + 2*J3^2*J8 + 2*J2*J4*J8 + 2*J2*J3*J9 + 6*J5*J9 + 4*J2^2*J10 + 2*J3*J11);

	// deg. 15
	Append(~RES,
	    J2^6*J3 + 4*J2^3*J3^3 + 5*J2*J3^3*J4 + 3*J2^2*J3*J4^2 + 6*J3*J4^3 + 3*J2^2*J3^2*J5 + 6*J3^2*J4*J5 +
	    5*J2*J4^2*J5 + 5*J2*J3*J5^2 + 6*J2^3*J3*J6 + 4*J4*J5*J6 + 2*J2*J3^2*J7 + 3*J2^2*J4*J7 + 5*J4^2*J7 +
	    4*J3*J5*J7 + J2*J6*J7 + J2^2*J3*J8 + 2*J3*J4*J8 + 3*J2*J5*J8 + 4*J7*J8 + J2^3*J9 + 2*J2*J4*J9 +
	    4*J2*J3*J10 + 5*J5*J10 + 2*J2^2*J11 + 2*J4*J11 + 4*J2*J13);

	// deg. 16
	Append(~RES,
	    J2^8 + 3*J2^2*J3^4 + 4*J3^4*J4 + 5*J2^4*J4^2 + 3*J2*J3^2*J4^2 + 5*J2^2*J4^3 + 6*J4^4 + 6*J2*J3^3*J5 +
	    6*J2^2*J3*J4*J5 + 5*J3*J4^2*J5 + 4*J3^2*J5^2 + J2*J4*J5^2 + J2^5*J6 + 4*J2^2*J3^2*J6 + 5*J2^3*J4*J6 +
	    3*J3^2*J4*J6 + 4*J2*J4^2*J6 + J2*J3*J5*J6 + 3*J5^2*J6 + 2*J2^2*J6^2 + J4*J6^2 + 3*J2^3*J3*J7 + 4*J3^3*J7 +
	    4*J2*J3*J4*J7 + 4*J2^2*J5*J7 + 5*J4*J5*J7 + 5*J2*J7^2 + 4*J2^4*J8 + 2*J2*J3^2*J8 + 4*J2^2*J4*J8 + J4^2*J8
	    + 2*J3*J5*J8 + 4*J2^2*J3*J9 + 3*J3*J4*J9 + 5*J2*J5*J9 + 6*J7*J9 + 4*J2^3*J10 + 6*J3^2*J10 + 2*J2*J4*J10 +
	    3*J2*J3*J11 + 5*J5*J11 + 5*J3*J13);

	// deg. 17
	Append(~RES,
	    J2^7*J3 + 2*J2*J3^5 + 2*J2^2*J3^3*J4 + 4*J3^3*J4^2 + 4*J2*J3*J4^3 + 5*J2^3*J3^2*J5 + 5*J2*J3^2*J4*J5 +
	    4*J2^2*J4^2*J5 + 4*J2^2*J3*J5^2 + J2*J3^3*J6 + 5*J2^2*J3*J4*J6 + 4*J3*J4^2*J6 + 5*J2^3*J5*J6 +
	    4*J2*J4*J5*J6 + J2*J3*J6^2 + 4*J5*J6^2 + J2^2*J3^2*J7 + 3*J2^3*J4*J7 + 2*J3^2*J4*J7 + J2*J4^2*J7 +
	    5*J2*J3*J5*J7 + 5*J2^2*J6*J7 + J3*J7^2 + 2*J2^3*J3*J8 + 4*J2*J3*J4*J8 + 5*J2^2*J5*J8 + 6*J3*J6*J8 +
	    4*J2*J7*J8 + 5*J2*J3^2*J9 + J2^2*J4*J9 + 4*J4^2*J9 + 5*J3*J5*J9 + 6*J2*J6*J9 + 6*J8*J9 + J2^2*J3*J10 +
	    3*J3*J4*J10 + 3*J2*J5*J10 + 4*J7*J10 + 2*J2^3*J11 + J3^2*J11 + 5*J2*J4*J11 + 2*J6*J11 + 2*J2^2*J13 +
	    6*J2*J15);

	// deg. 18
	Append(~RES,
	    J2^9 + J3^6 + 5*J2*J3^4*J4 + 6*J2^2*J3^2*J4^2 + 2*J2^3*J4^3 + 2*J3^2*J4^3 + 3*J2*J4^4 + 5*J2^2*J3^3*J5 +
	    3*J3^3*J4*J5 + 6*J2*J3*J4^2*J5 + 4*J2*J3^2*J5^2 + 5*J2^2*J4*J5^2 + J2^3*J3^2*J6 + 2*J2^4*J4*J6 +
	    3*J2*J3^2*J4*J6 + 2*J2^2*J4^2*J6 + 2*J4^3*J6 + 6*J2^2*J3*J5*J6 + J3*J4*J5*J6 + 3*J2*J5^2*J6 + J2^3*J6^2 +
	    6*J2*J4*J6^2 + 4*J6^3 + 3*J2*J3^3*J7 + J2^2*J3*J4*J7 + 5*J2^3*J5*J7 + 3*J3^2*J5*J7 + J2*J4*J5*J7 +
	    2*J2*J3*J6*J7 + J5*J6*J7 + 2*J2^2*J7^2 + 6*J2^5*J8 + 3*J2^3*J4*J8 + 4*J3^2*J4*J8 + 3*J2*J4^2*J8 +
	    5*J2*J3*J5*J8 + 6*J2^2*J6*J8 + 6*J4*J6*J8 + 6*J3*J7*J8 + 2*J2*J8^2 + 3*J3^3*J9 + 6*J2*J3*J4*J9 +
	    3*J4*J5*J9 + 3*J3*J6*J9 + 4*J2*J7*J9 + J9^2 + 6*J2^4*J10 + 2*J2*J3^2*J10 + J2^2*J4*J10 + 6*J4^2*J10 +
	    3*J3*J5*J10 + 4*J2*J6*J10 + 5*J2^2*J3*J11 + 2*J3*J4*J11 + 3*J2*J5*J11 + 5*J7*J11 + 6*J5*J13 + 4*J3*J15);

	// deg. 18
	Append(~RES,
	    3*J3^6 + J2^7*J4 + 6*J2*J3^4*J4 + J2^3*J4^3 + 4*J3^2*J4^3 + 6*J2*J4^4 + 6*J2^2*J3^3*J5 + 5*J2^3*J3*J4*J5 +
	    2*J3^3*J4*J5 + 5*J2*J3*J4^2*J5 + 6*J2*J3^2*J5^2 + 6*J2^2*J4*J5^2 + 4*J2^3*J3^2*J6 + J2^4*J4*J6 +
	    4*J2^2*J4^2*J6 + 6*J2^2*J3*J5*J6 + 5*J3*J4*J5*J6 + J2*J5^2*J6 + 5*J6^3 + 6*J2*J3^3*J7 + 3*J2^2*J3*J4*J7 +
	    2*J3*J4^2*J7 + J2^3*J5*J7 + 6*J3^2*J5*J7 + 4*J2*J4*J5*J7 + 6*J2*J3*J6*J7 + 4*J5*J6*J7 + 3*J2^2*J7^2 +
	    6*J4*J7^2 + 4*J2^5*J8 + 3*J2^2*J3^2*J8 + 2*J2^3*J4*J8 + 5*J3^2*J4*J8 + 6*J2*J4^2*J8 + J2*J3*J5*J8 +
	    6*J2^2*J6*J8 + 2*J4*J6*J8 + 3*J3*J7*J8 + 4*J2*J8^2 + 2*J3^3*J9 + 3*J2*J3*J4*J9 + 4*J2^2*J5*J9 + 5*J4*J5*J9
	    + 2*J3*J6*J9 + 5*J2*J7*J9 + 3*J9^2 + 4*J2^4*J10 + 6*J2*J3^2*J10 + 2*J2^2*J4*J10 + 5*J4^2*J10 + J3*J5*J10 +
	    2*J2*J6*J10 + 2*J3*J4*J11 + 5*J7*J11 + 2*J2*J3*J13 + 5*J5*J13 + 3*J3*J15);

	// deg. 19
	Append(~RES,
	    J2^8*J3 + J2*J3^3*J4^2 + 2*J2^2*J3*J4^3 + 5*J3*J4^4 + 3*J2*J3^4*J5 + 6*J2^3*J4^2*J5 + J3^2*J4^2*J5 +
	    3*J2*J4^3*J5 + 6*J3^3*J5^2 + 4*J2*J3*J4*J5^2 + 2*J2^2*J3^3*J6 + 6*J2^3*J3*J4*J6 + 3*J3^3*J4*J6 +
	    5*J2*J3*J4^2*J6 + J2^2*J4*J5*J6 + 4*J4^2*J5*J6 + 4*J3*J5^2*J6 + 6*J2^2*J3*J6^2 + 4*J3*J4*J6^2 +
	    2*J2*J5*J6^2 + 6*J3^4*J7 + 6*J2*J3^2*J4*J7 + 6*J2^2*J4^2*J7 + 3*J4^3*J7 + 3*J3*J4*J5*J7 + 2*J2*J5^2*J7 +
	    J2*J4*J6*J7 + 3*J6^2*J7 + 3*J2*J3*J7^2 + 3*J5*J7^2 + J2*J3^3*J8 + 3*J2^3*J5*J8 + 5*J3^2*J5*J8 +
	    6*J2*J4*J5*J8 + 2*J2*J3*J6*J8 + 2*J5*J6*J8 + 6*J2^2*J7*J8 + J4*J7*J8 + 6*J3*J8^2 + 5*J2^3*J4*J9 +
	    J3^2*J4*J9 + 4*J2*J4^2*J9 + J2*J3*J5*J9 + 4*J5^2*J9 + 3*J2^2*J6*J9 + J4*J6*J9 + 4*J3*J7*J9 + 3*J2*J8*J9 +
	    2*J3^3*J10 + 6*J2*J3*J4*J10 + 3*J2^2*J5*J10 + J4*J5*J10 + 6*J3*J6*J10 + 5*J2*J7*J10 + 2*J9*J10 +
	    5*J2*J3^2*J11 + 6*J2^2*J4*J11 + 2*J4^2*J11 + J3*J5*J11 + 6*J2*J6*J11 + J8*J11 + J2^3*J13 + 4*J3^2*J13 +
	    J6*J13 + 2*J2^2*J15 + 2*J4*J15);

	// deg. 20
	Append(~RES,
	    J2^10 + 4*J2*J3^6 + 5*J2^2*J3^4*J4 + 3*J3^4*J4^2 + 3*J2*J3^2*J4^3 + 2*J2^2*J4^4 + J4^5 + 2*J3^5*J5 +
	    2*J2*J3^3*J4*J5 + 6*J2^2*J3*J4^2*J5 + 6*J3*J4^3*J5 + 6*J2^2*J3^2*J5^2 + 5*J3^2*J4*J5^2 + 3*J2*J4^2*J5^2 +
	    6*J2*J3*J5^3 + J2*J3^4*J6 + J2^2*J3^2*J4*J6 + 5*J2^3*J4^2*J6 + 3*J2*J4^3*J6 + 4*J3^3*J5*J6 +
	    5*J2^2*J5^2*J6 + 4*J4*J5^2*J6 + 2*J2^4*J6^2 + 4*J2*J3^2*J6^2 + 4*J2^2*J4*J6^2 + 2*J4^2*J6^2 + 4*J3*J5*J6^2
	    + 5*J2*J6^3 + 3*J2^2*J3^3*J7 + 4*J3^3*J4*J7 + J2*J3*J4^2*J7 + J2^2*J4*J5*J7 + 4*J4^2*J5*J7 + 5*J3*J5^2*J7
	    + 6*J3*J4*J6*J7 + 2*J2*J5*J6*J7 + 5*J2^3*J7^2 + 3*J3^2*J7^2 + J2*J4*J7^2 + 4*J6*J7^2 + 6*J3^4*J8 +
	    J2^4*J4*J8 + J2*J3^2*J4*J8 + 4*J2^2*J4^2*J8 + J4^3*J8 + J2^2*J3*J5*J8 + 4*J3*J4*J5*J8 + 2*J2*J5^2*J8 +
	    5*J2^3*J6*J8 + J2*J4*J6*J8 + 6*J6^2*J8 + J2*J3*J7*J8 + 3*J5*J7*J8 + 5*J4*J8^2 + J2*J3^3*J9 + 4*J3*J4^2*J9
	    + 6*J2^3*J5*J9 + 2*J3^2*J5*J9 + J2*J4*J5*J9 + 6*J2*J3*J6*J9 + 2*J5*J6*J9 + 2*J2^2*J7*J9 + 6*J4*J7*J9 +
	    4*J3*J8*J9 + 3*J2*J9^2 + 2*J2^2*J3^2*J10 + 2*J2^3*J4*J10 + 5*J3^2*J4*J10 + J2*J4^2*J10 + 6*J2*J3*J5*J10 +
	    5*J5^2*J10 + 6*J2^2*J6*J10 + 3*J4*J6*J10 + J3*J7*J10 + 2*J2*J8*J10 + J2^3*J3*J11 + 5*J2*J3*J4*J11 +
	    4*J2^2*J5*J11 + 6*J4*J5*J11 + 6*J3*J6*J11 + 3*J2*J7*J11 + 3*J9*J11 + 5*J2^2*J3*J13 + 2*J2*J5*J13 +
	    2*J7*J13 + 5*J2^3*J14 + 4*J3^2*J14 + J2*J3*J15 + 6*J5*J15);

	// deg. 20
	Append(~RES,
	    2*J2*J3^6 + 4*J2^2*J3^4*J4 + J3^4*J4^2 + J2*J3^2*J4^3 + 5*J3^5*J5 + 3*J2*J3^3*J4*J5 + 5*J2^2*J3*J4^2*J5 +
	    2*J3^2*J4*J5^2 + 4*J2*J4^2*J5^2 + 5*J2*J3*J5^3 + 2*J2*J3^4*J6 + 5*J2^2*J3^2*J4*J6 + J3^3*J5*J6 +
	    2*J2*J3*J4*J5*J6 + 2*J2*J3^2*J6^2 + 5*J2^2*J3^3*J7 + 5*J3^3*J4*J7 + 2*J2*J3*J4^2*J7 + 5*J2*J3^2*J5*J7 +
	    5*J2^2*J4*J5*J7 + 4*J4^2*J5*J7 + J3*J5^2*J7 + J2*J5*J6*J7 + 2*J2^3*J7^2 + 3*J3^2*J7^2 + 3*J2*J4*J7^2 +
	    6*J6*J7^2 + J2^6*J8 + 3*J3^4*J8 + 3*J2^4*J4*J8 + J2*J3^2*J4*J8 + 5*J4^3*J8 + J2^2*J3*J5*J8 + 3*J3*J4*J5*J8
	    + 6*J2*J5^2*J8 + J2^3*J6*J8 + J3^2*J6*J8 + J2*J4*J6*J8 + 2*J6^2*J8 + 4*J2*J3*J7*J8 + 3*J5*J7*J8 +
	    4*J2^2*J8^2 + 2*J4*J8^2 + 4*J2*J3^3*J9 + 3*J2^2*J3*J4*J9 + J3*J4^2*J9 + 2*J2^3*J5*J9 + 6*J3^2*J5*J9 +
	    3*J2*J4*J5*J9 + 5*J2*J3*J6*J9 + 4*J5*J6*J9 + J2^2*J7*J9 + 6*J4*J7*J9 + 2*J3*J8*J9 + 4*J2*J9^2 +
	    2*J2^2*J3^2*J10 + 6*J3^2*J4*J10 + 6*J5^2*J10 + 2*J3*J7*J10 + 4*J2*J8*J10 + 6*J2^3*J3*J11 + 4*J3^3*J11 +
	    5*J2*J3*J4*J11 + 2*J2^2*J5*J11 + 6*J4*J5*J11 + J3*J6*J11 + 5*J2*J7*J11 + J9*J11 + 6*J3*J4*J13 + J2*J5*J13
	    + 6*J7*J13 + J5*J15);

	// deg. 21
	Append(~RES,
	    J2^9*J3 + 4*J3^7 + 2*J2*J3^5*J4 + 6*J2^2*J3^3*J4^2 + J3^3*J4^3 + 4*J2*J3*J4^4 + J2^2*J3^4*J5 + 6*J3^4*J4*J5 +
	    4*J2*J3^2*J4^2*J5 + J2^2*J4^3*J5 + 4*J4^4*J5 + 3*J2*J3^3*J5^2 + 2*J2^2*J3*J4*J5^2 + 5*J3*J4^2*J5^2 +
	    2*J3^2*J5^3 + 4*J2*J4*J5^3 + 3*J2*J3^3*J4*J6 + 2*J2^2*J3*J4^2*J6 + 3*J3*J4^3*J6 + 3*J2^2*J3^2*J5*J6 +
	    4*J2^3*J4*J5*J6 + 2*J3^2*J4*J5*J6 + 4*J2*J4^2*J5*J6 + 5*J2*J3*J5^2*J6 + 4*J5^3*J6 + 4*J2*J3*J4*J6^2 +
	    6*J2^2*J5*J6^2 + 2*J3*J6^3 + J2*J3^4*J7 + 2*J2^2*J3^2*J4*J7 + 5*J2^3*J4^2*J7 + 3*J3^2*J4^2*J7 +
	    3*J2*J4^3*J7 + 3*J3^3*J5*J7 + 6*J2*J3*J4*J5*J7 + J2^2*J5^2*J7 + 5*J4*J5^2*J7 + 2*J2*J3^2*J6*J7 +
	    3*J4^2*J6*J7 + 5*J2*J6^2*J7 + 6*J2^2*J3*J7^2 + 2*J3*J4*J7^2 + 5*J2*J5*J7^2 + 2*J7^3 + 6*J2^2*J3^3*J8 +
	    4*J3^3*J4*J8 + 4*J2*J3*J4^2*J8 + 2*J2*J3^2*J5*J8 + J2^2*J4*J5*J8 + 3*J4^2*J5*J8 + 4*J3*J5^2*J8 +
	    2*J2^2*J3*J6*J8 + 5*J2*J5*J6*J8 + 6*J2^3*J7*J8 + 2*J3^2*J7*J8 + 2*J2*J4*J7*J8 + J6*J7*J8 + 4*J2*J3*J8^2 +
	    5*J3^4*J9 + 6*J2*J3^2*J4*J9 + 6*J2^2*J4^2*J9 + J2^2*J3*J5*J9 + 3*J3*J4*J5*J9 + 3*J2^3*J6*J9 + 5*J3^2*J6*J9
	    + 5*J2*J4*J6*J9 + 3*J2*J3*J7*J9 + 2*J5*J7*J9 + 2*J2^2*J8*J9 + 4*J3*J9^2 + 2*J2*J3^3*J10 + J2^2*J3*J4*J10 +
	    4*J3*J4^2*J10 + J2^3*J5*J10 + 3*J3^2*J5*J10 + 5*J2*J4*J5*J10 + 4*J2*J3*J6*J10 + 2*J5*J6*J10 + J2^2*J7*J10
	    + 4*J4*J7*J10 + 5*J3*J8*J10 + J2^2*J3^2*J11 + J2^3*J4*J11 + 6*J3^2*J4*J11 + J2*J4^2*J11 + J2*J3*J5*J11 +
	    J5^2*J11 + 3*J2^2*J6*J11 + 5*J4*J6*J11 + J3*J7*J11 + 4*J2*J8*J11 + 6*J10*J11 + 3*J2*J3^2*J13 +
	    3*J2^2*J4*J13 + 4*J4^2*J13 + 3*J3*J5*J13 + 6*J2*J6*J13 + 3*J8*J13 + 3*J2^2*J3*J14 + 3*J3*J4*J14 +
	    3*J2*J5*J14 + 4*J2^3*J15 + 6*J3^2*J15 + J2*J4*J15 + J6*J15);

	// deg. 22
	Append(~RES,
	    J2^11 + 2*J2^2*J3^6 + 3*J3^6*J4 + 5*J2*J3^4*J4^2 + 2*J2^2*J3^2*J4^3 + J2^3*J4^4 + 3*J2*J4^5 + 4*J2*J3^5*J5 +
	    J2^2*J3^3*J4*J5 + 6*J3^3*J4^2*J5 + J3^4*J5^2 + 4*J2^2*J4^2*J5^2 + 4*J2^2*J3*J5^3 + 4*J2^2*J3^4*J6 +
	    J3^4*J4*J6 + 6*J2*J3^2*J4^2*J6 + 3*J2^2*J4^3*J6 + 3*J4^4*J6 + 5*J2*J3^3*J5*J6 + 5*J2^2*J3*J4*J5*J6 +
	    2*J3*J4^2*J5*J6 + 5*J2*J4*J5^2*J6 + 4*J2^2*J3^2*J6^2 + 2*J2^3*J4*J6^2 + 4*J3^2*J4*J6^2 + 4*J2*J4^2*J6^2 +
	    2*J5^2*J6^2 + 6*J2^2*J6^3 + J4*J6^3 + 5*J3^5*J7 + 2*J2*J3^3*J4*J7 + 2*J3*J4^3*J7 + 3*J2^2*J3^2*J5*J7 +
	    6*J3^2*J4*J5*J7 + 5*J2*J4^2*J5*J7 + 3*J2*J3*J5^2*J7 + 5*J3^3*J6*J7 + 5*J2*J3*J4*J6*J7 + 4*J2^2*J5*J6*J7 +
	    4*J4*J5*J6*J7 + 6*J3*J6^2*J7 + 6*J2*J3^2*J7^2 + J2^2*J4*J7^2 + J4^2*J7^2 + 3*J3*J5*J7^2 + 5*J2*J6*J7^2 +
	    6*J2*J3^4*J8 + 5*J2^2*J3^2*J4*J8 + 6*J3^2*J4^2*J8 + J2*J4^3*J8 + 2*J3^3*J5*J8 + 2*J2*J3*J4*J5*J8 +
	    J2^2*J5^2*J8 + J2^4*J6*J8 + 6*J2*J3^2*J6*J8 + 4*J4^2*J6*J8 + 6*J3*J5*J6*J8 + 5*J2*J6^2*J8 +
	    4*J2^2*J3*J7*J8 + 3*J3*J4*J7*J8 + 4*J7^2*J8 + 4*J2^3*J8^2 + 3*J3^2*J8^2 + 2*J2*J4*J8^2 + 6*J6*J8^2 +
	    6*J2*J3*J4^2*J9 + 5*J2*J3^2*J5*J9 + 2*J4^2*J5*J9 + 3*J2^2*J3*J6*J9 + J2*J5*J6*J9 + 3*J3^2*J7*J9 +
	    6*J2*J4*J7*J9 + J2*J3*J8*J9 + 6*J5*J8*J9 + 4*J2^2*J9^2 + 2*J4*J9^2 + 5*J3^4*J10 + 3*J2*J3^2*J4*J10 +
	    6*J2^2*J4^2*J10 + 4*J4^3*J10 + J2^2*J3*J5*J10 + 5*J3*J4*J5*J10 + J2*J5^2*J10 + J2^3*J6*J10 + 3*J3^2*J6*J10
	    + 3*J2*J4*J6*J10 + 5*J2*J3*J7*J10 + J5*J7*J10 + J2^2*J8*J10 + 4*J4*J8*J10 + 3*J3*J9*J10 + 3*J2*J10^2 +
	    4*J2*J3^3*J11 + 4*J2^2*J3*J4*J11 + 4*J3*J4^2*J11 + 6*J3^2*J5*J11 + 3*J2*J4*J5*J11 + 5*J2*J3*J6*J11 +
	    J5*J6*J11 + 6*J2^2*J7*J11 + 2*J4*J7*J11 + 5*J3*J8*J11 + 5*J2*J9*J11 + 6*J11^2 + 3*J2^2*J5*J13 + J4*J5*J13
	    + 5*J3*J6*J13 + 4*J2*J7*J13 + 5*J2^4*J14 + 5*J2^2*J4*J14 + 3*J3*J5*J14 + 4*J3*J4*J15 + 6*J2*J5*J15 +
	    5*J7*J15);

	// deg. 22
	Append(~RES,
	    6*J2^2*J3^6 + J2^9*J4 + 5*J2*J3^4*J4^2 + 2*J2^2*J3^2*J4^3 + 6*J2^3*J4^4 + 5*J3^2*J4^4 + 3*J2*J3^5*J5 +
	    2*J2^2*J3^3*J4*J5 + 5*J3^3*J4^2*J5 + J2*J3*J4^3*J5 + 4*J3^4*J5^2 + 2*J2*J3^2*J4*J5^2 + 5*J2^2*J4^2*J5^2 +
	    2*J2^2*J3*J5^3 + J2^2*J3^4*J6 + 4*J3^4*J4*J6 + 6*J2*J3^2*J4^2*J6 + 2*J4^4*J6 + 6*J2*J3^3*J5*J6 +
	    2*J3*J4^2*J5*J6 + 4*J3^2*J5^2*J6 + 2*J2*J4*J5^2*J6 + 4*J2^3*J4*J6^2 + 5*J3^2*J4*J6^2 + J5^2*J6^2 +
	    3*J2^2*J6^3 + J3^5*J7 + 5*J2*J3^3*J4*J7 + J2^2*J3*J4^2*J7 + 3*J3*J4^3*J7 + J2*J4^2*J5*J7 + 3*J2*J3*J5^2*J7
	    + 4*J3^3*J6*J7 + 6*J4*J5*J6*J7 + 2*J3*J6^2*J7 + J2*J3^2*J7^2 + 4*J2^2*J4*J7^2 + 5*J4^2*J7^2 + 2*J3*J5*J7^2
	    + 5*J2*J3^4*J8 + 5*J2^2*J3^2*J4*J8 + 3*J2^3*J4^2*J8 + 6*J3^2*J4^2*J8 + 4*J2*J4^3*J8 + 6*J3^3*J5*J8 +
	    6*J2^2*J5^2*J8 + 4*J2^4*J6*J8 + 4*J2*J3^2*J6*J8 + 3*J4^2*J6*J8 + 5*J2*J6^2*J8 + 5*J2^2*J3*J7*J8 +
	    J3*J4*J7*J8 + 4*J2*J5*J7*J8 + 5*J7^2*J8 + 4*J2^3*J8^2 + 5*J3^2*J8^2 + J2*J4*J8^2 + 4*J3^3*J4*J9 +
	    2*J2*J3*J4^2*J9 + 2*J2*J3^2*J5*J9 + 3*J2^2*J4*J5*J9 + 5*J4^2*J5*J9 + 2*J3*J5^2*J9 + J2^2*J3*J6*J9 +
	    2*J2^3*J7*J9 + J3^2*J7*J9 + 5*J2*J4*J7*J9 + 3*J6*J7*J9 + 6*J2*J3*J8*J9 + 2*J5*J8*J9 + 6*J2^2*J9^2 +
	    3*J4*J9^2 + 4*J3^4*J10 + 6*J2*J3^2*J4*J10 + 3*J2^2*J4^2*J10 + 4*J4^3*J10 + 6*J2^2*J3*J5*J10 +
	    3*J3*J4*J5*J10 + 5*J2*J5^2*J10 + 3*J2^3*J6*J10 + J3^2*J6*J10 + 5*J6^2*J10 + J2*J3*J7*J10 + 6*J5*J7*J10 +
	    6*J2^2*J8*J10 + 5*J4*J8*J10 + 3*J3*J9*J10 + 3*J2*J10^2 + J2*J3^3*J11 + 5*J2^2*J3*J4*J11 + 5*J3*J4^2*J11 +
	    6*J3^2*J5*J11 + 2*J2*J4*J5*J11 + 3*J2*J3*J6*J11 + 2*J5*J6*J11 + 2*J2^2*J7*J11 + 6*J4*J7*J11 + 3*J3*J8*J11
	    + 5*J2*J9*J11 + 4*J11^2 + 5*J3^3*J13 + 5*J2*J3*J4*J13 + 4*J2^2*J5*J13 + 2*J4*J5*J13 + 6*J3*J6*J13 +
	    2*J2*J7*J13 + 2*J9*J13 + J2^2*J3*J15 + 3*J3*J4*J15 + 3*J2*J5*J15 + 5*J7*J15);

	// deg. 23
	Append(~RES,
	    J2^10*J3 + J2*J3^7 + 3*J2^2*J3^5*J4 + 4*J3^5*J4^2 + J2*J3^3*J4^3 + 2*J2^2*J3*J4^4 + J3*J4^5 + J3^6*J5 +
	    J2*J3^4*J4*J5 + 6*J2^2*J3^2*J4^2*J5 + J3^2*J4^3*J5 + 6*J2*J4^4*J5 + 2*J3^3*J4*J5^2 + 3*J2*J3*J4^2*J5^2 +
	    6*J2*J3^5*J6 + 6*J2^2*J3^3*J4*J6 + 6*J3^3*J4^2*J6 + 6*J2*J3*J4^3*J6 + 4*J3^4*J5*J6 + 4*J2^2*J4^2*J5*J6 +
	    J4^3*J5*J6 + J3*J4*J5^2*J6 + 4*J2*J5^3*J6 + 4*J2^2*J3*J4*J6^2 + 6*J3*J4^2*J6^2 + 2*J3^2*J5*J6^2 +
	    2*J2*J4*J5*J6^2 + 6*J2*J3*J6^3 + 2*J5*J6^3 + 4*J2^2*J3^4*J7 + 3*J3^4*J4*J7 + 5*J2*J3^2*J4^2*J7 +
	    4*J2^2*J4^3*J7 + 3*J4^4*J7 + 6*J2^2*J3*J4*J5*J7 + 6*J3*J4^2*J5*J7 + 5*J3^2*J5^2*J7 + 2*J2*J4*J5^2*J7 +
	    6*J2^2*J3^2*J6*J7 + 4*J3^2*J4*J6*J7 + 6*J2*J3*J5*J6*J7 + 4*J5^2*J6*J7 + 6*J2^2*J6^2*J7 + 2*J4*J6^2*J7 +
	    J3^3*J7^2 + 5*J2*J3*J4*J7^2 + 5*J2^2*J5*J7^2 + 4*J4*J5*J7^2 + J2*J7^3 + 4*J3^5*J8 + 5*J2*J3^3*J4*J8 +
	    J2^2*J3*J4^2*J8 + 3*J3*J4^3*J8 + 3*J2^2*J3^2*J5*J8 + 4*J3^2*J4*J5*J8 + 6*J2*J4^2*J5*J8 + 6*J2*J3*J5^2*J8 +
	    2*J3^3*J6*J8 + 2*J2*J3*J4*J6*J8 + 3*J2^2*J5*J6*J8 + 6*J3*J6^2*J8 + 4*J2*J3^2*J7*J8 + 6*J2^2*J4*J7*J8 +
	    2*J4^2*J7*J8 + J3*J5*J7*J8 + 6*J2*J6*J7*J8 + J2*J5*J8^2 + 6*J7*J8^2 + J2*J3^4*J9 + 5*J2^2*J3^2*J4*J9 +
	    5*J3^2*J4^2*J9 + J2*J4^3*J9 + 5*J3^3*J5*J9 + 4*J2*J3*J4*J5*J9 + 3*J2^2*J5^2*J9 + 2*J4*J5^2*J9 +
	    3*J2*J3^2*J6*J9 + J2^2*J4*J6*J9 + 3*J4^2*J6*J9 + 3*J3*J5*J6*J9 + J2*J6^2*J9 + J2^2*J3*J7*J9 +
	    4*J2*J5*J7*J9 + 5*J7^2*J9 + 2*J3^2*J8*J9 + 2*J2*J4*J8*J9 + J6*J8*J9 + 4*J2*J3*J9^2 + 5*J5*J9^2 +
	    2*J2^2*J3^3*J10 + 2*J3^3*J4*J10 + J2*J3*J4^2*J10 + 5*J2*J3^2*J5*J10 + 6*J2^2*J4*J5*J10 + 4*J4^2*J5*J10 +
	    6*J3*J5^2*J10 + 4*J2^2*J3*J6*J10 + 3*J3*J4*J6*J10 + 5*J3^2*J7*J10 + 4*J2*J4*J7*J10 + 4*J6*J7*J10 +
	    6*J2*J3*J8*J10 + 3*J5*J8*J10 + 2*J2^2*J9*J10 + 2*J3*J10^2 + 4*J3^4*J11 + 5*J2*J3^2*J4*J11 +
	    3*J2^2*J4^2*J11 + J4^3*J11 + 3*J2^2*J3*J5*J11 + J3*J4*J5*J11 + 2*J2*J5^2*J11 + J2^3*J6*J11 + 5*J3^2*J6*J11
	    + 6*J2*J4*J6*J11 + 5*J6^2*J11 + 2*J2*J3*J7*J11 + 5*J5*J7*J11 + 2*J2^2*J8*J11 + 4*J4*J8*J11 + 6*J3*J9*J11 +
	    4*J2*J10*J11 + 5*J2^2*J3^2*J13 + 4*J3^2*J4*J13 + 6*J2*J3*J5*J13 + 5*J5^2*J13 + 2*J2^2*J6*J13 + 5*J4*J6*J13
	    + 6*J3*J7*J13 + J2*J8*J13 + 3*J10*J13 + 2*J3^3*J14 + 4*J2*J3*J4*J14 + 2*J2^2*J5*J14 + 2*J4*J5*J14 +
	    2*J3*J6*J14 + 5*J2*J7*J14 + 4*J2*J3^2*J15 + 2*J2^2*J4*J15 + 4*J4^2*J15 + J3*J5*J15 + 2*J2*J6*J15 +
	    5*J8*J15);

	// deg. 24
	Append(~RES,
	    J2^12 + 5*J3^8 + 6*J2*J3^6*J4 + 5*J2*J3^2*J4^4 + 4*J2^2*J4^5 + 3*J4^6 + J3^5*J4*J5 + 6*J2*J3^3*J4^2*J5 +
	    2*J2^2*J3*J4^3*J5 + 5*J2*J3^4*J5^2 + 6*J2^2*J3^2*J4*J5^2 + 6*J3^2*J4^2*J5^2 + 3*J2*J4^3*J5^2 + 4*J3^3*J5^3
	    + 2*J2*J3*J4*J5^3 + 6*J2^2*J5^4 + 5*J3^6*J6 + 5*J2*J3^4*J4*J6 + J2^2*J3^2*J4^2*J6 + 6*J2^3*J4^3*J6 +
	    6*J3^2*J4^3*J6 + 2*J2*J4^4*J6 + 4*J2^2*J3^3*J5*J6 + 4*J3^3*J4*J5*J6 + 5*J2*J3*J4^2*J5*J6 +
	    6*J2*J3^2*J5^2*J6 + 2*J2^2*J4*J5^2*J6 + 3*J4^2*J5^2*J6 + 3*J3*J5^3*J6 + J2*J3^2*J4*J6^2 + 2*J2^2*J4^2*J6^2
	    + 2*J4^3*J6^2 + 3*J2^2*J3*J5*J6^2 + 3*J3*J4*J5*J6^2 + J2*J5^2*J6^2 + J2^3*J6^3 + 6*J3^2*J6^3 +
	    3*J2*J4*J6^3 + 6*J6^4 + 3*J2^2*J3^3*J4*J7 + 2*J3^3*J4^2*J7 + 6*J2*J3^2*J4*J5*J7 + 4*J2^2*J4^2*J5*J7 +
	    4*J4^3*J5*J7 + 4*J2^2*J3*J5^2*J7 + 4*J3*J4*J5^2*J7 + 6*J2*J5^3*J7 + 2*J2*J3^3*J6*J7 + J2^2*J3*J4*J6*J7 +
	    2*J3*J4^2*J6*J7 + 3*J2*J4*J5*J6*J7 + 5*J2*J3*J6^2*J7 + 5*J5*J6^2*J7 + 3*J2^2*J3^2*J7^2 + J3^2*J4*J7^2 +
	    3*J2*J4^2*J7^2 + 2*J2*J3*J5*J7^2 + 3*J5^2*J7^2 + 3*J2^2*J6*J7^2 + J4*J6*J7^2 + 3*J3*J7^3 + 4*J2^2*J3^4*J8
	    + 4*J3^4*J4*J8 + 4*J2*J3^2*J4^2*J8 + J4^4*J8 + 6*J2*J3^3*J5*J8 + 3*J3*J4^2*J5*J8 + 4*J3^2*J5^2*J8 +
	    5*J2*J4*J5^2*J8 + 2*J2^2*J3^2*J6*J8 + 2*J2^3*J4*J6*J8 + 4*J3^2*J4*J6*J8 + 5*J2*J4^2*J6*J8 +
	    6*J2*J3*J5*J6*J8 + 5*J5^2*J6*J8 + 3*J4*J6^2*J8 + 4*J3^3*J7*J8 + J2*J3*J4*J7*J8 + 4*J2^2*J5*J7*J8 +
	    5*J4*J5*J7*J8 + 6*J3*J6*J7*J8 + 3*J2*J7^2*J8 + 4*J2^4*J8^2 + 2*J2^2*J4*J8^2 + 6*J4^2*J8^2 + J3*J5*J8^2 +
	    6*J2*J6*J8^2 + 4*J8^3 + J3^5*J9 + 5*J2*J3^3*J4*J9 + 3*J2^2*J3*J4^2*J9 + 4*J3*J4^3*J9 + 5*J3^2*J4*J5*J9 +
	    4*J2*J4^2*J5*J9 + 3*J2*J3*J5^2*J9 + 6*J5^3*J9 + 2*J3^3*J6*J9 + J2*J3*J4*J6*J9 + 6*J2^2*J5*J6*J9 +
	    6*J4*J5*J6*J9 + J3*J6^2*J9 + 5*J2*J3^2*J7*J9 + 5*J2^2*J4*J7*J9 + 3*J3*J5*J7*J9 + 3*J2*J6*J7*J9 +
	    5*J2^2*J3*J8*J9 + 4*J3*J4*J8*J9 + 5*J2*J5*J8*J9 + 5*J3^2*J9^2 + 4*J2*J4*J9^2 + 5*J6*J9^2 + J2*J3^4*J10 +
	    4*J2^2*J3^2*J4*J10 + J2*J4^3*J10 + 6*J3^3*J5*J10 + 4*J2*J3*J4*J5*J10 + 3*J2^2*J5^2*J10 + 4*J4*J5^2*J10 +
	    6*J2^4*J6*J10 + 4*J2^2*J4*J6*J10 + 4*J4^2*J6*J10 + 3*J3*J5*J6*J10 + 4*J2*J6^2*J10 + 2*J2^2*J3*J7*J10 +
	    2*J3*J4*J7*J10 + 3*J2*J5*J7*J10 + 3*J2^3*J8*J10 + 3*J3^2*J8*J10 + 4*J2*J4*J8*J10 + 2*J6*J8*J10 +
	    5*J2*J3*J9*J10 + 2*J5*J9*J10 + 3*J2^2*J10^2 + 6*J2^2*J3^3*J11 + 6*J3^3*J4*J11 + 3*J2*J3^2*J5*J11 +
	    J2^2*J4*J5*J11 + 5*J4^2*J5*J11 + J3*J5^2*J11 + 4*J2^2*J3*J6*J11 + 4*J3*J4*J6*J11 + 4*J2*J5*J6*J11 +
	    5*J3^2*J7*J11 + 5*J2*J4*J7*J11 + 6*J6*J7*J11 + 3*J2*J3*J8*J11 + 5*J5*J8*J11 + J2*J11^2 + 6*J2*J3^3*J13 +
	    6*J2^2*J3*J4*J13 + 2*J3*J4^2*J13 + J3^2*J5*J13 + 6*J2*J4*J5*J13 + 2*J5*J6*J13 + 4*J2^2*J7*J13 +
	    4*J4*J7*J13 + 6*J3*J8*J13 + 4*J2*J9*J13 + 3*J11*J13 + 5*J2^2*J3^2*J14 + 5*J2^3*J4*J14 + 2*J3^2*J4*J14 +
	    4*J2*J4^2*J14 + 2*J2*J3*J5*J14 + 2*J5^2*J14 + 2*J2^2*J6*J14 + 3*J3*J7*J14 + 3*J3^3*J15 + 2*J2*J3*J4*J15 +
	    6*J4*J5*J15 + 4*J3*J6*J15 + 5*J2*J7*J15);

	// deg. 24
	Append(~RES,
	    4*J3^8 + J2^10*J4 + 2*J2*J3^6*J4 + 5*J3^4*J4^3 + 5*J2*J3^2*J4^4 + 6*J2^2*J4^5 + 6*J4^6 + J3^5*J4*J5 +
	    2*J2^2*J3*J4^3*J5 + 5*J3*J4^4*J5 + 4*J2*J3^4*J5^2 + J2^2*J3^2*J4*J5^2 + 3*J3^2*J4^2*J5^2 + 2*J2*J4^3*J5^2
	    + J3^3*J5^3 + 3*J2*J3*J4*J5^3 + 5*J2^2*J5^4 + 5*J3^6*J6 + J2*J3^4*J4*J6 + J2^2*J3^2*J4^2*J6 +
	    3*J2^3*J4^3*J6 + 4*J3^2*J4^3*J6 + 4*J2*J4^4*J6 + 2*J2^2*J3^3*J5*J6 + 4*J3^3*J4*J5*J6 + 2*J2*J3*J4^2*J5*J6
	    + 3*J2*J3^2*J5^2*J6 + 6*J2^2*J4*J5^2*J6 + 5*J4^2*J5^2*J6 + 2*J3*J5^3*J6 + 6*J2*J3^2*J4*J6^2 +
	    4*J2^2*J4^2*J6^2 + 6*J4^3*J6^2 + 6*J2^2*J3*J5*J6^2 + 2*J3*J4*J5*J6^2 + 6*J2*J5^2*J6^2 + J2^3*J6^3 +
	    2*J3^2*J6^3 + J2*J4*J6^3 + 6*J6^4 + 2*J2*J3^5*J7 + 6*J2^2*J3^3*J4*J7 + 4*J3^3*J4^2*J7 + 2*J2*J3*J4^3*J7 +
	    6*J3^4*J5*J7 + J2*J3^2*J4*J5*J7 + 2*J2^2*J4^2*J5*J7 + J4^3*J5*J7 + 2*J2^2*J3*J5^2*J7 + 4*J3*J4*J5^2*J7 +
	    3*J2*J5^3*J7 + J2^2*J3*J4*J6*J7 + J3*J4^2*J6*J7 + J3^2*J5*J6*J7 + J2*J4*J5*J6*J7 + J2*J3*J6^2*J7 +
	    3*J5*J6^2*J7 + J2^2*J3^2*J7^2 + 5*J3^2*J4*J7^2 + 2*J2*J4^2*J7^2 + 5*J2*J3*J5*J7^2 + 3*J5^2*J7^2 +
	    6*J2^2*J6*J7^2 + 4*J4*J6*J7^2 + J3*J7^3 + J2^2*J3^4*J8 + 6*J3^4*J4*J8 + 6*J2*J3^2*J4^2*J8 + J2^2*J4^3*J8 +
	    3*J4^4*J8 + 3*J2*J3^3*J5*J8 + 2*J2^2*J3*J4*J5*J8 + 4*J3*J4^2*J5*J8 + 5*J3^2*J5^2*J8 + 3*J2^2*J3^2*J6*J8 +
	    5*J2^3*J4*J6*J8 + 5*J3^2*J4*J6*J8 + 5*J2*J4^2*J6*J8 + J2*J3*J5*J6*J8 + 2*J5^2*J6*J8 + 3*J2^2*J6^2*J8 +
	    5*J3^3*J7*J8 + 2*J2*J3*J4*J7*J8 + J3*J6*J7*J8 + 4*J2^4*J8^2 + 2*J2*J3^2*J8^2 + J2^2*J4*J8^2 + 4*J4^2*J8^2
	    + 3*J3*J5*J8^2 + 2*J2*J6*J8^2 + J8^3 + 5*J3^5*J9 + 4*J2*J3^3*J4*J9 + J3*J4^3*J9 + 4*J2^2*J3^2*J5*J9 +
	    J3^2*J4*J5*J9 + 4*J2*J4^2*J5*J9 + 5*J2*J3*J5^2*J9 + J5^3*J9 + 6*J3^3*J6*J9 + 3*J2*J3*J4*J6*J9 +
	    2*J2^2*J5*J6*J9 + 4*J4*J5*J6*J9 + J3*J6^2*J9 + 4*J2*J3^2*J7*J9 + 5*J2^2*J4*J7*J9 + 5*J4^2*J7*J9 +
	    J3*J5*J7*J9 + 6*J2*J6*J7*J9 + 3*J2^2*J3*J8*J9 + 4*J3*J4*J8*J9 + 2*J7*J8*J9 + 4*J3^2*J9^2 + 4*J2*J4*J9^2 +
	    5*J6*J9^2 + 5*J2*J3^4*J10 + 3*J2^3*J4^2*J10 + 2*J3^2*J4^2*J10 + 4*J2*J4^3*J10 + 2*J3^3*J5*J10 +
	    J2*J3*J4*J5*J10 + J2^2*J5^2*J10 + 3*J4*J5^2*J10 + 4*J2^4*J6*J10 + 4*J2*J3^2*J6*J10 + 5*J2^2*J4*J6*J10 +
	    6*J4^2*J6*J10 + J2*J6^2*J10 + 6*J2^2*J3*J7*J10 + 5*J2*J5*J7*J10 + 5*J7^2*J10 + 2*J2^3*J8*J10 +
	    6*J3^2*J8*J10 + 3*J2*J4*J8*J10 + J2*J3*J9*J10 + 4*J5*J9*J10 + 4*J2^2*J10^2 + 2*J4*J10^2 + 6*J2^2*J3^3*J11
	    + 3*J3^3*J4*J11 + 3*J2*J3^2*J5*J11 + J2^2*J4*J5*J11 + 2*J4^2*J5*J11 + 4*J2^2*J3*J6*J11 + 2*J3*J4*J6*J11 +
	    J2*J5*J6*J11 + 5*J3^2*J7*J11 + 3*J2*J4*J7*J11 + 4*J6*J7*J11 + 3*J5*J8*J11 + 3*J4*J9*J11 + 6*J3*J10*J11 +
	    6*J2*J3^3*J13 + 6*J2^2*J3*J4*J13 + 3*J3^2*J5*J13 + 2*J2*J3*J6*J13 + 5*J5*J6*J13 + 3*J2^2*J7*J13 +
	    2*J4*J7*J13 + J3*J8*J13 + J11*J13 + 2*J2^2*J3^2*J14 + 6*J2^3*J4*J14 + 6*J3^2*J4*J14 + 6*J2*J4^2*J14 +
	    5*J2*J3*J5*J14 + 4*J5^2*J14 + 6*J2^2*J6*J14 + 2*J3*J7*J14 + 2*J3^3*J15 + 5*J2*J3*J4*J15 + 4*J2^2*J5*J15 +
	    J4*J5*J15 + 3*J3*J6*J15 + 5*J9*J15);

	// deg. 26
	Append(~RES,
	    J2^13 + 6*J3^6*J4^2 + 3*J2*J3^4*J4^3 + 2*J2*J4^6 + 5*J3^7*J5 + 3*J2*J3^5*J4*J5 + 4*J3^3*J4^3*J5 +
	    J2*J3*J4^4*J5 + J3^4*J4*J5^2 + 6*J2*J3^2*J4^2*J5^2 + 5*J2^2*J4^3*J5^2 + 3*J4^4*J5^2 + 5*J2*J3^3*J5^3 +
	    4*J2^2*J3*J4*J5^3 + 3*J3^2*J5^4 + 6*J2*J4*J5^4 + 6*J2*J3^6*J6 + 2*J3^4*J4^2*J6 + 5*J2*J3^2*J4^3*J6 +
	    4*J2^2*J4^4*J6 + 5*J3^5*J5*J6 + 5*J2*J3^3*J4*J5*J6 + 4*J2^2*J3*J4^2*J5*J6 + 2*J3*J4^3*J5*J6 +
	    5*J2^2*J3^2*J5^2*J6 + 5*J3^2*J4*J5^2*J6 + 4*J2*J4^2*J5^2*J6 + 5*J2*J3*J5^3*J6 + 2*J5^4*J6 + J2*J3^4*J6^2 +
	    2*J2^3*J4^2*J6^2 + 4*J3^2*J4^2*J6^2 + 5*J2*J4^3*J6^2 + J3^3*J5*J6^2 + J2*J3*J4*J5*J6^2 + 4*J2^2*J5^2*J6^2
	    + 3*J4*J5^2*J6^2 + 3*J2*J3^2*J6^3 + J2^2*J4*J6^3 + J4^2*J6^3 + 3*J2*J6^4 + 3*J3^5*J4*J7 +
	    2*J2*J3^3*J4^2*J7 + 3*J3*J4^4*J7 + 3*J2*J3^4*J5*J7 + J3^2*J4^2*J5*J7 + 5*J2*J4^3*J5*J7 + J3^3*J5^2*J7 +
	    2*J2*J3*J4*J5^2*J7 + 2*J2^2*J5^3*J7 + 3*J3^3*J4*J6*J7 + 2*J2*J3*J4^2*J6*J7 + 2*J2*J3^2*J5*J6*J7 +
	    6*J2^2*J4*J5*J6*J7 + 5*J4^2*J5*J6*J7 + 4*J3*J5^2*J6*J7 + J3*J4*J6^2*J7 + 5*J2*J5*J6^2*J7 + 2*J3^4*J7^2 +
	    2*J2*J3^2*J4*J7^2 + 6*J2^2*J4^2*J7^2 + J4^3*J7^2 + 4*J3*J4*J5*J7^2 + 2*J2*J5^2*J7^2 + 2*J3^2*J6*J7^2 +
	    5*J2*J4*J6*J7^2 + J6^2*J7^2 + 3*J2*J3*J7^3 + 2*J3^6*J8 + 4*J2*J3^4*J4*J8 + 5*J2^3*J4^3*J8 + 2*J3^2*J4^3*J8
	    + J2*J4^4*J8 + 2*J2*J3*J4^2*J5*J8 + 3*J2*J3^2*J5^2*J8 + 2*J2^2*J4*J5^2*J8 + 4*J4^2*J5^2*J8 + 4*J3*J5^3*J8
	    + 2*J3^4*J6*J8 + 6*J2*J3^2*J4*J6*J8 + 3*J2^2*J4^2*J6*J8 + J4^3*J6*J8 + 3*J2^2*J3*J5*J6*J8 +
	    6*J3*J4*J5*J6*J8 + 5*J2*J5^2*J6*J8 + 4*J2^3*J6^2*J8 + 6*J3^2*J6^2*J8 + 6*J2*J4*J6^2*J8 + 4*J2*J3^3*J7*J8 +
	    2*J2^2*J3*J4*J7*J8 + 6*J3*J4^2*J7*J8 + 2*J3^2*J5*J7*J8 + 4*J2*J4*J5*J7*J8 + 2*J2*J3*J6*J7*J8 +
	    J2^2*J7^2*J8 + J4*J7^2*J8 + 6*J2^2*J3^2*J8^2 + 3*J2^3*J4*J8^2 + 4*J3^2*J4*J8^2 + 6*J2*J4^2*J8^2 +
	    4*J2*J3*J5*J8^2 + 6*J2^2*J6*J8^2 + J4*J6*J8^2 + 6*J2*J8^3 + J2*J3^5*J9 + 6*J2*J3*J4^3*J9 + 5*J3^4*J5*J9 +
	    4*J2*J3^2*J4*J5*J9 + J4^3*J5*J9 + 3*J2^2*J3*J5^2*J9 + 6*J3*J4*J5^2*J9 + J2*J5^3*J9 + 5*J2*J3^3*J6*J9 +
	    2*J2^2*J3*J4*J6*J9 + 3*J3*J4^2*J6*J9 + 5*J2*J4*J5*J6*J9 + 2*J2^2*J3^2*J7*J9 + 2*J3^2*J4*J7*J9 +
	    2*J2*J4^2*J7*J9 + 2*J2*J3*J5*J7*J9 + 5*J5^2*J7*J9 + 2*J4*J6*J7*J9 + 2*J3*J7^2*J9 + 4*J3^3*J8*J9 +
	    6*J2*J3*J4*J8*J9 + J2^2*J5*J8*J9 + 5*J4*J5*J8*J9 + 6*J3*J6*J8*J9 + 4*J2*J7*J8*J9 + 2*J2^2*J4*J9^2 +
	    6*J4^2*J9^2 + J3*J5*J9^2 + 3*J2*J6*J9^2 + 6*J8*J9^2 + 4*J3^4*J4*J10 + 4*J2*J3^2*J4^2*J10 + J2^2*J4^3*J10 +
	    5*J2*J3^3*J5*J10 + 6*J2^2*J3*J4*J5*J10 + J3*J4^2*J5*J10 + 3*J2*J4*J5^2*J10 + 5*J2^2*J3^2*J6*J10 +
	    5*J2^3*J4*J6*J10 + 5*J2*J4^2*J6*J10 + 6*J2*J3*J5*J6*J10 + J2^2*J6^2*J10 + 4*J4*J6^2*J10 + J2*J3*J4*J7*J10
	    + 2*J2^2*J5*J7*J10 + 6*J4*J5*J7*J10 + 6*J3*J6*J7*J10 + 5*J2^4*J8*J10 + 2*J2*J3^2*J8*J10 + J2^2*J4*J8*J10 +
	    3*J4^2*J8*J10 + 6*J3*J5*J8*J10 + 6*J2*J6*J8*J10 + 2*J8^2*J10 + 3*J2^2*J3*J9*J10 + 3*J3*J4*J9*J10 +
	    J2*J5*J9*J10 + J7*J9*J10 + 6*J2^3*J10^2 + 5*J3^2*J10^2 + 6*J2*J4*J10^2 + 5*J6*J10^2 + 4*J3^5*J11 +
	    5*J2^2*J3*J4^2*J11 + J3*J4^3*J11 + J3^2*J4*J5*J11 + 5*J2*J4^2*J5*J11 + 6*J2*J3*J5^2*J11 + 6*J5^3*J11 +
	    5*J3^3*J6*J11 + 3*J2*J3*J4*J6*J11 + J2^2*J5*J6*J11 + 5*J4*J5*J6*J11 + J3*J6^2*J11 + 4*J2*J3^2*J7*J11 +
	    5*J2^2*J4*J7*J11 + 6*J4^2*J7*J11 + 4*J3*J5*J7*J11 + 4*J2*J6*J7*J11 + 3*J2^2*J3*J8*J11 + 4*J3*J4*J8*J11 +
	    J7*J8*J11 + 3*J3^2*J9*J11 + 2*J6*J9*J11 + 4*J2*J3*J10*J11 + J5*J10*J11 + J2^2*J11^2 + 6*J4*J11^2 +
	    6*J3^3*J4*J13 + 2*J2*J3*J4^2*J13 + 3*J4^2*J5*J13 + 6*J3*J5^2*J13 + 4*J3*J4*J6*J13 + 3*J3^2*J7*J13 +
	    4*J2*J4*J7*J13 + 6*J6*J7*J13 + 2*J2*J3*J8*J13 + 5*J5*J8*J13 + 5*J2^2*J9*J13 + 6*J4*J9*J13 + 3*J2*J11*J13 +
	    3*J13^2 + 6*J2*J3^2*J4*J14 + 4*J2^2*J4^2*J14 + 5*J4^3*J14 + 5*J2^2*J3*J5*J14 + 4*J3*J4*J5*J14 +
	    6*J2*J5^2*J14 + J2^3*J6*J14 + 3*J3^2*J6*J14 + 5*J2*J4*J6*J14 + 5*J2*J3*J7*J14 + 2*J5*J7*J14 +
	    3*J2^2*J8*J14 + 4*J3*J9*J14 + 5*J2*J3^3*J15 + 4*J2^2*J3*J4*J15 + 6*J3^2*J5*J15 + J2*J4*J5*J15 +
	    5*J2*J3*J6*J15 + 4*J5*J6*J15 + J2^2*J7*J15 + 4*J4*J7*J15 + 5*J3*J8*J15 + 3*J2*J9*J15 + 5*J11*J15);

	// deg. 26
	Append(~RES,
	    6*J2*J3^8 + 6*J3^6*J4^2 + 5*J2*J3^4*J4^3 + 5*J2^3*J4^5 + 2*J3^2*J4^5 + 5*J2*J4^6 + 5*J3^7*J5 + 3*J2*J3^5*J4*J5
	    + 6*J3^3*J4^3*J5 + 3*J2*J3*J4^4*J5 + J3^4*J4*J5^2 + 3*J2*J3^2*J4^2*J5^2 + J2^2*J4^3*J5^2 + 5*J4^4*J5^2 +
	    J2*J3^3*J5^3 + 5*J2^2*J3*J4*J5^3 + J3*J4^2*J5^3 + 6*J3^2*J5^4 + 5*J2*J4*J5^4 + J2^10*J6 + J2*J3^6*J6 +
	    2*J3^4*J4^2*J6 + J2^2*J4^4*J6 + 6*J4^5*J6 + J3^5*J5*J6 + J2*J3^3*J4*J5*J6 + 2*J2^2*J3*J4^2*J5*J6 +
	    3*J2^2*J3^2*J5^2*J6 + 4*J3^2*J4*J5^2*J6 + J2*J4^2*J5^2*J6 + J2*J3*J5^3*J6 + 5*J5^4*J6 + 3*J2*J3^4*J6^2 +
	    J2^3*J4^2*J6^2 + 3*J3^2*J4^2*J6^2 + 3*J2*J4^3*J6^2 + 5*J3^3*J5*J6^2 + 3*J2*J3*J4*J5*J6^2 + J2^2*J5^2*J6^2
	    + 5*J4*J5^2*J6^2 + 6*J2*J3^2*J6^3 + 6*J2^2*J4*J6^3 + 4*J4^2*J6^3 + 6*J3*J5*J6^3 + 5*J2*J6^4 + J3^5*J4*J7 +
	    6*J2*J3^3*J4^2*J7 + 5*J2^2*J3*J4^3*J7 + 5*J2*J3^4*J5*J7 + 6*J3^2*J4^2*J5*J7 + 2*J2*J4^3*J5*J7 +
	    2*J3^3*J5^2*J7 + J2*J3*J4*J5^2*J7 + 5*J2^2*J5^3*J7 + J4*J5^3*J7 + 5*J3^3*J4*J6*J7 + 4*J2*J3*J4^2*J6*J7 +
	    5*J2*J3^2*J5*J6*J7 + 5*J2^2*J4*J5*J6*J7 + J4^2*J5*J6*J7 + 4*J3*J5^2*J6*J7 + 6*J2^2*J3*J6^2*J7 +
	    6*J3*J4*J6^2*J7 + 6*J2*J5*J6^2*J7 + 2*J3^4*J7^2 + 4*J2*J3^2*J4*J7^2 + 6*J2^2*J4^2*J7^2 + 3*J2^2*J3*J5*J7^2
	    + 6*J3*J4*J5*J7^2 + 6*J2*J5^2*J7^2 + 5*J3^2*J6*J7^2 + J2*J4*J6*J7^2 + 3*J6^2*J7^2 + 4*J2*J3*J7^3 +
	    5*J5*J7^3 + 2*J3^6*J8 + 5*J2*J3^4*J4*J8 + 5*J2^3*J4^3*J8 + 3*J2*J4^4*J8 + 4*J3^3*J4*J5*J8 +
	    2*J2*J3*J4^2*J5*J8 + 5*J2*J3^2*J5^2*J8 + 4*J2^2*J4*J5^2*J8 + 2*J4^2*J5^2*J8 + 5*J3*J5^3*J8 + 6*J3^4*J6*J8
	    + 4*J2*J3^2*J4*J6*J8 + 2*J2^2*J4^2*J6*J8 + 3*J4^3*J6*J8 + 4*J3*J4*J5*J6*J8 + 4*J2*J5^2*J6*J8 +
	    5*J2^3*J6^2*J8 + 4*J3^2*J6^2*J8 + J2*J4*J6^2*J8 + 2*J6^3*J8 + 2*J2*J3^3*J7*J8 + 6*J2^2*J3*J4*J7*J8 +
	    6*J3^2*J5*J7*J8 + 6*J2*J4*J5*J7*J8 + 6*J2*J3*J6*J7*J8 + J5*J6*J7*J8 + 5*J2^2*J7^2*J8 + 2*J4*J7^2*J8 +
	    4*J2^3*J4*J8^2 + 4*J3^2*J4*J8^2 + 6*J2*J4^2*J8^2 + 6*J2*J3*J5*J8^2 + 3*J2^2*J6*J8^2 + 3*J4*J6*J8^2 +
	    4*J3*J7*J8^2 + 2*J2*J8^3 + 2*J2*J3^5*J9 + 4*J3^3*J4^2*J9 + 4*J2*J3*J4^3*J9 + J3^4*J5*J9 +
	    3*J2*J3^2*J4*J5*J9 + 6*J2^2*J4^2*J5*J9 + 5*J4^3*J5*J9 + 3*J2^2*J3*J5^2*J9 + 6*J3*J4*J5^2*J9 + 2*J2*J5^3*J9
	    + 6*J2*J3^3*J6*J9 + 3*J2^2*J3*J4*J6*J9 + 3*J3*J4^2*J6*J9 + 6*J3^2*J5*J6*J9 + 6*J2*J4*J5*J6*J9 +
	    4*J2*J3*J6^2*J9 + 4*J5*J6^2*J9 + 5*J2^2*J3^2*J7*J9 + 5*J3^2*J4*J7*J9 + 5*J5^2*J7*J9 + J2^2*J6*J7*J9 +
	    6*J4*J6*J7*J9 + 4*J3*J7^2*J9 + 5*J3^3*J8*J9 + J2*J3*J4*J8*J9 + 2*J2^2*J5*J8*J9 + 2*J4*J5*J8*J9 +
	    3*J2*J7*J8*J9 + 4*J2^2*J4*J9^2 + 3*J4^2*J9^2 + J2*J6*J9^2 + 5*J8*J9^2 + 2*J3^4*J4*J10 + 6*J2*J3^2*J4^2*J10
	    + 4*J2^2*J4^3*J10 + 4*J2^2*J3*J4*J5*J10 + 4*J3*J4^2*J5*J10 + 5*J3^2*J5^2*J10 + 4*J2^3*J4*J6*J10 +
	    4*J3^2*J4*J6*J10 + 3*J2*J4^2*J6*J10 + 5*J2*J3*J5*J6*J10 + J5^2*J6*J10 + 5*J2^2*J6^2*J10 + J4*J6^2*J10 +
	    3*J3^3*J7*J10 + 2*J2^2*J5*J7*J10 + J4*J5*J7*J10 + 5*J3*J6*J7*J10 + 5*J2*J7^2*J10 + 4*J2^4*J8*J10 +
	    4*J2*J3^2*J8*J10 + 4*J4^2*J8*J10 + J3*J5*J8*J10 + 4*J2*J6*J8*J10 + 5*J8^2*J10 + 5*J2^2*J3*J9*J10 +
	    4*J3*J4*J9*J10 + 5*J2*J5*J9*J10 + J7*J9*J10 + 3*J3^2*J10^2 + 4*J2*J4*J10^2 + 4*J6*J10^2 + 2*J3^5*J11 +
	    5*J2*J3^3*J4*J11 + 3*J2^2*J3*J4^2*J11 + 6*J3*J4^3*J11 + 5*J3^2*J4*J5*J11 + 4*J2*J3*J5^2*J11 + 3*J5^3*J11 +
	    6*J3^3*J6*J11 + J2*J3*J4*J6*J11 + 3*J2^2*J5*J6*J11 + 5*J4*J5*J6*J11 + 3*J3*J6^2*J11 + 3*J2*J3^2*J7*J11 +
	    2*J2^2*J4*J7*J11 + 6*J4^2*J7*J11 + 3*J3*J5*J7*J11 + 4*J2*J6*J7*J11 + 6*J2^2*J3*J8*J11 + 3*J3*J4*J8*J11 +
	    J2*J5*J8*J11 + 2*J7*J8*J11 + 4*J3^2*J9*J11 + 4*J2*J4*J9*J11 + 5*J6*J9*J11 + 4*J2*J3*J10*J11 + 4*J5*J10*J11
	    + 5*J4*J11^2 + 5*J3^3*J4*J13 + 6*J2*J3*J4^2*J13 + 3*J2*J3^2*J5*J13 + J2^2*J4*J5*J13 + 5*J4^2*J5*J13 +
	    6*J2^2*J3*J6*J13 + 6*J3*J4*J6*J13 + J2*J5*J6*J13 + J3^2*J7*J13 + 2*J2*J4*J7*J13 + J6*J7*J13 +
	    3*J2*J3*J8*J13 + 6*J5*J8*J13 + 4*J4*J9*J13 + 2*J3*J10*J13 + 5*J2*J11*J13 + 2*J13^2 + J3^4*J14 +
	    2*J2*J3^2*J4*J14 + 4*J2^2*J4^2*J14 + J4^3*J14 + 5*J3*J4*J5*J14 + 2*J2*J5^2*J14 + 6*J2^3*J6*J14 +
	    4*J3^2*J6*J14 + 6*J2*J4*J6*J14 + 3*J2*J3*J7*J14 + 5*J5*J7*J14 + 6*J2*J3^3*J15 + 2*J2^2*J3*J4*J15 +
	    4*J3*J4^2*J15 + J3^2*J5*J15 + 2*J2*J4*J5*J15 + 3*J2*J3*J6*J15 + 2*J5*J6*J15 + 4*J2^2*J7*J15 + 4*J4*J7*J15
	    + J3*J8*J15 + 3*J2*J9*J15 + 4*J11*J15);

	// deg. 28
	Append(~RES,
	    J2^14 + 4*J3^8*J4 + 4*J2*J3^6*J4^2 + 4*J2*J3^2*J4^5 + 6*J2^2*J4^6 + 6*J4^7 + 3*J3^5*J4^2*J5 +
	    6*J2*J3^3*J4^3*J5 + 3*J3*J4^5*J5 + J3^6*J5^2 + J2*J3^4*J4*J5^2 + 4*J3^2*J4^3*J5^2 + J2*J4^4*J5^2 +
	    6*J3^3*J4*J5^3 + 5*J2*J3*J4^2*J5^3 + 5*J2*J3^2*J5^4 + J2^2*J4*J5^4 + 4*J3^6*J4*J6 + J2*J3^4*J4^2*J6 +
	    4*J3^2*J4^4*J6 + 3*J2*J4^5*J6 + 2*J3^3*J4^2*J5*J6 + J2*J3*J4^3*J5*J6 + 3*J3^4*J5^2*J6 +
	    6*J2*J3^2*J4*J5^2*J6 + J2^2*J4^2*J5^2*J6 + 4*J2^2*J3*J5^3*J6 + J3*J4*J5^3*J6 + 6*J2*J5^4*J6 +
	    6*J3^4*J4*J6^2 + 4*J2*J3^2*J4^2*J6^2 + 2*J2^2*J4^3*J6^2 + 6*J4^4*J6^2 + 4*J2*J3^3*J5*J6^2 +
	    3*J2^2*J3*J4*J5*J6^2 + 2*J3*J4^2*J5*J6^2 + 2*J3^2*J5^2*J6^2 + 2*J2*J4*J5^2*J6^2 + 3*J2^3*J4*J6^3 +
	    3*J3^2*J4*J6^3 + 6*J2*J4^2*J6^3 + J2*J3*J5*J6^3 + 5*J5^2*J6^3 + J2^2*J6^4 + 2*J4*J6^4 + 6*J3^7*J7 +
	    3*J2*J3^5*J4*J7 + 3*J2*J3*J4^4*J7 + 2*J3^4*J4*J5*J7 + 2*J2*J3^2*J4^2*J5*J7 + 6*J2^2*J4^3*J5*J7 +
	    3*J4^4*J5*J7 + 6*J2*J3^3*J5^2*J7 + 2*J2^2*J3*J4*J5^2*J7 + 3*J3*J4^2*J5^2*J7 + 5*J3^2*J5^3*J7 +
	    4*J2*J4*J5^3*J7 + J3^5*J6*J7 + 3*J2*J3^3*J4*J6*J7 + 6*J2^2*J3*J4^2*J6*J7 + 3*J3*J4^3*J6*J7 +
	    5*J3^2*J4*J5*J6*J7 + 4*J2*J4^2*J5*J6*J7 + 5*J2*J3*J5^2*J6*J7 + 5*J5^3*J6*J7 + 6*J3^3*J6^2*J7 +
	    J2*J3*J4*J6^2*J7 + 2*J2^2*J5*J6^2*J7 + 3*J4*J5*J6^2*J7 + 3*J3*J6^3*J7 + 4*J2*J3^4*J7^2 + 5*J3^2*J4^2*J7^2
	    + 3*J2*J4^3*J7^2 + 2*J3^3*J5*J7^2 + 3*J2*J3*J4*J5*J7^2 + 3*J2^2*J5^2*J7^2 + J4*J5^2*J7^2 +
	    3*J2*J3^2*J6*J7^2 + 6*J2^2*J4*J6*J7^2 + 5*J4^2*J6*J7^2 + 6*J3*J5*J6*J7^2 + 2*J2*J6^2*J7^2 + 2*J2^2*J3*J7^3
	    + 6*J3*J4*J7^3 + J2*J5*J7^3 + 4*J7^4 + 5*J2*J3^6*J8 + 5*J3^4*J4^2*J8 + 6*J2*J3^2*J4^3*J8 + 5*J2^2*J4^4*J8
	    + 2*J4^5*J8 + 5*J3^5*J5*J8 + 4*J2*J3^3*J4*J5*J8 + 5*J3^2*J4*J5^2*J8 + 2*J2*J4^2*J5^2*J8 + 5*J2*J3*J5^3*J8
	    + 3*J2^3*J4^2*J6*J8 + 4*J3^2*J4^2*J6*J8 + 4*J2*J4^3*J6*J8 + 4*J3^3*J5*J6*J8 + 3*J2*J3*J4*J5*J6*J8 +
	    6*J4*J5^2*J6*J8 + 5*J2*J3^2*J6^2*J8 + 3*J2^2*J4*J6^2*J8 + 6*J4^2*J6^2*J8 + 4*J3*J5*J6^2*J8 + 3*J2*J6^3*J8
	    + 3*J3^3*J4*J7*J8 + 4*J2*J3*J4^2*J7*J8 + 3*J2*J3^2*J5*J7*J8 + 6*J2^2*J4*J5*J7*J8 + 6*J4^2*J5*J7*J8 +
	    J2^2*J3*J6*J7*J8 + 4*J3*J4*J6*J7*J8 + 2*J2*J5*J6*J7*J8 + 3*J3^2*J7^2*J8 + 4*J2*J4*J7^2*J8 + 6*J6*J7^2*J8 +
	    2*J3^4*J8^2 + 2*J2*J3^2*J4*J8^2 + 5*J4^3*J8^2 + 2*J2^2*J3*J5*J8^2 + 5*J3*J4*J5*J8^2 + J2*J5^2*J8^2 +
	    4*J3^2*J6*J8^2 + 3*J6^2*J8^2 + 6*J2*J3*J7*J8^2 + J4*J8^3 + 3*J3^5*J4*J9 + 2*J2*J3^3*J4^2*J9 + 6*J3*J4^4*J9
	    + 3*J2*J3^4*J5*J9 + 6*J3^2*J4^2*J5*J9 + 3*J3^3*J5^2*J9 + 3*J2^2*J5^3*J9 + 4*J4*J5^3*J9 + J3^3*J4*J6*J9 +
	    5*J2*J3*J4^2*J6*J9 + 6*J2^2*J4*J5*J6*J9 + 2*J4^2*J5*J6*J9 + 5*J3*J5^2*J6*J9 + 3*J2^2*J3*J6^2*J9 +
	    2*J3^4*J7*J9 + 5*J2*J3^2*J4*J7*J9 + J2^2*J4^2*J7*J9 + 6*J4^3*J7*J9 + 4*J2^2*J3*J5*J7*J9 + 5*J3*J4*J5*J7*J9
	    + 4*J3^2*J6*J7*J9 + 5*J2*J4*J6*J7*J9 + 4*J6^2*J7*J9 + 6*J2*J3*J7^2*J9 + 3*J5*J7^2*J9 + 3*J2*J4*J5*J8*J9 +
	    3*J2*J3*J6*J8*J9 + 4*J5*J6*J8*J9 + 3*J2^2*J7*J8*J9 + 4*J4*J7*J8*J9 + 4*J3^2*J4*J9^2 + 5*J2*J4^2*J9^2 +
	    5*J5^2*J9^2 + 4*J2^2*J6*J9^2 + 4*J4*J6*J9^2 + 2*J3*J7*J9^2 + J3^6*J10 + 3*J2*J3^4*J4*J10 + 6*J3^2*J4^3*J10
	    + 5*J2*J4^4*J10 + 4*J3^3*J4*J5*J10 + 6*J2*J3*J4^2*J5*J10 + 4*J2*J3^2*J5^2*J10 + 3*J2^2*J4*J5^2*J10 +
	    J4^2*J5^2*J10 + 5*J2*J3^2*J4*J6*J10 + J2^2*J4^2*J6*J10 + 3*J4^3*J6*J10 + 4*J2^2*J3*J5*J6*J10 +
	    6*J2^3*J6^2*J10 + 3*J3^2*J6^2*J10 + 5*J2*J4*J6^2*J10 + 3*J6^3*J10 + 2*J2*J3^3*J7*J10 + 6*J2^2*J3*J4*J7*J10
	    + J3^2*J5*J7*J10 + 4*J2*J4*J5*J7*J10 + 5*J2*J3*J6*J7*J10 + 4*J5*J6*J7*J10 + 2*J2^2*J7^2*J10 +
	    6*J4*J7^2*J10 + 3*J2^3*J4*J8*J10 + 2*J3^2*J4*J8*J10 + 4*J5^2*J8*J10 + J2^2*J6*J8*J10 + 2*J3*J7*J8*J10 +
	    3*J2*J8^2*J10 + 5*J3^3*J9*J10 + 2*J2*J3*J4*J9*J10 + 6*J2^2*J5*J9*J10 + J4*J5*J9*J10 + 4*J3*J6*J9*J10 +
	    3*J2*J7*J9*J10 + 5*J9^2*J10 + 4*J2^2*J4*J10^2 + J4^2*J10^2 + 6*J3*J5*J10^2 + 3*J2*J6*J10^2 + 3*J8*J10^2 +
	    J3^3*J4^2*J11 + J2*J3*J4^3*J11 + 2*J3^4*J5*J11 + 3*J2*J3^2*J4*J5*J11 + 5*J2^2*J4^2*J5*J11 + 3*J4^3*J5*J11
	    + 5*J2^2*J3*J5^2*J11 + 6*J3*J4*J5^2*J11 + 6*J2*J3^3*J6*J11 + 4*J2^2*J3*J4*J6*J11 + 3*J3*J4^2*J6*J11 +
	    2*J3^2*J5*J6*J11 + 5*J2*J3*J6^2*J11 + 5*J5*J6^2*J11 + 6*J3^2*J4*J7*J11 + 6*J2*J4^2*J7*J11 +
	    J2*J3*J5*J7*J11 + 4*J2^2*J6*J7*J11 + 3*J4*J6*J7*J11 + 5*J3^3*J8*J11 + 4*J2*J3*J4*J8*J11 + 3*J2^2*J5*J8*J11
	    + 2*J3*J6*J8*J11 + 5*J2*J7*J8*J11 + 5*J2*J3^2*J9*J11 + 4*J4^2*J9*J11 + J3*J5*J9*J11 + 2*J2*J6*J9*J11 +
	    4*J8*J9*J11 + 5*J2^2*J3*J10*J11 + 3*J3*J4*J10*J11 + 2*J2*J5*J10*J11 + 2*J7*J10*J11 + 5*J3^2*J11^2 +
	    5*J2*J4*J11^2 + J6*J11^2 + 4*J3^5*J13 + J2*J3^3*J4*J13 + 5*J2^2*J3*J4^2*J13 + J3*J4^3*J13 +
	    5*J3^2*J4*J5*J13 + 2*J2*J4^2*J5*J13 + 2*J2*J3*J5^2*J13 + 4*J5^3*J13 + J3^3*J6*J13 + 4*J2^2*J5*J6*J13 +
	    2*J3*J6^2*J13 + 3*J2*J3^2*J7*J13 + 2*J2^2*J4*J7*J13 + 4*J4^2*J7*J13 + 5*J3*J5*J7*J13 + 5*J2*J6*J7*J13 +
	    6*J2^2*J3*J8*J13 + 4*J3*J4*J8*J13 + 4*J2*J5*J8*J13 + J7*J8*J13 + 6*J3^2*J9*J13 + 6*J2*J4*J9*J13 +
	    3*J6*J9*J13 + 2*J2*J3*J10*J13 + 6*J5*J10*J13 + 4*J2^2*J11*J13 + J4*J11*J13 + 5*J2*J13^2 + 5*J2*J3^4*J14 +
	    5*J3^2*J4^2*J14 + J2*J4^3*J14 + 6*J3^3*J5*J14 + 3*J2*J3*J4*J5*J14 + 5*J2^2*J5^2*J14 + 3*J4*J5^2*J14 +
	    5*J4^2*J6*J14 + 4*J3*J5*J6*J14 + 4*J2*J6^2*J14 + 6*J2^2*J3*J7*J14 + 6*J3*J4*J7*J14 + 6*J2*J5*J7*J14 +
	    J7^2*J14 + 5*J2^3*J8*J14 + 2*J2*J4*J8*J14 + 6*J2*J3*J9*J14 + 5*J5*J9*J14 + 3*J2^2*J10*J14 + 5*J3*J11*J14 +
	    5*J3^3*J4*J15 + 3*J2*J3*J4^2*J15 + 4*J2*J3^2*J5*J15 + 3*J2^2*J4*J5*J15 + 5*J4^2*J5*J15 + 2*J3*J5^2*J15 +
	    4*J3*J4*J6*J15 + 2*J2*J5*J6*J15 + 3*J3^2*J7*J15 + 5*J2*J4*J7*J15 + 2*J6*J7*J15 + J5*J8*J15 + 2*J2^2*J9*J15
	    + 4*J4*J9*J15 + 5*J3*J10*J15 + 5*J2*J11*J15 + J13*J15);

	// deg. 30
	Append(~RES,
	    J2^15 + 6*J3^10 + 2*J2*J3^8*J4 + 2*J2*J3^4*J4^4 + 2*J3^2*J4^6 + 4*J2*J4^7 + 4*J3^7*J4*J5 + 2*J2*J3^5*J4^2*J5 +
	    J2*J3*J4^5*J5 + 6*J3^4*J4^2*J5^2 + J2*J3^2*J4^3*J5^2 + 4*J2^2*J4^4*J5^2 + 6*J3^5*J5^3 + 2*J2*J3^3*J4*J5^3
	    + 6*J3*J4^3*J5^3 + 4*J3^2*J4*J5^4 + J2*J4^2*J5^4 + 2*J3^8*J6 + 5*J2*J3^2*J4^4*J6 + 2*J2^2*J4^5*J6 +
	    4*J4^6*J6 + J3^5*J4*J5*J6 + 6*J2*J3^3*J4^2*J5*J6 + 2*J3*J4^4*J5*J6 + 3*J2*J3^4*J5^2*J6 + J3^2*J4^2*J5^2*J6
	    + 6*J2*J4^3*J5^2*J6 + J3^3*J5^3*J6 + 5*J2*J3*J4*J5^3*J6 + 4*J4*J5^4*J6 + 4*J3^2*J4^3*J6^2 + 3*J2*J4^4*J6^2
	    + 6*J3^3*J4*J5*J6^2 + 5*J2*J3^2*J5^2*J6^2 + J2^2*J4*J5^2*J6^2 + 6*J4^2*J5^2*J6^2 + 2*J3*J5^3*J6^2 +
	    3*J3^4*J6^3 + J2*J3^2*J4*J6^3 + 3*J2^2*J4^2*J6^3 + J3*J4*J5*J6^3 + 6*J2^3*J6^4 + J3^2*J6^4 + 4*J2*J3^7*J7
	    + 4*J3^5*J4^2*J7 + J3*J4^5*J7 + 2*J3^6*J5*J7 + 6*J3^2*J4^3*J5*J7 + 3*J2*J4^4*J5*J7 + 4*J3^3*J4*J5^2*J7 +
	    2*J2^2*J4*J5^3*J7 + J4^2*J5^3*J7 + 3*J3*J5^4*J7 + 5*J2*J3^5*J6*J7 + 6*J3^3*J4^2*J6*J7 + 4*J2*J3*J4^3*J6*J7
	    + 5*J3^4*J5*J6*J7 + 3*J2*J3^2*J4*J5*J6*J7 + 3*J2^2*J4^2*J5*J6*J7 + 6*J4^3*J5*J6*J7 + 5*J3*J4*J5^2*J6*J7 +
	    5*J2*J5^3*J6*J7 + 2*J2*J3^3*J6^2*J7 + 5*J3*J4^2*J6^2*J7 + 6*J3^2*J5*J6^2*J7 + 2*J2*J4*J5*J6^2*J7 +
	    6*J5*J6^3*J7 + 6*J3^4*J4*J7^2 + 6*J2*J3^2*J4^2*J7^2 + 5*J2^2*J4^3*J7^2 + 3*J4^4*J7^2 + 2*J2*J3^3*J5*J7^2 +
	    2*J3^2*J5^2*J7^2 + 3*J2*J4*J5^2*J7^2 + 4*J2*J4^2*J6*J7^2 + 4*J2*J3*J5*J6*J7^2 + 5*J2^2*J6^2*J7^2 +
	    4*J4*J6^2*J7^2 + 2*J3^3*J7^3 + 5*J2*J3*J4*J7^3 + 2*J2^2*J5*J7^3 + 5*J4*J5*J7^3 + 4*J3*J6*J7^3 + 4*J2*J7^4
	    + J2*J3^4*J4^2*J8 + 5*J2*J4^5*J8 + 5*J2*J3^5*J5*J8 + 3*J3^3*J4^2*J5*J8 + 3*J2*J3*J4^3*J5*J8 + J3^4*J5^2*J8
	    + 5*J2*J3^2*J4*J5^2*J8 + 5*J2^2*J4^2*J5^2*J8 + 2*J4^3*J5^2*J8 + 6*J3*J4*J5^3*J8 + 4*J2*J5^4*J8 +
	    3*J3^4*J4*J6*J8 + 4*J2^2*J4^3*J6*J8 + J4^4*J6*J8 + 2*J2*J3^3*J5*J6*J8 + 2*J3*J4^2*J5*J6*J8 +
	    J3^2*J5^2*J6*J8 + 6*J2*J4*J5^2*J6*J8 + 6*J3^2*J4*J6^2*J8 + 3*J2*J4^2*J6^2*J8 + J2*J3*J5*J6^2*J8 +
	    6*J5^2*J6^2*J8 + 4*J2^2*J6^3*J8 + 2*J4*J6^3*J8 + J3^5*J7*J8 + 5*J2*J3^3*J4*J7*J8 + 4*J3*J4^3*J7*J8 +
	    J3^2*J4*J5*J7*J8 + J2*J4^2*J5*J7*J8 + 4*J2*J3*J5^2*J7*J8 + J5^3*J7*J8 + 3*J3^3*J6*J7*J8 +
	    3*J2*J3*J4*J6*J7*J8 + 3*J2^2*J5*J6*J7*J8 + 5*J4*J5*J6*J7*J8 + 6*J3*J6^2*J7*J8 + 6*J2*J3^2*J7^2*J8 +
	    3*J2^2*J4*J7^2*J8 + 6*J4^2*J7^2*J8 + 4*J2*J6*J7^2*J8 + 3*J2*J3^4*J8^2 + 6*J2^3*J4^2*J8^2 + J3^2*J4^2*J8^2
	    + 3*J2*J4^3*J8^2 + J3^3*J5*J8^2 + 6*J2^2*J5^2*J8^2 + 5*J4*J5^2*J8^2 + 6*J2^2*J4*J6*J8^2 + 3*J3*J5*J6*J8^2
	    + 2*J2*J6^2*J8^2 + 5*J3*J4*J7*J8^2 + 4*J2*J5*J7*J8^2 + 4*J2^3*J8^3 + 4*J6*J8^3 + J3^7*J9 + 4*J2*J3^5*J4*J9
	    + 3*J3^3*J4^3*J9 + 6*J2*J3*J4^4*J9 + 5*J3^4*J4*J5*J9 + 2*J2*J3^2*J4^2*J5*J9 + J2^2*J4^3*J5*J9 +
	    3*J4^4*J5*J9 + 4*J2*J3^3*J5^2*J9 + J3*J4^2*J5^2*J9 + J3^2*J5^3*J9 + J2*J4*J5^3*J9 + 3*J3^5*J6*J9 +
	    3*J2*J3^3*J4*J6*J9 + 6*J3*J4^3*J6*J9 + 3*J2*J4^2*J5*J6*J9 + 5*J2*J3*J5^2*J6*J9 + 6*J3^3*J6^2*J9 +
	    4*J2^2*J5*J6^2*J9 + 2*J3*J6^3*J9 + 5*J2*J3^4*J7*J9 + 6*J3^2*J4^2*J7*J9 + 2*J2*J4^3*J7*J9 + 5*J3^3*J5*J7*J9
	    + 6*J2*J3*J4*J5*J7*J9 + 3*J2^2*J5^2*J7*J9 + 5*J4*J5^2*J7*J9 + 3*J2*J3^2*J6*J7*J9 + 2*J2^2*J4*J6*J7*J9 +
	    J2*J6^2*J7*J9 + 4*J2^2*J3*J7^2*J9 + 3*J3*J4*J7^2*J9 + 5*J2*J5*J7^2*J9 + J7^3*J9 + 5*J3^3*J4*J8*J9 +
	    3*J2*J3*J4^2*J8*J9 + 2*J2*J3^2*J5*J8*J9 + 6*J2^2*J4*J5*J8*J9 + 4*J4^2*J5*J8*J9 + 3*J3*J5^2*J8*J9 +
	    5*J3*J4*J6*J8*J9 + 4*J2*J5*J6*J8*J9 + 4*J3^2*J7*J8*J9 + 2*J2*J4*J7*J8*J9 + J6*J7*J8*J9 + 3*J2*J3*J8^2*J9 +
	    J5*J8^2*J9 + 4*J3^4*J9^2 + 5*J2*J3^2*J4*J9^2 + 5*J2^2*J4^2*J9^2 + 6*J4^3*J9^2 + 2*J3*J4*J5*J9^2 +
	    3*J2*J5^2*J9^2 + 6*J2*J3*J7*J9^2 + 3*J5*J7*J9^2 + 3*J2^2*J8*J9^2 + 4*J3*J9^3 + 2*J2*J3^6*J10 +
	    3*J3^4*J4^2*J10 + 2*J2*J3^2*J4^3*J10 + 6*J2^2*J4^4*J10 + J4^5*J10 + J3^5*J5*J10 + 6*J2*J3^3*J4*J5*J10 +
	    6*J3*J4^3*J5*J10 + 3*J3^2*J4*J5^2*J10 + J2*J4^2*J5^2*J10 + 5*J2*J3*J5^3*J10 + 6*J5^4*J10 +
	    6*J3^2*J4^2*J6*J10 + 2*J2*J4^3*J6*J10 + 6*J3^3*J5*J6*J10 + 5*J2*J3*J4*J5*J6*J10 + 3*J2^2*J5^2*J6*J10 +
	    4*J4*J5^2*J6*J10 + 3*J2*J3^2*J6^2*J10 + 6*J2^2*J4*J6^2*J10 + 4*J4^2*J6^2*J10 + 3*J3*J5*J6^2*J10 +
	    4*J2*J6^3*J10 + 5*J3^3*J4*J7*J10 + 2*J2*J3*J4^2*J7*J10 + 6*J2*J3^2*J5*J7*J10 + 3*J2^2*J4*J5*J7*J10 +
	    6*J4^2*J5*J7*J10 + 6*J3*J5^2*J7*J10 + 3*J2^2*J3*J6*J7*J10 + 4*J3*J4*J6*J7*J10 + 6*J2*J5*J6*J7*J10 +
	    J3^2*J7^2*J10 + 3*J2*J4*J7^2*J10 + 3*J6*J7^2*J10 + 4*J2*J3^2*J4*J8*J10 + 4*J2^2*J4^2*J8*J10 +
	    2*J4^3*J8*J10 + 5*J3*J4*J5*J8*J10 + 2*J2*J5^2*J8*J10 + 3*J2^3*J6*J8*J10 + 4*J3^2*J6*J8*J10 +
	    2*J2*J4*J6*J8*J10 + 5*J6^2*J8*J10 + 3*J2*J3*J7*J8*J10 + J5*J7*J8*J10 + 2*J2^2*J8^2*J10 + 2*J4*J8^2*J10 +
	    2*J2*J3^3*J9*J10 + 6*J3*J4^2*J9*J10 + 3*J3^2*J5*J9*J10 + 4*J2*J4*J5*J9*J10 + 3*J5*J6*J9*J10 +
	    J2^2*J7*J9*J10 + 6*J3*J8*J9*J10 + 4*J2*J9^2*J10 + 6*J2^3*J4*J10^2 + 4*J3^2*J4*J10^2 + 4*J2*J4^2*J10^2 +
	    6*J2*J3*J5*J10^2 + 3*J5^2*J10^2 + 2*J2^2*J6*J10^2 + 2*J4*J6*J10^2 + 3*J3*J7*J10^2 + 3*J2*J8*J10^2 +
	    4*J10^3 + 2*J3*J4^4*J11 + 2*J2*J3^4*J5*J11 + J2*J4^3*J5*J11 + 4*J2*J3*J4*J5^2*J11 + 5*J2^2*J5^3*J11 +
	    3*J4*J5^3*J11 + J3^3*J4*J6*J11 + 3*J2*J3*J4^2*J6*J11 + 6*J2*J3^2*J5*J6*J11 + J2^2*J4*J5*J6*J11 +
	    2*J3*J5^2*J6*J11 + 2*J2^2*J3*J6^2*J11 + 6*J3*J4*J6^2*J11 + 5*J2*J5*J6^2*J11 + 2*J3^4*J7*J11 +
	    6*J2*J3^2*J4*J7*J11 + 2*J2^2*J4^2*J7*J11 + 6*J4^3*J7*J11 + 4*J3*J4*J5*J7*J11 + 5*J2*J5^2*J7*J11 +
	    2*J3^2*J6*J7*J11 + 5*J2*J4*J6*J7*J11 + 2*J5*J7^2*J11 + 3*J2*J3^3*J8*J11 + J3*J4^2*J8*J11 +
	    3*J3^2*J5*J8*J11 + 3*J2*J4*J5*J8*J11 + 3*J2*J3*J6*J8*J11 + 2*J2^2*J7*J8*J11 + 5*J4*J7*J8*J11 +
	    3*J3*J8^2*J11 + 2*J3^2*J4*J9*J11 + J2*J4^2*J9*J11 + 6*J2*J3*J5*J9*J11 + 3*J5^2*J9*J11 + 2*J2^2*J6*J9*J11 +
	    3*J4*J6*J9*J11 + J3*J7*J9*J11 + 4*J2*J8*J9*J11 + J3^3*J10*J11 + J2*J3*J4*J10*J11 + 2*J3*J6*J10*J11 +
	    5*J9*J10*J11 + J2^2*J4*J11^2 + 3*J4^2*J11^2 + 2*J2*J6*J11^2 + 5*J8*J11^2 + 3*J2*J3^5*J13 + 2*J3^3*J4^2*J13
	    + J2*J3*J4^3*J13 + 4*J3^4*J5*J13 + 4*J2*J3^2*J4*J5*J13 + 3*J4^3*J5*J13 + 5*J2*J5^3*J13 + 2*J2*J3^3*J6*J13
	    + 2*J3*J4^2*J6*J13 + 3*J3^2*J5*J6*J13 + 4*J2*J4*J5*J6*J13 + 4*J2*J3*J6^2*J13 + 4*J5*J6^2*J13 +
	    5*J3^2*J4*J7*J13 + 4*J2*J4^2*J7*J13 + 6*J2*J3*J5*J7*J13 + 5*J5^2*J7*J13 + 5*J2^2*J6*J7*J13 +
	    5*J4*J6*J7*J13 + J3^3*J8*J13 + 6*J2^2*J5*J8*J13 + 3*J3*J6*J8*J13 + J2*J7*J8*J13 + J2*J3^2*J9*J13 +
	    2*J2^2*J4*J9*J13 + 2*J4^2*J9*J13 + 5*J3*J5*J9*J13 + 3*J8*J9*J13 + 5*J2^2*J3*J10*J13 + 3*J3*J4*J10*J13 +
	    4*J2*J5*J10*J13 + 3*J7*J10*J13 + 6*J6*J11*J13 + 5*J4*J13^2 + 4*J3^4*J4*J14 + 2*J2*J3^2*J4^2*J14 +
	    5*J4^4*J14 + 6*J2*J3^3*J5*J14 + 3*J3*J4^2*J5*J14 + 4*J2*J4*J5^2*J14 + 4*J2^3*J4*J6*J14 + 2*J3^2*J4*J6*J14
	    + 6*J2*J4^2*J6*J14 + 4*J2*J3*J5*J6*J14 + 4*J5^2*J6*J14 + 3*J2^2*J6^2*J14 + 5*J4*J6^2*J14 + J3^3*J7*J14 +
	    2*J2*J3*J4*J7*J14 + J2^2*J5*J7*J14 + 2*J4*J5*J7*J14 + 2*J3*J6*J7*J14 + 3*J2*J7^2*J14 + 5*J2*J3^2*J8*J14 +
	    4*J2^2*J4*J8*J14 + J4^2*J8*J14 + 5*J3*J5*J8*J14 + 4*J2*J6*J8*J14 + 6*J2^2*J3*J9*J14 + 4*J3*J4*J9*J14 +
	    3*J2*J5*J9*J14 + 2*J7*J9*J14 + 3*J2^3*J10*J14 + 3*J3^2*J10*J14 + 5*J2*J4*J10*J14 + 2*J2*J3*J11*J14 +
	    5*J3*J13*J14 + 2*J3^5*J15 + 5*J2*J3^3*J4*J15 + 3*J3*J4^3*J15 + 5*J3^2*J4*J5*J15 + J2*J4^2*J5*J15 +
	    5*J2*J3*J5^2*J15 + 6*J5^3*J15 + 6*J3^3*J6*J15 + 6*J2^2*J5*J6*J15 + J3*J6^2*J15 + 3*J2*J3^2*J7*J15 +
	    5*J2^2*J4*J7*J15 + 5*J4^2*J7*J15 + 3*J3*J5*J7*J15 + 4*J2*J6*J7*J15 + J2^2*J3*J8*J15 + J3*J4*J8*J15 +
	    2*J7*J8*J15 + 3*J3^2*J9*J15 + 5*J2*J4*J9*J15 + 3*J6*J9*J15 + 6*J2*J3*J10*J15 + 5*J5*J10*J15 + J2^2*J11*J15
	    + 4*J4*J11*J15 + 6*J2*J13*J15 + 4*J15^2);

	if ratsolve eq false then return RES; end if;

	V := Variety(ideal<P7 | RES>);
	for v in V do
	    NJI := [ FF!v[1], J3, J4, J5, J6, FF!v[2], FF!v[3], FF!v[4], J10, FF!v[5], FF!v[6], J14, FF!v[7]];
	    NJI := WPSNormalize([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15], NJI);
	    JIs := ShiodaInvariantsAppend(JIs, NJI);
	end for;
    end for;

    return JIs;

end function;
