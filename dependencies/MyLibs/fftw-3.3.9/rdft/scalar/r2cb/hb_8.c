/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Thu Dec 10 07:06:27 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hb_8 -include rdft/scalar/hb.h */

/*
 * This function contains 66 FP additions, 36 FP multiplications,
 * (or, 44 additions, 14 multiplications, 22 fused multiply/add),
 * 33 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hb.h"

static void hb_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 14, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E T7, T1i, T1n, Tk, TD, TV, T1b, TQ, Te, T1e, T1o, T1j, TE, TF, TR;
	       E Tv, TW;
	       {
		    E T3, Tg, TC, T19, T6, Tz, Tj, T1a;
		    {
			 E T1, T2, TA, TB;
			 T1 = cr[0];
			 T2 = ci[WS(rs, 3)];
			 T3 = T1 + T2;
			 Tg = T1 - T2;
			 TA = ci[WS(rs, 7)];
			 TB = cr[WS(rs, 4)];
			 TC = TA + TB;
			 T19 = TA - TB;
		    }
		    {
			 E T4, T5, Th, Ti;
			 T4 = cr[WS(rs, 2)];
			 T5 = ci[WS(rs, 1)];
			 T6 = T4 + T5;
			 Tz = T4 - T5;
			 Th = ci[WS(rs, 5)];
			 Ti = cr[WS(rs, 6)];
			 Tj = Th + Ti;
			 T1a = Th - Ti;
		    }
		    T7 = T3 + T6;
		    T1i = T3 - T6;
		    T1n = T19 - T1a;
		    Tk = Tg - Tj;
		    TD = Tz + TC;
		    TV = TC - Tz;
		    T1b = T19 + T1a;
		    TQ = Tg + Tj;
	       }
	       {
		    E Ta, Tl, Tt, T1d, Td, Tq, To, T1c, Tp, Tu;
		    {
			 E T8, T9, Tr, Ts;
			 T8 = cr[WS(rs, 1)];
			 T9 = ci[WS(rs, 2)];
			 Ta = T8 + T9;
			 Tl = T8 - T9;
			 Tr = ci[WS(rs, 4)];
			 Ts = cr[WS(rs, 7)];
			 Tt = Tr + Ts;
			 T1d = Tr - Ts;
		    }
		    {
			 E Tb, Tc, Tm, Tn;
			 Tb = ci[0];
			 Tc = cr[WS(rs, 3)];
			 Td = Tb + Tc;
			 Tq = Tb - Tc;
			 Tm = ci[WS(rs, 6)];
			 Tn = cr[WS(rs, 5)];
			 To = Tm + Tn;
			 T1c = Tm - Tn;
		    }
		    Te = Ta + Td;
		    T1e = T1c + T1d;
		    T1o = Ta - Td;
		    T1j = T1d - T1c;
		    TE = Tl + To;
		    TF = Tq + Tt;
		    TR = TE + TF;
		    Tp = Tl - To;
		    Tu = Tq - Tt;
		    Tv = Tp + Tu;
		    TW = Tp - Tu;
	       }
	       cr[0] = T7 + Te;
	       ci[0] = T1b + T1e;
	       {
		    E TS, TX, TT, TY, TP, TU;
		    TS = FNMS(KP707106781, TR, TQ);
		    TX = FMA(KP707106781, TW, TV);
		    TP = W[4];
		    TT = TP * TS;
		    TY = TP * TX;
		    TU = W[5];
		    cr[WS(rs, 3)] = FNMS(TU, TX, TT);
		    ci[WS(rs, 3)] = FMA(TU, TS, TY);
	       }
	       {
		    E T1s, T1v, T1t, T1w, T1r, T1u;
		    T1s = T1i + T1j;
		    T1v = T1o + T1n;
		    T1r = W[2];
		    T1t = T1r * T1s;
		    T1w = T1r * T1v;
		    T1u = W[3];
		    cr[WS(rs, 2)] = FNMS(T1u, T1v, T1t);
		    ci[WS(rs, 2)] = FMA(T1u, T1s, T1w);
	       }
	       {
		    E T10, T13, T11, T14, TZ, T12;
		    T10 = FMA(KP707106781, TR, TQ);
		    T13 = FNMS(KP707106781, TW, TV);
		    TZ = W[12];
		    T11 = TZ * T10;
		    T14 = TZ * T13;
		    T12 = W[13];
		    cr[WS(rs, 7)] = FNMS(T12, T13, T11);
		    ci[WS(rs, 7)] = FMA(T12, T10, T14);
	       }
	       {
		    E T1f, T15, T17, T18, T1g, T16;
		    T1f = T1b - T1e;
		    T16 = T7 - Te;
		    T15 = W[6];
		    T17 = T15 * T16;
		    T18 = W[7];
		    T1g = T18 * T16;
		    cr[WS(rs, 4)] = FNMS(T18, T1f, T17);
		    ci[WS(rs, 4)] = FMA(T15, T1f, T1g);
	       }
	       {
		    E T1k, T1p, T1l, T1q, T1h, T1m;
		    T1k = T1i - T1j;
		    T1p = T1n - T1o;
		    T1h = W[10];
		    T1l = T1h * T1k;
		    T1q = T1h * T1p;
		    T1m = W[11];
		    cr[WS(rs, 6)] = FNMS(T1m, T1p, T1l);
		    ci[WS(rs, 6)] = FMA(T1m, T1k, T1q);
	       }
	       {
		    E TH, TN, TJ, TL, TM, TO, Tf, Tx, Ty, TI, TG, TK, Tw;
		    TG = TE - TF;
		    TH = FNMS(KP707106781, TG, TD);
		    TN = FMA(KP707106781, TG, TD);
		    TK = FMA(KP707106781, Tv, Tk);
		    TJ = W[0];
		    TL = TJ * TK;
		    TM = W[1];
		    TO = TM * TK;
		    Tw = FNMS(KP707106781, Tv, Tk);
		    Tf = W[8];
		    Tx = Tf * Tw;
		    Ty = W[9];
		    TI = Ty * Tw;
		    cr[WS(rs, 5)] = FNMS(Ty, TH, Tx);
		    ci[WS(rs, 5)] = FMA(Tf, TH, TI);
		    cr[WS(rs, 1)] = FNMS(TM, TN, TL);
		    ci[WS(rs, 1)] = FMA(TJ, TN, TO);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     { TW_FULL, 1, 8 },
     { TW_NEXT, 1, 0 }
};

static const hc2hc_desc desc = { 8, "hb_8", twinstr, &GENUS, { 44, 14, 22, 0 } };

void X(codelet_hb_8) (planner *p) {
     X(khc2hc_register) (p, hb_8, &desc);
}
#else

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hb_8 -include rdft/scalar/hb.h */

/*
 * This function contains 66 FP additions, 32 FP multiplications,
 * (or, 52 additions, 18 multiplications, 14 fused multiply/add),
 * 30 stack variables, 1 constants, and 32 memory accesses
 */
#include "rdft/scalar/hb.h"

static void hb_8(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 14); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 14, MAKE_VOLATILE_STRIDE(16, rs)) {
	       E T7, T18, T1c, To, Ty, TM, TY, TC, Te, TZ, T10, Tv, Tz, TP, TS;
	       E TD;
	       {
		    E T3, TK, Tn, TL, T6, TW, Tk, TX;
		    {
			 E T1, T2, Tl, Tm;
			 T1 = cr[0];
			 T2 = ci[WS(rs, 3)];
			 T3 = T1 + T2;
			 TK = T1 - T2;
			 Tl = ci[WS(rs, 5)];
			 Tm = cr[WS(rs, 6)];
			 Tn = Tl - Tm;
			 TL = Tl + Tm;
		    }
		    {
			 E T4, T5, Ti, Tj;
			 T4 = cr[WS(rs, 2)];
			 T5 = ci[WS(rs, 1)];
			 T6 = T4 + T5;
			 TW = T4 - T5;
			 Ti = ci[WS(rs, 7)];
			 Tj = cr[WS(rs, 4)];
			 Tk = Ti - Tj;
			 TX = Ti + Tj;
		    }
		    T7 = T3 + T6;
		    T18 = TK + TL;
		    T1c = TX - TW;
		    To = Tk + Tn;
		    Ty = T3 - T6;
		    TM = TK - TL;
		    TY = TW + TX;
		    TC = Tk - Tn;
	       }
	       {
		    E Ta, TN, Tu, TR, Td, TQ, Tr, TO;
		    {
			 E T8, T9, Ts, Tt;
			 T8 = cr[WS(rs, 1)];
			 T9 = ci[WS(rs, 2)];
			 Ta = T8 + T9;
			 TN = T8 - T9;
			 Ts = ci[WS(rs, 4)];
			 Tt = cr[WS(rs, 7)];
			 Tu = Ts - Tt;
			 TR = Ts + Tt;
		    }
		    {
			 E Tb, Tc, Tp, Tq;
			 Tb = ci[0];
			 Tc = cr[WS(rs, 3)];
			 Td = Tb + Tc;
			 TQ = Tb - Tc;
			 Tp = ci[WS(rs, 6)];
			 Tq = cr[WS(rs, 5)];
			 Tr = Tp - Tq;
			 TO = Tp + Tq;
		    }
		    Te = Ta + Td;
		    TZ = TN + TO;
		    T10 = TQ + TR;
		    Tv = Tr + Tu;
		    Tz = Tu - Tr;
		    TP = TN - TO;
		    TS = TQ - TR;
		    TD = Ta - Td;
	       }
	       cr[0] = T7 + Te;
	       ci[0] = To + Tv;
	       {
		    E Tg, Tw, Tf, Th;
		    Tg = T7 - Te;
		    Tw = To - Tv;
		    Tf = W[6];
		    Th = W[7];
		    cr[WS(rs, 4)] = FNMS(Th, Tw, Tf * Tg);
		    ci[WS(rs, 4)] = FMA(Th, Tg, Tf * Tw);
	       }
	       {
		    E TG, TI, TF, TH;
		    TG = Ty + Tz;
		    TI = TD + TC;
		    TF = W[2];
		    TH = W[3];
		    cr[WS(rs, 2)] = FNMS(TH, TI, TF * TG);
		    ci[WS(rs, 2)] = FMA(TF, TI, TH * TG);
	       }
	       {
		    E TA, TE, Tx, TB;
		    TA = Ty - Tz;
		    TE = TC - TD;
		    Tx = W[10];
		    TB = W[11];
		    cr[WS(rs, 6)] = FNMS(TB, TE, Tx * TA);
		    ci[WS(rs, 6)] = FMA(Tx, TE, TB * TA);
	       }
	       {
		    E T1a, T1g, T1e, T1i, T19, T1d;
		    T19 = KP707106781 * (TZ + T10);
		    T1a = T18 - T19;
		    T1g = T18 + T19;
		    T1d = KP707106781 * (TP - TS);
		    T1e = T1c + T1d;
		    T1i = T1c - T1d;
		    {
			 E T17, T1b, T1f, T1h;
			 T17 = W[4];
			 T1b = W[5];
			 cr[WS(rs, 3)] = FNMS(T1b, T1e, T17 * T1a);
			 ci[WS(rs, 3)] = FMA(T17, T1e, T1b * T1a);
			 T1f = W[12];
			 T1h = W[13];
			 cr[WS(rs, 7)] = FNMS(T1h, T1i, T1f * T1g);
			 ci[WS(rs, 7)] = FMA(T1f, T1i, T1h * T1g);
		    }
	       }
	       {
		    E TU, T14, T12, T16, TT, T11;
		    TT = KP707106781 * (TP + TS);
		    TU = TM - TT;
		    T14 = TM + TT;
		    T11 = KP707106781 * (TZ - T10);
		    T12 = TY - T11;
		    T16 = TY + T11;
		    {
			 E TJ, TV, T13, T15;
			 TJ = W[8];
			 TV = W[9];
			 cr[WS(rs, 5)] = FNMS(TV, T12, TJ * TU);
			 ci[WS(rs, 5)] = FMA(TV, TU, TJ * T12);
			 T13 = W[0];
			 T15 = W[1];
			 cr[WS(rs, 1)] = FNMS(T15, T16, T13 * T14);
			 ci[WS(rs, 1)] = FMA(T15, T14, T13 * T16);
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     { TW_FULL, 1, 8 },
     { TW_NEXT, 1, 0 }
};

static const hc2hc_desc desc = { 8, "hb_8", twinstr, &GENUS, { 52, 18, 14, 0 } };

void X(codelet_hb_8) (planner *p) {
     X(khc2hc_register) (p, hb_8, &desc);
}
#endif
