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
/* Generated on Thu Dec 10 07:06:49 EST 2020 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2cdft.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 20 -dif -name hc2cbdft_20 -include rdft/scalar/hc2cb.h */

/*
 * This function contains 286 FP additions, 148 FP multiplications,
 * (or, 176 additions, 38 multiplications, 110 fused multiply/add),
 * 104 stack variables, 4 constants, and 80 memory accesses
 */
#include "rdft/scalar/hc2cb.h"

static void hc2cbdft_20(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP618033988, +0.618033988749894848204586834365638117720309180);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 38); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 38, MAKE_VOLATILE_STRIDE(80, rs)) {
	       E T27, T2o, T3T, T41, T2p, T40, T1N, T2Q, T1w, T2L, T4n, T59, T4A, T5e, T24;
	       E T2m, T2h, T2Z, T3P, T4J, T3W, T3Y, T7, TC, T2c, T2d, T3y, T3F, T3G, T3H;
	       E T46, T4d, T4e, T4f, T4r, T4u, T4v, T4w, T1E, T1H, T1I, T1J, TJ, T16, T17;
	       E T18;
	       {
		    E T3, T1A, TI, T25, T6, TF, T1D, T26, Te, T47, T4k, TO, T1e, T3z, T3M;
		    E T1S, Tt, T4a, T4h, TZ, T1p, T3C, T3J, T1Z, TA, T4b, T4i, T14, T1u, T3D;
		    E T3K, T22, Tl, T48, T4l, TT, T1j, T3A, T3N, T1V;
		    {
			 E T1, T2, TG, TH;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 9)];
			 T3 = T1 + T2;
			 T1A = T1 - T2;
			 TG = Ip[0];
			 TH = Im[WS(rs, 9)];
			 TI = TG + TH;
			 T25 = TG - TH;
		    }
		    {
			 E T4, T5, T1B, T1C;
			 T4 = Rp[WS(rs, 5)];
			 T5 = Rm[WS(rs, 4)];
			 T6 = T4 + T5;
			 TF = T4 - T5;
			 T1B = Ip[WS(rs, 5)];
			 T1C = Im[WS(rs, 4)];
			 T1D = T1B + T1C;
			 T26 = T1B - T1C;
		    }
		    {
			 E Ta, T1a, TN, T1Q, Td, TK, T1d, T1R;
			 {
			      E T8, T9, TL, TM;
			      T8 = Rp[WS(rs, 4)];
			      T9 = Rm[WS(rs, 5)];
			      Ta = T8 + T9;
			      T1a = T8 - T9;
			      TL = Ip[WS(rs, 4)];
			      TM = Im[WS(rs, 5)];
			      TN = TL + TM;
			      T1Q = TL - TM;
			 }
			 {
			      E Tb, Tc, T1b, T1c;
			      Tb = Rp[WS(rs, 9)];
			      Tc = Rm[0];
			      Td = Tb + Tc;
			      TK = Tb - Tc;
			      T1b = Ip[WS(rs, 9)];
			      T1c = Im[0];
			      T1d = T1b + T1c;
			      T1R = T1b - T1c;
			 }
			 Te = Ta + Td;
			 T47 = TN - TK;
			 T4k = T1a + T1d;
			 TO = TK + TN;
			 T1e = T1a - T1d;
			 T3z = Ta - Td;
			 T3M = T1Q - T1R;
			 T1S = T1Q + T1R;
		    }
		    {
			 E Tp, T1l, TY, T1X, Ts, TV, T1o, T1Y;
			 {
			      E Tn, To, TW, TX;
			      Tn = Rp[WS(rs, 8)];
			      To = Rm[WS(rs, 1)];
			      Tp = Tn + To;
			      T1l = Tn - To;
			      TW = Ip[WS(rs, 8)];
			      TX = Im[WS(rs, 1)];
			      TY = TW + TX;
			      T1X = TW - TX;
			 }
			 {
			      E Tq, Tr, T1m, T1n;
			      Tq = Rm[WS(rs, 6)];
			      Tr = Rp[WS(rs, 3)];
			      Ts = Tq + Tr;
			      TV = Tq - Tr;
			      T1m = Im[WS(rs, 6)];
			      T1n = Ip[WS(rs, 3)];
			      T1o = T1m + T1n;
			      T1Y = T1n - T1m;
			 }
			 Tt = Tp + Ts;
			 T4a = TY - TV;
			 T4h = T1l - T1o;
			 TZ = TV + TY;
			 T1p = T1l + T1o;
			 T3C = Tp - Ts;
			 T3J = T1X - T1Y;
			 T1Z = T1X + T1Y;
		    }
		    {
			 E Tw, T1q, T13, T20, Tz, T10, T1t, T21;
			 {
			      E Tu, Tv, T11, T12;
			      Tu = Rm[WS(rs, 7)];
			      Tv = Rp[WS(rs, 2)];
			      Tw = Tu + Tv;
			      T1q = Tu - Tv;
			      T11 = Im[WS(rs, 7)];
			      T12 = Ip[WS(rs, 2)];
			      T13 = T11 + T12;
			      T20 = T12 - T11;
			 }
			 {
			      E Tx, Ty, T1r, T1s;
			      Tx = Rm[WS(rs, 2)];
			      Ty = Rp[WS(rs, 7)];
			      Tz = Tx + Ty;
			      T10 = Tx - Ty;
			      T1r = Im[WS(rs, 2)];
			      T1s = Ip[WS(rs, 7)];
			      T1t = T1r + T1s;
			      T21 = T1s - T1r;
			 }
			 TA = Tw + Tz;
			 T4b = T10 + T13;
			 T4i = T1q - T1t;
			 T14 = T10 - T13;
			 T1u = T1q + T1t;
			 T3D = Tw - Tz;
			 T3K = T20 - T21;
			 T22 = T20 + T21;
		    }
		    {
			 E Th, T1f, TS, T1T, Tk, TP, T1i, T1U;
			 {
			      E Tf, Tg, TQ, TR;
			      Tf = Rm[WS(rs, 3)];
			      Tg = Rp[WS(rs, 6)];
			      Th = Tf + Tg;
			      T1f = Tf - Tg;
			      TQ = Im[WS(rs, 3)];
			      TR = Ip[WS(rs, 6)];
			      TS = TQ + TR;
			      T1T = TR - TQ;
			 }
			 {
			      E Ti, Tj, T1g, T1h;
			      Ti = Rp[WS(rs, 1)];
			      Tj = Rm[WS(rs, 8)];
			      Tk = Ti + Tj;
			      TP = Ti - Tj;
			      T1g = Ip[WS(rs, 1)];
			      T1h = Im[WS(rs, 8)];
			      T1i = T1g + T1h;
			      T1U = T1g - T1h;
			 }
			 Tl = Th + Tk;
			 T48 = TP + TS;
			 T4l = T1f + T1i;
			 TT = TP - TS;
			 T1j = T1f - T1i;
			 T3A = Th - Tk;
			 T3N = T1T - T1U;
			 T1V = T1T + T1U;
		    }
		    T27 = T25 + T26;
		    T2o = Tt - TA;
		    T3T = T25 - T26;
		    T41 = T3z - T3A;
		    T2p = Te - Tl;
		    {
			 E T1L, T1M, T1k, T1v;
			 T40 = T3C - T3D;
			 T1L = TO - TT;
			 T1M = TZ - T14;
			 T1N = FMA(KP618033988, T1M, T1L);
			 T2Q = FNMS(KP618033988, T1L, T1M);
			 T1k = T1e - T1j;
			 T1v = T1p - T1u;
			 T1w = FMA(KP618033988, T1v, T1k);
			 T2L = FNMS(KP618033988, T1k, T1v);
			 {
			      E T4j, T4m, T4y, T4z;
			      T4j = T4h - T4i;
			      T4m = T4k - T4l;
			      T4n = FNMS(KP618033988, T4m, T4j);
			      T59 = FMA(KP618033988, T4j, T4m);
			      T4y = T4a + T4b;
			      T4z = T47 + T48;
			      T4A = FNMS(KP618033988, T4z, T4y);
			      T5e = FMA(KP618033988, T4y, T4z);
			 }
		    }
		    {
			 E T3L, T3O, T4s, T4t;
			 {
			      E T1W, T23, T2f, T2g;
			      T1W = T1S + T1V;
			      T23 = T1Z + T22;
			      T24 = T1W + T23;
			      T2m = T1W - T23;
			      T2f = T1Z - T22;
			      T2g = T1S - T1V;
			      T2h = FNMS(KP618033988, T2g, T2f);
			      T2Z = FMA(KP618033988, T2f, T2g);
			 }
			 T3L = T3J - T3K;
			 T3O = T3M - T3N;
			 T3P = FNMS(KP618033988, T3O, T3L);
			 T4J = FMA(KP618033988, T3L, T3O);
			 {
			      E T3U, T3V, Tm, TB;
			      T3U = T3M + T3N;
			      T3V = T3J + T3K;
			      T3W = T3U + T3V;
			      T3Y = T3U - T3V;
			      T7 = T3 + T6;
			      Tm = Te + Tl;
			      TB = Tt + TA;
			      TC = Tm + TB;
			      T2c = FNMS(KP250000000, TC, T7);
			      T2d = Tm - TB;
			 }
			 {
			      E T3B, T3E, T49, T4c;
			      T3y = T3 - T6;
			      T3B = T3z + T3A;
			      T3E = T3C + T3D;
			      T3F = T3B + T3E;
			      T3G = FNMS(KP250000000, T3F, T3y);
			      T3H = T3B - T3E;
			      T46 = TI - TF;
			      T49 = T47 - T48;
			      T4c = T4a - T4b;
			      T4d = T49 + T4c;
			      T4e = FNMS(KP250000000, T4d, T46);
			      T4f = T49 - T4c;
			 }
			 T4r = T1A + T1D;
			 T4s = T4k + T4l;
			 T4t = T4h + T4i;
			 T4u = T4s + T4t;
			 T4v = FNMS(KP250000000, T4u, T4r);
			 T4w = T4s - T4t;
			 {
			      E T1F, T1G, TU, T15;
			      T1E = T1A - T1D;
			      T1F = T1e + T1j;
			      T1G = T1p + T1u;
			      T1H = T1F + T1G;
			      T1I = FNMS(KP250000000, T1H, T1E);
			      T1J = T1F - T1G;
			      TJ = TF + TI;
			      TU = TO + TT;
			      T15 = TZ + T14;
			      T16 = TU + T15;
			      T17 = FNMS(KP250000000, T16, TJ);
			      T18 = TU - T15;
			 }
		    }
	       }
	       {
		    E TD, T28, T3o, T3r, T3p, T3v, T2r, T3l, T2H, T35, T2b, T2j, T2k, T2z, T2D;
		    E T2F, T2G, T2T, T2X, T31, T32, T3d, T3h, T3j, T3k, T3t, T1x, T2u, T1O, T2x;
		    E T1y, T29, T2v, T2B, T2M, T38, T2R, T3b, T2N, T2V, T39, T3f, T3n, T1P, T2a;
		    E T1z;
		    TD = T7 + TC;
		    T28 = T24 + T27;
		    T3o = TJ + T16;
		    T3r = T1H + T1E;
		    T3n = W[8];
		    T3p = T3n * T3o;
		    T3v = T3n * T3r;
		    {
			 E T2q, T34, T2n, T33, T2l;
			 T2q = FNMS(KP618033988, T2p, T2o);
			 T34 = FMA(KP618033988, T2o, T2p);
			 T2l = FNMS(KP250000000, T24, T27);
			 T2n = FNMS(KP559016994, T2m, T2l);
			 T33 = FMA(KP559016994, T2m, T2l);
			 T2r = FMA(KP951056516, T2q, T2n);
			 T3l = FNMS(KP951056516, T34, T33);
			 T2H = FNMS(KP951056516, T2q, T2n);
			 T35 = FMA(KP951056516, T34, T33);
		    }
		    {
			 E T2i, T2E, T2e, T30, T3i, T2Y;
			 T2e = FNMS(KP559016994, T2d, T2c);
			 T2i = FNMS(KP951056516, T2h, T2e);
			 T2E = FMA(KP951056516, T2h, T2e);
			 T2b = W[14];
			 T2j = T2b * T2i;
			 T2k = W[15];
			 T2z = T2k * T2i;
			 T2D = W[22];
			 T2F = T2D * T2E;
			 T2G = W[23];
			 T2T = T2G * T2E;
			 T2Y = FMA(KP559016994, T2d, T2c);
			 T30 = FNMS(KP951056516, T2Z, T2Y);
			 T3i = FMA(KP951056516, T2Z, T2Y);
			 T2X = W[30];
			 T31 = T2X * T30;
			 T32 = W[31];
			 T3d = T32 * T30;
			 T3h = W[6];
			 T3j = T3h * T3i;
			 T3k = W[7];
			 T3t = T3k * T3i;
		    }
		    {
			 E T19, T1K, TE, T2t;
			 T19 = FMA(KP559016994, T18, T17);
			 T1x = FMA(KP951056516, T1w, T19);
			 T2u = FNMS(KP951056516, T1w, T19);
			 T1K = FMA(KP559016994, T1J, T1I);
			 T1O = FNMS(KP951056516, T1N, T1K);
			 T2x = FMA(KP951056516, T1N, T1K);
			 TE = W[0];
			 T1y = TE * T1x;
			 T29 = TE * T1O;
			 T2t = W[16];
			 T2v = T2t * T2u;
			 T2B = T2t * T2x;
		    }
		    {
			 E T2K, T2P, T2J, T37;
			 T2K = FNMS(KP559016994, T18, T17);
			 T2M = FMA(KP951056516, T2L, T2K);
			 T38 = FNMS(KP951056516, T2L, T2K);
			 T2P = FNMS(KP559016994, T1J, T1I);
			 T2R = FNMS(KP951056516, T2Q, T2P);
			 T3b = FMA(KP951056516, T2Q, T2P);
			 T2J = W[24];
			 T2N = T2J * T2M;
			 T2V = T2J * T2R;
			 T37 = W[32];
			 T39 = T37 * T38;
			 T3f = T37 * T3b;
		    }
		    T1z = W[1];
		    T1P = FMA(T1z, T1O, T1y);
		    T2a = FNMS(T1z, T1x, T29);
		    Rp[0] = TD - T1P;
		    Ip[0] = T28 + T2a;
		    Rm[0] = TD + T1P;
		    Im[0] = T2a - T28;
		    {
			 E T3m, T3u, T3s, T3w, T3q;
			 T3m = FNMS(T3k, T3l, T3j);
			 T3u = FMA(T3h, T3l, T3t);
			 T3q = W[9];
			 T3s = FMA(T3q, T3r, T3p);
			 T3w = FNMS(T3q, T3o, T3v);
			 Rp[WS(rs, 2)] = T3m - T3s;
			 Ip[WS(rs, 2)] = T3u + T3w;
			 Rm[WS(rs, 2)] = T3m + T3s;
			 Im[WS(rs, 2)] = T3w - T3u;
		    }
		    {
			 E T2s, T2A, T2y, T2C, T2w;
			 T2s = FNMS(T2k, T2r, T2j);
			 T2A = FMA(T2b, T2r, T2z);
			 T2w = W[17];
			 T2y = FMA(T2w, T2x, T2v);
			 T2C = FNMS(T2w, T2u, T2B);
			 Rp[WS(rs, 4)] = T2s - T2y;
			 Ip[WS(rs, 4)] = T2A + T2C;
			 Rm[WS(rs, 4)] = T2s + T2y;
			 Im[WS(rs, 4)] = T2C - T2A;
		    }
		    {
			 E T2I, T2U, T2S, T2W, T2O;
			 T2I = FNMS(T2G, T2H, T2F);
			 T2U = FMA(T2D, T2H, T2T);
			 T2O = W[25];
			 T2S = FMA(T2O, T2R, T2N);
			 T2W = FNMS(T2O, T2M, T2V);
			 Rp[WS(rs, 6)] = T2I - T2S;
			 Ip[WS(rs, 6)] = T2U + T2W;
			 Rm[WS(rs, 6)] = T2I + T2S;
			 Im[WS(rs, 6)] = T2W - T2U;
		    }
		    {
			 E T36, T3e, T3c, T3g, T3a;
			 T36 = FNMS(T32, T35, T31);
			 T3e = FMA(T2X, T35, T3d);
			 T3a = W[33];
			 T3c = FMA(T3a, T3b, T39);
			 T3g = FNMS(T3a, T38, T3f);
			 Rp[WS(rs, 8)] = T36 - T3c;
			 Ip[WS(rs, 8)] = T3e + T3g;
			 Rm[WS(rs, 8)] = T36 + T3c;
			 Im[WS(rs, 8)] = T3g - T3e;
		    }
	       }
	       {
		    E T55, T51, T53, T54, T5h, T5I, T5L, T5J, T5P, T43, T5F, T4P, T5p, T3x, T3R;
		    E T3S, T4D, T5l, T5n, T5o, T5x, T4H, T4L, T4M, T4X, T5B, T5D, T5E, T5N, T4o;
		    E T4S, T4B, T4V, T4p, T4F, T4T, T4Z, T5a, T5s, T5f, T5v, T5b, T5j, T5t, T5z;
		    E T52, T5H;
		    T55 = T3W + T3T;
		    T52 = T3y + T3F;
		    T51 = W[18];
		    T53 = T51 * T52;
		    T54 = W[19];
		    T5h = T54 * T52;
		    T5I = T46 + T4d;
		    T5L = T4u + T4r;
		    T5H = W[28];
		    T5J = T5H * T5I;
		    T5P = T5H * T5L;
		    {
			 E T42, T4O, T3Z, T4N, T3X;
			 T42 = FNMS(KP618033988, T41, T40);
			 T4O = FMA(KP618033988, T40, T41);
			 T3X = FNMS(KP250000000, T3W, T3T);
			 T3Z = FNMS(KP559016994, T3Y, T3X);
			 T4N = FMA(KP559016994, T3Y, T3X);
			 T43 = FNMS(KP951056516, T42, T3Z);
			 T5F = FNMS(KP951056516, T4O, T4N);
			 T4P = FMA(KP951056516, T4O, T4N);
			 T5p = FMA(KP951056516, T42, T3Z);
		    }
		    {
			 E T3Q, T5m, T3I, T4K, T5C, T4I;
			 T3I = FNMS(KP559016994, T3H, T3G);
			 T3Q = FMA(KP951056516, T3P, T3I);
			 T5m = FNMS(KP951056516, T3P, T3I);
			 T3x = W[2];
			 T3R = T3x * T3Q;
			 T3S = W[3];
			 T4D = T3S * T3Q;
			 T5l = W[34];
			 T5n = T5l * T5m;
			 T5o = W[35];
			 T5x = T5o * T5m;
			 T4I = FMA(KP559016994, T3H, T3G);
			 T4K = FNMS(KP951056516, T4J, T4I);
			 T5C = FMA(KP951056516, T4J, T4I);
			 T4H = W[10];
			 T4L = T4H * T4K;
			 T4M = W[11];
			 T4X = T4M * T4K;
			 T5B = W[26];
			 T5D = T5B * T5C;
			 T5E = W[27];
			 T5N = T5E * T5C;
		    }
		    {
			 E T4g, T4x, T45, T4R;
			 T4g = FNMS(KP559016994, T4f, T4e);
			 T4o = FMA(KP951056516, T4n, T4g);
			 T4S = FNMS(KP951056516, T4n, T4g);
			 T4x = FNMS(KP559016994, T4w, T4v);
			 T4B = FNMS(KP951056516, T4A, T4x);
			 T4V = FMA(KP951056516, T4A, T4x);
			 T45 = W[4];
			 T4p = T45 * T4o;
			 T4F = T45 * T4B;
			 T4R = W[12];
			 T4T = T4R * T4S;
			 T4Z = T4R * T4V;
		    }
		    {
			 E T58, T5d, T57, T5r;
			 T58 = FMA(KP559016994, T4f, T4e);
			 T5a = FMA(KP951056516, T59, T58);
			 T5s = FNMS(KP951056516, T59, T58);
			 T5d = FMA(KP559016994, T4w, T4v);
			 T5f = FNMS(KP951056516, T5e, T5d);
			 T5v = FMA(KP951056516, T5e, T5d);
			 T57 = W[20];
			 T5b = T57 * T5a;
			 T5j = T57 * T5f;
			 T5r = W[36];
			 T5t = T5r * T5s;
			 T5z = T5r * T5v;
		    }
		    {
			 E T44, T4E, T4C, T4G, T4q;
			 T44 = FNMS(T3S, T43, T3R);
			 T4E = FMA(T3x, T43, T4D);
			 T4q = W[5];
			 T4C = FMA(T4q, T4B, T4p);
			 T4G = FNMS(T4q, T4o, T4F);
			 Rp[WS(rs, 1)] = T44 - T4C;
			 Ip[WS(rs, 1)] = T4E + T4G;
			 Rm[WS(rs, 1)] = T44 + T4C;
			 Im[WS(rs, 1)] = T4G - T4E;
		    }
		    {
			 E T5G, T5O, T5M, T5Q, T5K;
			 T5G = FNMS(T5E, T5F, T5D);
			 T5O = FMA(T5B, T5F, T5N);
			 T5K = W[29];
			 T5M = FMA(T5K, T5L, T5J);
			 T5Q = FNMS(T5K, T5I, T5P);
			 Rp[WS(rs, 7)] = T5G - T5M;
			 Ip[WS(rs, 7)] = T5O + T5Q;
			 Rm[WS(rs, 7)] = T5G + T5M;
			 Im[WS(rs, 7)] = T5Q - T5O;
		    }
		    {
			 E T4Q, T4Y, T4W, T50, T4U;
			 T4Q = FNMS(T4M, T4P, T4L);
			 T4Y = FMA(T4H, T4P, T4X);
			 T4U = W[13];
			 T4W = FMA(T4U, T4V, T4T);
			 T50 = FNMS(T4U, T4S, T4Z);
			 Rp[WS(rs, 3)] = T4Q - T4W;
			 Ip[WS(rs, 3)] = T4Y + T50;
			 Rm[WS(rs, 3)] = T4Q + T4W;
			 Im[WS(rs, 3)] = T50 - T4Y;
		    }
		    {
			 E T56, T5i, T5g, T5k, T5c;
			 T56 = FNMS(T54, T55, T53);
			 T5i = FMA(T51, T55, T5h);
			 T5c = W[21];
			 T5g = FMA(T5c, T5f, T5b);
			 T5k = FNMS(T5c, T5a, T5j);
			 Rp[WS(rs, 5)] = T56 - T5g;
			 Ip[WS(rs, 5)] = T5i + T5k;
			 Rm[WS(rs, 5)] = T56 + T5g;
			 Im[WS(rs, 5)] = T5k - T5i;
		    }
		    {
			 E T5q, T5y, T5w, T5A, T5u;
			 T5q = FNMS(T5o, T5p, T5n);
			 T5y = FMA(T5l, T5p, T5x);
			 T5u = W[37];
			 T5w = FMA(T5u, T5v, T5t);
			 T5A = FNMS(T5u, T5s, T5z);
			 Rp[WS(rs, 9)] = T5q - T5w;
			 Ip[WS(rs, 9)] = T5y + T5A;
			 Rm[WS(rs, 9)] = T5q + T5w;
			 Im[WS(rs, 9)] = T5A - T5y;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     { TW_FULL, 1, 20 },
     { TW_NEXT, 1, 0 }
};

static const hc2c_desc desc = { 20, "hc2cbdft_20", twinstr, &GENUS, { 176, 38, 110, 0 } };

void X(codelet_hc2cbdft_20) (planner *p) {
     X(khc2c_register) (p, hc2cbdft_20, &desc, HC2C_VIA_DFT);
}
#else

/* Generated by: ../../../genfft/gen_hc2cdft.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 20 -dif -name hc2cbdft_20 -include rdft/scalar/hc2cb.h */

/*
 * This function contains 286 FP additions, 124 FP multiplications,
 * (or, 224 additions, 62 multiplications, 62 fused multiply/add),
 * 89 stack variables, 4 constants, and 80 memory accesses
 */
#include "rdft/scalar/hc2cb.h"

static void hc2cbdft_20(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 38); m < me; m = m + 1, Rp = Rp + ms, Ip = Ip + ms, Rm = Rm - ms, Im = Im - ms, W = W + 38, MAKE_VOLATILE_STRIDE(80, rs)) {
	       E T7, T3N, T4a, T16, T1G, T3g, T3D, T26, T1k, T3A, T3B, T1v, T2e, T48, T47;
	       E T2d, T1L, T43, T40, T1K, T2l, T3t, T2m, T3w, T3n, T3p, TC, T2b, T4d, T4f;
	       E T23, T2j, T1B, T1H, T3U, T3W, T3G, T3I, T11, T17;
	       {
		    E T3, T1C, T15, T24, T6, T12, T1F, T25;
		    {
			 E T1, T2, T13, T14;
			 T1 = Rp[0];
			 T2 = Rm[WS(rs, 9)];
			 T3 = T1 + T2;
			 T1C = T1 - T2;
			 T13 = Ip[0];
			 T14 = Im[WS(rs, 9)];
			 T15 = T13 + T14;
			 T24 = T13 - T14;
		    }
		    {
			 E T4, T5, T1D, T1E;
			 T4 = Rp[WS(rs, 5)];
			 T5 = Rm[WS(rs, 4)];
			 T6 = T4 + T5;
			 T12 = T4 - T5;
			 T1D = Ip[WS(rs, 5)];
			 T1E = Im[WS(rs, 4)];
			 T1F = T1D + T1E;
			 T25 = T1D - T1E;
		    }
		    T7 = T3 + T6;
		    T3N = T15 - T12;
		    T4a = T1C + T1F;
		    T16 = T12 + T15;
		    T1G = T1C - T1F;
		    T3g = T3 - T6;
		    T3D = T24 - T25;
		    T26 = T24 + T25;
	       }
	       {
		    E Te, T3O, T3Y, TJ, T1e, T3h, T3r, T1R, TA, T3S, T42, TZ, T1u, T3l, T3v;
		    E T21, Tl, T3P, T3Z, TO, T1j, T3i, T3s, T1U, Tt, T3R, T41, TU, T1p, T3k;
		    E T3u, T1Y;
		    {
			 E Ta, T1a, TI, T1P, Td, TF, T1d, T1Q;
			 {
			      E T8, T9, TG, TH;
			      T8 = Rp[WS(rs, 4)];
			      T9 = Rm[WS(rs, 5)];
			      Ta = T8 + T9;
			      T1a = T8 - T9;
			      TG = Ip[WS(rs, 4)];
			      TH = Im[WS(rs, 5)];
			      TI = TG + TH;
			      T1P = TG - TH;
			 }
			 {
			      E Tb, Tc, T1b, T1c;
			      Tb = Rp[WS(rs, 9)];
			      Tc = Rm[0];
			      Td = Tb + Tc;
			      TF = Tb - Tc;
			      T1b = Ip[WS(rs, 9)];
			      T1c = Im[0];
			      T1d = T1b + T1c;
			      T1Q = T1b - T1c;
			 }
			 Te = Ta + Td;
			 T3O = TI - TF;
			 T3Y = T1a + T1d;
			 TJ = TF + TI;
			 T1e = T1a - T1d;
			 T3h = Ta - Td;
			 T3r = T1P - T1Q;
			 T1R = T1P + T1Q;
		    }
		    {
			 E Tw, T1q, TY, T1Z, Tz, TV, T1t, T20;
			 {
			      E Tu, Tv, TW, TX;
			      Tu = Rm[WS(rs, 7)];
			      Tv = Rp[WS(rs, 2)];
			      Tw = Tu + Tv;
			      T1q = Tu - Tv;
			      TW = Im[WS(rs, 7)];
			      TX = Ip[WS(rs, 2)];
			      TY = TW + TX;
			      T1Z = TX - TW;
			 }
			 {
			      E Tx, Ty, T1r, T1s;
			      Tx = Rm[WS(rs, 2)];
			      Ty = Rp[WS(rs, 7)];
			      Tz = Tx + Ty;
			      TV = Tx - Ty;
			      T1r = Im[WS(rs, 2)];
			      T1s = Ip[WS(rs, 7)];
			      T1t = T1r + T1s;
			      T20 = T1s - T1r;
			 }
			 TA = Tw + Tz;
			 T3S = TV + TY;
			 T42 = T1q - T1t;
			 TZ = TV - TY;
			 T1u = T1q + T1t;
			 T3l = Tw - Tz;
			 T3v = T1Z - T20;
			 T21 = T1Z + T20;
		    }
		    {
			 E Th, T1f, TN, T1S, Tk, TK, T1i, T1T;
			 {
			      E Tf, Tg, TL, TM;
			      Tf = Rm[WS(rs, 3)];
			      Tg = Rp[WS(rs, 6)];
			      Th = Tf + Tg;
			      T1f = Tf - Tg;
			      TL = Im[WS(rs, 3)];
			      TM = Ip[WS(rs, 6)];
			      TN = TL + TM;
			      T1S = TM - TL;
			 }
			 {
			      E Ti, Tj, T1g, T1h;
			      Ti = Rp[WS(rs, 1)];
			      Tj = Rm[WS(rs, 8)];
			      Tk = Ti + Tj;
			      TK = Ti - Tj;
			      T1g = Ip[WS(rs, 1)];
			      T1h = Im[WS(rs, 8)];
			      T1i = T1g + T1h;
			      T1T = T1g - T1h;
			 }
			 Tl = Th + Tk;
			 T3P = TK + TN;
			 T3Z = T1f + T1i;
			 TO = TK - TN;
			 T1j = T1f - T1i;
			 T3i = Th - Tk;
			 T3s = T1S - T1T;
			 T1U = T1S + T1T;
		    }
		    {
			 E Tp, T1l, TT, T1W, Ts, TQ, T1o, T1X;
			 {
			      E Tn, To, TR, TS;
			      Tn = Rp[WS(rs, 8)];
			      To = Rm[WS(rs, 1)];
			      Tp = Tn + To;
			      T1l = Tn - To;
			      TR = Ip[WS(rs, 8)];
			      TS = Im[WS(rs, 1)];
			      TT = TR + TS;
			      T1W = TR - TS;
			 }
			 {
			      E Tq, Tr, T1m, T1n;
			      Tq = Rm[WS(rs, 6)];
			      Tr = Rp[WS(rs, 3)];
			      Ts = Tq + Tr;
			      TQ = Tq - Tr;
			      T1m = Im[WS(rs, 6)];
			      T1n = Ip[WS(rs, 3)];
			      T1o = T1m + T1n;
			      T1X = T1n - T1m;
			 }
			 Tt = Tp + Ts;
			 T3R = TT - TQ;
			 T41 = T1l - T1o;
			 TU = TQ + TT;
			 T1p = T1l + T1o;
			 T3k = Tp - Ts;
			 T3u = T1W - T1X;
			 T1Y = T1W + T1X;
		    }
		    T1k = T1e - T1j;
		    T3A = T3h - T3i;
		    T3B = T3k - T3l;
		    T1v = T1p - T1u;
		    T2e = T1Y - T21;
		    T48 = T3R + T3S;
		    T47 = T3O + T3P;
		    T2d = T1R - T1U;
		    T1L = TU - TZ;
		    T43 = T41 - T42;
		    T40 = T3Y - T3Z;
		    T1K = TJ - TO;
		    T2l = Te - Tl;
		    T3t = T3r - T3s;
		    T2m = Tt - TA;
		    T3w = T3u - T3v;
		    {
			 E T3j, T3m, Tm, TB;
			 T3j = T3h + T3i;
			 T3m = T3k + T3l;
			 T3n = T3j + T3m;
			 T3p = KP559016994 * (T3j - T3m);
			 Tm = Te + Tl;
			 TB = Tt + TA;
			 TC = Tm + TB;
			 T2b = KP559016994 * (Tm - TB);
		    }
		    {
			 E T4b, T4c, T3Q, T3T;
			 T4b = T3Y + T3Z;
			 T4c = T41 + T42;
			 T4d = T4b + T4c;
			 T4f = KP559016994 * (T4b - T4c);
			 {
			      E T1V, T22, T1z, T1A;
			      T1V = T1R + T1U;
			      T22 = T1Y + T21;
			      T23 = T1V + T22;
			      T2j = KP559016994 * (T1V - T22);
			      T1z = T1e + T1j;
			      T1A = T1p + T1u;
			      T1B = KP559016994 * (T1z - T1A);
			      T1H = T1z + T1A;
			 }
			 T3Q = T3O - T3P;
			 T3T = T3R - T3S;
			 T3U = T3Q + T3T;
			 T3W = KP559016994 * (T3Q - T3T);
			 {
			      E T3E, T3F, TP, T10;
			      T3E = T3r + T3s;
			      T3F = T3u + T3v;
			      T3G = T3E + T3F;
			      T3I = KP559016994 * (T3E - T3F);
			      TP = TJ + TO;
			      T10 = TU + TZ;
			      T11 = KP559016994 * (TP - T10);
			      T17 = TP + T10;
			 }
		    }
	       }
	       {
		    E TD, T27, T3c, T3e, T2o, T36, T2A, T2U, T1N, T2Z, T2t, T2J, T1x, T2X, T2r;
		    E T2F, T2g, T34, T2y, T2Q;
		    TD = T7 + TC;
		    T27 = T23 + T26;
		    {
			 E T39, T3b, T38, T3a;
			 T39 = T16 + T17;
			 T3b = T1H + T1G;
			 T38 = W[8];
			 T3a = W[9];
			 T3c = FMA(T38, T39, T3a * T3b);
			 T3e = FNMS(T3a, T39, T38 * T3b);
		    }
		    {
			 E T2n, T2S, T2k, T2T, T2i;
			 T2n = FNMS(KP951056516, T2m, KP587785252 * T2l);
			 T2S = FMA(KP951056516, T2l, KP587785252 * T2m);
			 T2i = FNMS(KP250000000, T23, T26);
			 T2k = T2i - T2j;
			 T2T = T2j + T2i;
			 T2o = T2k - T2n;
			 T36 = T2T - T2S;
			 T2A = T2n + T2k;
			 T2U = T2S + T2T;
		    }
		    {
			 E T1M, T2H, T1J, T2I, T1I;
			 T1M = FMA(KP951056516, T1K, KP587785252 * T1L);
			 T2H = FNMS(KP951056516, T1L, KP587785252 * T1K);
			 T1I = FNMS(KP250000000, T1H, T1G);
			 T1J = T1B + T1I;
			 T2I = T1I - T1B;
			 T1N = T1J - T1M;
			 T2Z = T2I - T2H;
			 T2t = T1M + T1J;
			 T2J = T2H + T2I;
		    }
		    {
			 E T1w, T2E, T19, T2D, T18;
			 T1w = FMA(KP951056516, T1k, KP587785252 * T1v);
			 T2E = FNMS(KP951056516, T1v, KP587785252 * T1k);
			 T18 = FNMS(KP250000000, T17, T16);
			 T19 = T11 + T18;
			 T2D = T18 - T11;
			 T1x = T19 + T1w;
			 T2X = T2D + T2E;
			 T2r = T19 - T1w;
			 T2F = T2D - T2E;
		    }
		    {
			 E T2f, T2P, T2c, T2O, T2a;
			 T2f = FNMS(KP951056516, T2e, KP587785252 * T2d);
			 T2P = FMA(KP951056516, T2d, KP587785252 * T2e);
			 T2a = FNMS(KP250000000, TC, T7);
			 T2c = T2a - T2b;
			 T2O = T2b + T2a;
			 T2g = T2c + T2f;
			 T34 = T2O + T2P;
			 T2y = T2c - T2f;
			 T2Q = T2O - T2P;
		    }
		    {
			 E T1O, T28, TE, T1y;
			 TE = W[0];
			 T1y = W[1];
			 T1O = FMA(TE, T1x, T1y * T1N);
			 T28 = FNMS(T1y, T1x, TE * T1N);
			 Rp[0] = TD - T1O;
			 Ip[0] = T27 + T28;
			 Rm[0] = TD + T1O;
			 Im[0] = T28 - T27;
		    }
		    {
			 E T37, T3d, T33, T35;
			 T33 = W[6];
			 T35 = W[7];
			 T37 = FNMS(T35, T36, T33 * T34);
			 T3d = FMA(T35, T34, T33 * T36);
			 Rp[WS(rs, 2)] = T37 - T3c;
			 Ip[WS(rs, 2)] = T3d + T3e;
			 Rm[WS(rs, 2)] = T37 + T3c;
			 Im[WS(rs, 2)] = T3e - T3d;
		    }
		    {
			 E T2p, T2v, T2u, T2w;
			 {
			      E T29, T2h, T2q, T2s;
			      T29 = W[14];
			      T2h = W[15];
			      T2p = FNMS(T2h, T2o, T29 * T2g);
			      T2v = FMA(T2h, T2g, T29 * T2o);
			      T2q = W[16];
			      T2s = W[17];
			      T2u = FMA(T2q, T2r, T2s * T2t);
			      T2w = FNMS(T2s, T2r, T2q * T2t);
			 }
			 Rp[WS(rs, 4)] = T2p - T2u;
			 Ip[WS(rs, 4)] = T2v + T2w;
			 Rm[WS(rs, 4)] = T2p + T2u;
			 Im[WS(rs, 4)] = T2w - T2v;
		    }
		    {
			 E T2B, T2L, T2K, T2M;
			 {
			      E T2x, T2z, T2C, T2G;
			      T2x = W[22];
			      T2z = W[23];
			      T2B = FNMS(T2z, T2A, T2x * T2y);
			      T2L = FMA(T2z, T2y, T2x * T2A);
			      T2C = W[24];
			      T2G = W[25];
			      T2K = FMA(T2C, T2F, T2G * T2J);
			      T2M = FNMS(T2G, T2F, T2C * T2J);
			 }
			 Rp[WS(rs, 6)] = T2B - T2K;
			 Ip[WS(rs, 6)] = T2L + T2M;
			 Rm[WS(rs, 6)] = T2B + T2K;
			 Im[WS(rs, 6)] = T2M - T2L;
		    }
		    {
			 E T2V, T31, T30, T32;
			 {
			      E T2N, T2R, T2W, T2Y;
			      T2N = W[30];
			      T2R = W[31];
			      T2V = FNMS(T2R, T2U, T2N * T2Q);
			      T31 = FMA(T2R, T2Q, T2N * T2U);
			      T2W = W[32];
			      T2Y = W[33];
			      T30 = FMA(T2W, T2X, T2Y * T2Z);
			      T32 = FNMS(T2Y, T2X, T2W * T2Z);
			 }
			 Rp[WS(rs, 8)] = T2V - T30;
			 Ip[WS(rs, 8)] = T31 + T32;
			 Rm[WS(rs, 8)] = T2V + T30;
			 Im[WS(rs, 8)] = T32 - T31;
		    }
	       }
	       {
		    E T4F, T4P, T5c, T5e, T3y, T54, T4o, T4S, T4h, T4Z, T4x, T4N, T45, T4X, T4v;
		    E T4J, T3K, T56, T4s, T4U;
		    {
			 E T4C, T4E, T4B, T4D;
			 T4C = T3g + T3n;
			 T4E = T3G + T3D;
			 T4B = W[18];
			 T4D = W[19];
			 T4F = FNMS(T4D, T4E, T4B * T4C);
			 T4P = FMA(T4D, T4C, T4B * T4E);
		    }
		    {
			 E T59, T5b, T58, T5a;
			 T59 = T3N + T3U;
			 T5b = T4d + T4a;
			 T58 = W[28];
			 T5a = W[29];
			 T5c = FMA(T58, T59, T5a * T5b);
			 T5e = FNMS(T5a, T59, T58 * T5b);
		    }
		    {
			 E T3x, T4n, T3q, T4m, T3o;
			 T3x = FNMS(KP951056516, T3w, KP587785252 * T3t);
			 T4n = FMA(KP951056516, T3t, KP587785252 * T3w);
			 T3o = FNMS(KP250000000, T3n, T3g);
			 T3q = T3o - T3p;
			 T4m = T3p + T3o;
			 T3y = T3q - T3x;
			 T54 = T4m + T4n;
			 T4o = T4m - T4n;
			 T4S = T3q + T3x;
		    }
		    {
			 E T49, T4M, T4g, T4L, T4e;
			 T49 = FNMS(KP951056516, T48, KP587785252 * T47);
			 T4M = FMA(KP951056516, T47, KP587785252 * T48);
			 T4e = FNMS(KP250000000, T4d, T4a);
			 T4g = T4e - T4f;
			 T4L = T4f + T4e;
			 T4h = T49 + T4g;
			 T4Z = T4M + T4L;
			 T4x = T4g - T49;
			 T4N = T4L - T4M;
		    }
		    {
			 E T44, T4I, T3X, T4H, T3V;
			 T44 = FNMS(KP951056516, T43, KP587785252 * T40);
			 T4I = FMA(KP951056516, T40, KP587785252 * T43);
			 T3V = FNMS(KP250000000, T3U, T3N);
			 T3X = T3V - T3W;
			 T4H = T3W + T3V;
			 T45 = T3X - T44;
			 T4X = T4H - T4I;
			 T4v = T3X + T44;
			 T4J = T4H + T4I;
		    }
		    {
			 E T3C, T4q, T3J, T4r, T3H;
			 T3C = FNMS(KP951056516, T3B, KP587785252 * T3A);
			 T4q = FMA(KP951056516, T3A, KP587785252 * T3B);
			 T3H = FNMS(KP250000000, T3G, T3D);
			 T3J = T3H - T3I;
			 T4r = T3I + T3H;
			 T3K = T3C + T3J;
			 T56 = T4r - T4q;
			 T4s = T4q + T4r;
			 T4U = T3J - T3C;
		    }
		    {
			 E T4O, T4Q, T4G, T4K;
			 T4G = W[20];
			 T4K = W[21];
			 T4O = FMA(T4G, T4J, T4K * T4N);
			 T4Q = FNMS(T4K, T4J, T4G * T4N);
			 Rp[WS(rs, 5)] = T4F - T4O;
			 Ip[WS(rs, 5)] = T4P + T4Q;
			 Rm[WS(rs, 5)] = T4F + T4O;
			 Im[WS(rs, 5)] = T4Q - T4P;
		    }
		    {
			 E T57, T5d, T53, T55;
			 T53 = W[26];
			 T55 = W[27];
			 T57 = FNMS(T55, T56, T53 * T54);
			 T5d = FMA(T55, T54, T53 * T56);
			 Rp[WS(rs, 7)] = T57 - T5c;
			 Ip[WS(rs, 7)] = T5d + T5e;
			 Rm[WS(rs, 7)] = T57 + T5c;
			 Im[WS(rs, 7)] = T5e - T5d;
		    }
		    {
			 E T3L, T4j, T4i, T4k;
			 {
			      E T3f, T3z, T3M, T46;
			      T3f = W[2];
			      T3z = W[3];
			      T3L = FNMS(T3z, T3K, T3f * T3y);
			      T4j = FMA(T3z, T3y, T3f * T3K);
			      T3M = W[4];
			      T46 = W[5];
			      T4i = FMA(T3M, T45, T46 * T4h);
			      T4k = FNMS(T46, T45, T3M * T4h);
			 }
			 Rp[WS(rs, 1)] = T3L - T4i;
			 Ip[WS(rs, 1)] = T4j + T4k;
			 Rm[WS(rs, 1)] = T3L + T4i;
			 Im[WS(rs, 1)] = T4k - T4j;
		    }
		    {
			 E T4t, T4z, T4y, T4A;
			 {
			      E T4l, T4p, T4u, T4w;
			      T4l = W[10];
			      T4p = W[11];
			      T4t = FNMS(T4p, T4s, T4l * T4o);
			      T4z = FMA(T4p, T4o, T4l * T4s);
			      T4u = W[12];
			      T4w = W[13];
			      T4y = FMA(T4u, T4v, T4w * T4x);
			      T4A = FNMS(T4w, T4v, T4u * T4x);
			 }
			 Rp[WS(rs, 3)] = T4t - T4y;
			 Ip[WS(rs, 3)] = T4z + T4A;
			 Rm[WS(rs, 3)] = T4t + T4y;
			 Im[WS(rs, 3)] = T4A - T4z;
		    }
		    {
			 E T4V, T51, T50, T52;
			 {
			      E T4R, T4T, T4W, T4Y;
			      T4R = W[34];
			      T4T = W[35];
			      T4V = FNMS(T4T, T4U, T4R * T4S);
			      T51 = FMA(T4T, T4S, T4R * T4U);
			      T4W = W[36];
			      T4Y = W[37];
			      T50 = FMA(T4W, T4X, T4Y * T4Z);
			      T52 = FNMS(T4Y, T4X, T4W * T4Z);
			 }
			 Rp[WS(rs, 9)] = T4V - T50;
			 Ip[WS(rs, 9)] = T51 + T52;
			 Rm[WS(rs, 9)] = T4V + T50;
			 Im[WS(rs, 9)] = T52 - T51;
		    }
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     { TW_FULL, 1, 20 },
     { TW_NEXT, 1, 0 }
};

static const hc2c_desc desc = { 20, "hc2cbdft_20", twinstr, &GENUS, { 224, 62, 62, 0 } };

void X(codelet_hc2cbdft_20) (planner *p) {
     X(khc2c_register) (p, hc2cbdft_20, &desc, HC2C_VIA_DFT);
}
#endif
