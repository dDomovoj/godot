/* Copyright (c) 2007-2008 CSIRO
   Copyright (c) 2007-2008 Xiph.Org Foundation
   Written by Jean-Marc Valin */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* This is a simple MDCT implementation that uses a N/4 complex FFT
   to do most of the work. It should be relatively straightforward to
   plug in pretty much and FFT here.

   This replaces the Vorbis FFT (and uses the exact same API), which
   was a bit too messy and that was ending up duplicating code
   (might as well use the same FFT everywhere).

   The algorithm is similar to (and inspired from) Fabrice Bellard's
   MDCT implementation in FFMPEG, but has differences in signs, ordering
   and scaling in many places.
*/
#ifndef __MDCT_MIPSR1_H__
#define __MDCT_MIPSR1_H__

#ifndef SKIP_CONFIG_H
#include "opus/opus_config.h"
#endif

#include "opus/celt/mdct.h"
#include "opus/celt/kiss_fft.h"
#include "opus/celt/_kiss_fft_guts.h"
#include <math.h>
#include "opus/celt/os_support.h"
#include "opus/celt/mathops.h"
#include "opus/celt/stack_alloc.h"

/* Forward MDCT trashes the input array */
#define OVERRIDE_clt_mdct_forward
void clt_mdct_forward(const mdct_lookup *l, kiss_fft_scalar *in, kiss_fft_scalar * OPUS_RESTRICT out,
      const opus_val16 *window, int overlap, int shift, int stride, int arch)
{
   int i;
   int N, N2, N4;
   VARDECL(kiss_fft_scalar, f);
   VARDECL(kiss_fft_cpx, f2);
   const kiss_fft_state *st = l->kfft[shift];
   const kiss_twiddle_scalar *trig;
   opus_val16 scale;
#ifdef OPUS_FIXED_POINT
   /* Allows us to scale with MULT16_32_Q16(), which is faster than
      MULT16_32_Q15() on ARM. */
   int scale_shift = st->scale_shift-1;
#endif

    (void)arch;

   SAVE_STACK;
   scale = st->scale;

   N = l->n;
   trig = l->trig;
   for (i=0;i<shift;i++)
   {
      N >>= 1;
      trig += N;
   }
   N2 = N>>1;
   N4 = N>>2;

   ALLOC(f, N2, kiss_fft_scalar);
   ALLOC(f2, N4, kiss_fft_cpx);

   /* Consider the input to be composed of four blocks: [a, b, c, d] */
   /* Window, shuffle, fold */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_scalar * OPUS_RESTRICT xp1 = in+(overlap>>1);
      const kiss_fft_scalar * OPUS_RESTRICT xp2 = in+N2-1+(overlap>>1);
      kiss_fft_scalar * OPUS_RESTRICT yp = f;
      const opus_val16 * OPUS_RESTRICT wp1 = window+(overlap>>1);
      const opus_val16 * OPUS_RESTRICT wp2 = window+(overlap>>1)-1;
      for(i=0;i<((overlap+3)>>2);i++)
      {
         /* Real part arranged as -d-cR, Imag part arranged as -b+aR*/
          *yp++ = S_MUL_ADD(*wp2, xp1[N2],*wp1,*xp2);
          *yp++ = S_MUL_SUB(*wp1, *xp1,*wp2, xp2[-N2]);
         xp1+=2;
         xp2-=2;
         wp1+=2;
         wp2-=2;
      }
      wp1 = window;
      wp2 = window+overlap-1;
      for(;i<N4-((overlap+3)>>2);i++)
      {
         /* Real part arranged as a-bR, Imag part arranged as -c-dR */
         *yp++ = *xp2;
         *yp++ = *xp1;
         xp1+=2;
         xp2-=2;
      }
      for(;i<N4;i++)
      {
         /* Real part arranged as a-bR, Imag part arranged as -c-dR */
          *yp++ =  S_MUL_SUB(*wp2, *xp2, *wp1, xp1[-N2]);
          *yp++ = S_MUL_ADD(*wp2, *xp1, *wp1, xp2[N2]);
         xp1+=2;
         xp2-=2;
         wp1+=2;
         wp2-=2;
      }
   }
   /* Pre-rotation */
   {
      kiss_fft_scalar * OPUS_RESTRICT yp = f;
      const kiss_twiddle_scalar *t = &trig[0];
      for(i=0;i<N4;i++)
      {
         kiss_fft_cpx yc;
         kiss_twiddle_scalar t0, t1;
         kiss_fft_scalar re, im, yr, yi;
         t0 = t[i];
         t1 = t[N4+i];
         re = *yp++;
         im = *yp++;

         yr = S_MUL_SUB(re,t0,im,t1);
         yi = S_MUL_ADD(im,t0,re,t1);

         yc.r = yr;
         yc.i = yi;
         yc.r = PSHR32(MULT16_32_Q16(scale, yc.r), scale_shift);
         yc.i = PSHR32(MULT16_32_Q16(scale, yc.i), scale_shift);
         f2[st->bitrev[i]] = yc;
      }
   }

   /* N/4 complex FFT, does not downscale anymore */
   opus_fft_impl(st, f2);

   /* Post-rotate */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_cpx * OPUS_RESTRICT fp = f2;
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out;
      kiss_fft_scalar * OPUS_RESTRICT yp2 = out+stride*(N2-1);
      const kiss_twiddle_scalar *t = &trig[0];
      /* Temp pointers to make it really clear to the compiler what we're doing */
      for(i=0;i<N4;i++)
      {
         kiss_fft_scalar yr, yi;
         yr = S_MUL_SUB(fp->i,t[N4+i] , fp->r,t[i]);
         yi = S_MUL_ADD(fp->r,t[N4+i] ,fp->i,t[i]);
         *yp1 = yr;
         *yp2 = yi;
         fp++;
         yp1 += 2*stride;
         yp2 -= 2*stride;
      }
   }
   RESTORE_STACK;
}

#define OVERRIDE_clt_mdct_backward
void clt_mdct_backward(const mdct_lookup *l, kiss_fft_scalar *in, kiss_fft_scalar * OPUS_RESTRICT out,
      const opus_val16 * OPUS_RESTRICT window, int overlap, int shift, int stride, int arch)
{
   int i;
   int N, N2, N4;
   const kiss_twiddle_scalar *trig;

    (void)arch;

   N = l->n;
   trig = l->trig;
   for (i=0;i<shift;i++)
   {
      N >>= 1;
      trig += N;
   }
   N2 = N>>1;
   N4 = N>>2;

   /* Pre-rotate */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_scalar * OPUS_RESTRICT xp1 = in;
      const kiss_fft_scalar * OPUS_RESTRICT xp2 = in+stride*(N2-1);
      kiss_fft_scalar * OPUS_RESTRICT yp = out+(overlap>>1);
      const kiss_twiddle_scalar * OPUS_RESTRICT t = &trig[0];
      const opus_int16 * OPUS_RESTRICT bitrev = l->kfft[shift]->bitrev;
      for(i=0;i<N4;i++)
      {
         int rev;
         kiss_fft_scalar yr, yi;
         rev = *bitrev++;
         yr = S_MUL_ADD(*xp2, t[i] , *xp1, t[N4+i]);
         yi = S_MUL_SUB(*xp1, t[i] , *xp2, t[N4+i]);
         /* We swap real and imag because we use an FFT instead of an IFFT. */
         yp[2*rev+1] = yr;
         yp[2*rev] = yi;
         /* Storing the pre-rotation directly in the bitrev order. */
         xp1+=2*stride;
         xp2-=2*stride;
      }
   }

   opus_fft_impl(l->kfft[shift], (kiss_fft_cpx*)(out+(overlap>>1)));

   /* Post-rotate and de-shuffle from both ends of the buffer at once to make
      it in-place. */
   {
      kiss_fft_scalar * OPUS_RESTRICT yp0 = out+(overlap>>1);
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out+(overlap>>1)+N2-2;
      const kiss_twiddle_scalar *t = &trig[0];
      /* Loop to (N4+1)>>1 to handle odd N4. When N4 is odd, the
         middle pair will be computed twice. */
      for(i=0;i<(N4+1)>>1;i++)
      {
         kiss_fft_scalar re, im, yr, yi;
         kiss_twiddle_scalar t0, t1;
         /* We swap real and imag because we're using an FFT instead of an IFFT. */
         re = yp0[1];
         im = yp0[0];
         t0 = t[i];
         t1 = t[N4+i];
         /* We'd scale up by 2 here, but instead it's done when mixing the windows */
         yr = S_MUL_ADD(re,t0 , im,t1);
         yi = S_MUL_SUB(re,t1 , im,t0);
         /* We swap real and imag because we're using an FFT instead of an IFFT. */
         re = yp1[1];
         im = yp1[0];
         yp0[0] = yr;
         yp1[1] = yi;

         t0 = t[(N4-i-1)];
         t1 = t[(N2-i-1)];
         /* We'd scale up by 2 here, but instead it's done when mixing the windows */
         yr = S_MUL_ADD(re,t0,im,t1);
         yi = S_MUL_SUB(re,t1,im,t0);
         yp1[0] = yr;
         yp0[1] = yi;
         yp0 += 2;
         yp1 -= 2;
      }
   }

   /* Mirror on both sides for TDAC */
   {
      kiss_fft_scalar * OPUS_RESTRICT xp1 = out+overlap-1;
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out;
      const opus_val16 * OPUS_RESTRICT wp1 = window;
      const opus_val16 * OPUS_RESTRICT wp2 = window+overlap-1;

      for(i = 0; i < overlap/2; i++)
      {
         kiss_fft_scalar x1, x2;
         x1 = *xp1;
         x2 = *yp1;
         *yp1++ = MULT16_32_Q15(*wp2, x2) - MULT16_32_Q15(*wp1, x1);
         *xp1-- = MULT16_32_Q15(*wp1, x2) + MULT16_32_Q15(*wp2, x1);
         wp1++;
         wp2--;
      }
   }
}
#endif /* __MDCT_MIPSR1_H__ */
