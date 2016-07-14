// Copyright 2011 Google Inc. All Rights Reserved.
//
// Use of this source code is governed by a BSD-style license
// that can be found in the COPYING file in the root of the source
// tree. An additional intellectual property rights grant can be found
// in the file PATENTS. All contributing project authors may
// be found in the AUTHORS file in the root of the source tree.
// -----------------------------------------------------------------------------
//
// Selecting filter level
//
// Author: somnath@google.com (Somnath Banerjee)

#include <assert.h>
#include "./vp8enci.h"
#include "../dsp/dsp.h"

// This table gives, for a given sharpness, the filtering strength to be
// used (at least) in order to filter a given edge step delta.
// This is constructed by brute force inspection: for all delta, we iterate
// over all possible filtering strength / thresh until needs_filter() returns
// true.
#define MAX_DELTA_SIZE 64
static const uint8_t kLevelsFromDelta[8][MAX_DELTA_SIZE] = {
  { 0,   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 },
  { 0,  1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 13, 14, 15, 17, 18,
    20, 21, 23, 24, 26, 27, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42,
    44, 45, 47, 48, 50, 51, 53, 54, 56, 57, 59, 60, 62, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 13, 14, 16, 17, 19,
    20, 22, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43,
    44, 46, 47, 49, 50, 52, 53, 55, 56, 58, 59, 61, 62, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 13, 15, 16, 18, 19,
    21, 22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 40, 42, 43,
    45, 46, 48, 49, 51, 52, 54, 55, 57, 58, 60, 61, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 14, 15, 17, 18, 20,
    21, 23, 24, 26, 27, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42, 44,
    45, 47, 48, 50, 51, 53, 54, 56, 57, 59, 60, 62, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  4,  5,  7,  8,  9, 11, 12, 13, 15, 16, 17, 19, 20,
    22, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43, 44,
    46, 47, 49, 50, 52, 53, 55, 56, 58, 59, 61, 62, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  4,  5,  7,  8,  9, 11, 12, 13, 15, 16, 18, 19, 21,
    22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 40, 42, 43, 45,
    46, 48, 49, 51, 52, 54, 55, 57, 58, 60, 61, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 },
  {  0,  1,  2,  4,  5,  7,  8,  9, 11, 12, 14, 15, 17, 18, 20, 21,
    23, 24, 26, 27, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42, 44, 45,
    47, 48, 50, 51, 53, 54, 56, 57, 59, 60, 62, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 }
};

int VP8FilterStrengthFromDelta(int sharpness, int delta) {
  const int pos = (delta < MAX_DELTA_SIZE) ? delta : MAX_DELTA_SIZE - 1;
  assert(sharpness >= 0 && sharpness <= 7);
  return kLevelsFromDelta[sharpness][pos];
}

//------------------------------------------------------------------------------
// Paragraph 15.4: compute the inner-edge filtering strength

static int GetILevel(int sharpness, int level) {
  if (sharpness > 0) {
    if (sharpness > 4) {
      level >>= 2;
    } else {
      level >>= 1;
    }
    if (level > 9 - sharpness) {
      level = 9 - sharpness;
    }
  }
  if (level < 1) level = 1;
  return level;
}

static void DoFilter(const VP8EncIterator* const it, int level) {
  const VP8Encoder* const enc = it->enc_;
  const int ilevel = GetILevel(enc->config_->filter_sharpness, level);
  const int limit = 2 * level + ilevel;

  uint8_t* const y_dst = it->yuv_out2_ + Y_OFF_ENC;
  uint8_t* const u_dst = it->yuv_out2_ + U_OFF_ENC;
  uint8_t* const v_dst = it->yuv_out2_ + V_OFF_ENC;

  // copy current block to yuv_out2_
  memcpy(y_dst, it->yuv_out_, YUV_SIZE_ENC * sizeof(uint8_t));

  if (enc->filter_hdr_.simple_ == 1) {   // simple
    VP8SimpleHFilter16i(y_dst, BPS, limit);
    VP8SimpleVFilter16i(y_dst, BPS, limit);
  } else {    // complex
    const int hev_thresh = (level >= 40) ? 2 : (level >= 15) ? 1 : 0;
    VP8HFilter16i(y_dst, BPS, limit, ilevel, hev_thresh);
    VP8HFilter8i(u_dst, v_dst, BPS, limit, ilevel, hev_thresh);
    VP8VFilter16i(y_dst, BPS, limit, ilevel, hev_thresh);
    VP8VFilter8i(u_dst, v_dst, BPS, limit, ilevel, hev_thresh);
  }
}

//------------------------------------------------------------------------------
// SSIM metric

static const double kMinValue = 1.e-10;  // minimal threshold

void VP8SSIMAddStats(const VP8DistoStats* const src, VP8DistoStats* const dst) {
  dst->w   += src->w;
  dst->xm  += src->xm;
  dst->ym  += src->ym;
  dst->xxm += src->xxm;
  dst->xym += src->xym;
  dst->yym += src->yym;
}

double VP8SSIMGet(const VP8DistoStats* const stats) {
  const double xmxm = stats->xm * stats->xm;
  const double ymym = stats->ym * stats->ym;
  const double xmym = stats->xm * stats->ym;
  const double w2 = stats->w * stats->w;
  double sxx = stats->xxm * stats->w - xmxm;
  double syy = stats->yym * stats->w - ymym;
  double sxy = stats->xym * stats->w - xmym;
  double C1, C2;
  double fnum;
  double fden;
  // small errors are possible, due to rounding. Clamp to zero.
  if (sxx < 0.) sxx = 0.;
  if (syy < 0.) syy = 0.;
  C1 = 6.5025 * w2;
  C2 = 58.5225 * w2;
  fnum = (2 * xmym + C1) * (2 * sxy + C2);
  fden = (xmxm + ymym + C1) * (sxx + syy + C2);
  return (fden != 0.) ? fnum / fden : kMinValue;
}

double VP8SSIMGetSquaredError(const VP8DistoStats* const s) {
  if (s->w > 0.) {
    const double iw2 = 1. / (s->w * s->w);
    const double sxx = s->xxm * s->w - s->xm * s->xm;
    const double syy = s->yym * s->w - s->ym * s->ym;
    const double sxy = s->xym * s->w - s->xm * s->ym;
    const double SSE = iw2 * (sxx + syy - 2. * sxy);
    if (SSE > kMinValue) return SSE;
  }
  return kMinValue;
}

#define LIMIT(A, M)  ((A) > (M) ? (M) : (A))
static void VP8SSIMAccumulateRow(const uint8_t* src1, int stride1,
                                 const uint8_t* src2, int stride2,
                                 int y, int W, int H,
                                 VP8DistoStats* const stats) {
  int x = 0;
  const int w0 = LIMIT(VP8_SSIM_KERNEL, W);
  for (x = 0; x < w0; ++x) {
    VP8SSIMAccumulateClipped(src1, stride1, src2, stride2, x, y, W, H, stats);
  }
  for (; x <= W - 8 + VP8_SSIM_KERNEL; ++x) {
    VP8SSIMAccumulate(
        src1 + (y - VP8_SSIM_KERNEL) * stride1 + (x - VP8_SSIM_KERNEL), stride1,
        src2 + (y - VP8_SSIM_KERNEL) * stride2 + (x - VP8_SSIM_KERNEL), stride2,
        stats);
  }
  for (; x < W; ++x) {
    VP8SSIMAccumulateClipped(src1, stride1, src2, stride2, x, y, W, H, stats);
  }
}

void VP8SSIMAccumulatePlane(const uint8_t* src1, int stride1,
                            const uint8_t* src2, int stride2,
                            int W, int H, VP8DistoStats* const stats) {
  int x, y;
  const int h0 = LIMIT(VP8_SSIM_KERNEL, H);
  const int h1 = LIMIT(VP8_SSIM_KERNEL, H - VP8_SSIM_KERNEL);
  for (y = 0; y < h0; ++y) {
    for (x = 0; x < W; ++x) {
      VP8SSIMAccumulateClipped(src1, stride1, src2, stride2, x, y, W, H, stats);
    }
  }
  for (; y < h1; ++y) {
    VP8SSIMAccumulateRow(src1, stride1, src2, stride2, y, W, H, stats);
  }
  for (; y < H; ++y) {
    for (x = 0; x < W; ++x) {
      VP8SSIMAccumulateClipped(src1, stride1, src2, stride2, x, y, W, H, stats);
    }
  }
}
#undef LIMIT

static double GetMBSSIM(const uint8_t* yuv1, const uint8_t* yuv2) {
  int x, y;
  VP8DistoStats s = { .0, .0, .0, .0, .0, .0 };

  // compute SSIM in a 10 x 10 window
  for (y = VP8_SSIM_KERNEL; y < 16 - VP8_SSIM_KERNEL; y++) {
    for (x = VP8_SSIM_KERNEL; x < 16 - VP8_SSIM_KERNEL; x++) {
      VP8SSIMAccumulateClipped(yuv1 + Y_OFF_ENC, BPS, yuv2 + Y_OFF_ENC, BPS,
                               x, y, 16, 16, &s);
    }
  }
  for (x = 1; x < 7; x++) {
    for (y = 1; y < 7; y++) {
      VP8SSIMAccumulateClipped(yuv1 + U_OFF_ENC, BPS, yuv2 + U_OFF_ENC, BPS,
                               x, y, 8, 8, &s);
      VP8SSIMAccumulateClipped(yuv1 + V_OFF_ENC, BPS, yuv2 + V_OFF_ENC, BPS,
                               x, y, 8, 8, &s);
    }
  }
  return VP8SSIMGet(&s);
}

//------------------------------------------------------------------------------
// Exposed APIs: Encoder should call the following 3 functions to adjust
// loop filter strength

void VP8InitFilter(VP8EncIterator* const it) {
  if (it->lf_stats_ != NULL) {
    int s, i;
    for (s = 0; s < NUM_MB_SEGMENTS; s++) {
      for (i = 0; i < MAX_LF_LEVELS; i++) {
        (*it->lf_stats_)[s][i] = 0;
      }
    }
    VP8SSIMDspInit();
  }
}

void VP8StoreFilterStats(VP8EncIterator* const it) {
  int d;
  VP8Encoder* const enc = it->enc_;
  const int s = it->mb_->segment_;
  const int level0 = enc->dqm_[s].fstrength_;

  // explore +/-quant range of values around level0
  const int delta_min = -enc->dqm_[s].quant_;
  const int delta_max = enc->dqm_[s].quant_;
  const int step_size = (delta_max - delta_min >= 4) ? 4 : 1;

  if (it->lf_stats_ == NULL) return;

  // NOTE: Currently we are applying filter only across the sublock edges
  // There are two reasons for that.
  // 1. Applying filter on macro block edges will change the pixels in
  // the left and top macro blocks. That will be hard to restore
  // 2. Macro Blocks on the bottom and right are not yet compressed. So we
  // cannot apply filter on the right and bottom macro block edges.
  if (it->mb_->type_ == 1 && it->mb_->skip_) return;

  // Always try filter level  zero
  (*it->lf_stats_)[s][0] += GetMBSSIM(it->yuv_in_, it->yuv_out_);

  for (d = delta_min; d <= delta_max; d += step_size) {
    const int level = level0 + d;
    if (level <= 0 || level >= MAX_LF_LEVELS) {
      continue;
    }
    DoFilter(it, level);
    (*it->lf_stats_)[s][level] += GetMBSSIM(it->yuv_in_, it->yuv_out2_);
  }
}

void VP8AdjustFilterStrength(VP8EncIterator* const it) {
  VP8Encoder* const enc = it->enc_;
  if (it->lf_stats_ != NULL) {
    int s;
    for (s = 0; s < NUM_MB_SEGMENTS; s++) {
      int i, best_level = 0;
      // Improvement over filter level 0 should be at least 1e-5 (relatively)
      double best_v = 1.00001 * (*it->lf_stats_)[s][0];
      for (i = 1; i < MAX_LF_LEVELS; i++) {
        const double v = (*it->lf_stats_)[s][i];
        if (v > best_v) {
          best_v = v;
          best_level = i;
        }
      }
      enc->dqm_[s].fstrength_ = best_level;
    }
  } else if (enc->config_->filter_strength > 0) {
    int max_level = 0;
    int s;
    for (s = 0; s < NUM_MB_SEGMENTS; s++) {
      VP8SegmentInfo* const dqm = &enc->dqm_[s];
      // this '>> 3' accounts for some inverse WHT scaling
      const int delta = (dqm->max_edge_ * dqm->y2_.q_[1]) >> 3;
      const int level =
          VP8FilterStrengthFromDelta(enc->filter_hdr_.sharpness_, delta);
      if (level > dqm->fstrength_) {
        dqm->fstrength_ = level;
      }
      if (max_level < dqm->fstrength_) {
        max_level = dqm->fstrength_;
      }
    }
    enc->filter_hdr_.level_ = max_level;
  }
}

// -----------------------------------------------------------------------------
