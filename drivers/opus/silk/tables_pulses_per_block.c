/***********************************************************************
Copyright (c) 2006-2011, Skype Limited. All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
- Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
- Neither the name of Internet Society, IETF or IETF Trust, nor the
names of specific contributors, may be used to endorse or promote
products derived from this software without specific prior written
permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/

#ifdef OPUS_HAVE_CONFIG_H
#include "opus_config.h"
#endif

#include "tables.h"

const opus_uint8 silk_max_pulses_table[ 4 ] = {
         8,     10,     12,     16
};

const opus_uint8 silk_pulses_per_block_iCDF[ 10 ][ 18 ] = {
{
       125,     51,     26,     18,     15,     12,     11,     10,
         9,      8,      7,      6,      5,      4,      3,      2,
         1,      0
},
{
       198,    105,     45,     22,     15,     12,     11,     10,
         9,      8,      7,      6,      5,      4,      3,      2,
         1,      0
},
{
       213,    162,    116,     83,     59,     43,     32,     24,
        18,     15,     12,      9,      7,      6,      5,      3,
         2,      0
},
{
       239,    187,    116,     59,     28,     16,     11,     10,
         9,      8,      7,      6,      5,      4,      3,      2,
         1,      0
},
{
       250,    229,    188,    135,     86,     51,     30,     19,
        13,     10,      8,      6,      5,      4,      3,      2,
         1,      0
},
{
       249,    235,    213,    185,    156,    128,    103,     83,
        66,     53,     42,     33,     26,     21,     17,     13,
        10,      0
},
{
       254,    249,    235,    206,    164,    118,     77,     46,
        27,     16,     10,      7,      5,      4,      3,      2,
         1,      0
},
{
       255,    253,    249,    239,    220,    191,    156,    119,
        85,     57,     37,     23,     15,     10,      6,      4,
         2,      0
},
{
       255,    253,    251,    246,    237,    223,    203,    179,
       152,    124,     98,     75,     55,     40,     29,     21,
        15,      0
},
{
       255,    254,    253,    247,    220,    162,    106,     67,
        42,     28,     18,     12,      9,      6,      4,      3,
         2,      0
}
};

const opus_uint8 silk_pulses_per_block_BITS_Q5[ 9 ][ 18 ] = {
{
        31,     57,    107,    160,    205,    205,    255,    255,
       255,    255,    255,    255,    255,    255,    255,    255,
       255,    255
},
{
        69,     47,     67,    111,    166,    205,    255,    255,
       255,    255,    255,    255,    255,    255,    255,    255,
       255,    255
},
{
        82,     74,     79,     95,    109,    128,    145,    160,
       173,    205,    205,    205,    224,    255,    255,    224,
       255,    224
},
{
       125,     74,     59,     69,     97,    141,    182,    255,
       255,    255,    255,    255,    255,    255,    255,    255,
       255,    255
},
{
       173,    115,     85,     73,     76,     92,    115,    145,
       173,    205,    224,    224,    255,    255,    255,    255,
       255,    255
},
{
       166,    134,    113,    102,    101,    102,    107,    118,
       125,    138,    145,    155,    166,    182,    192,    192,
       205,    150
},
{
       224,    182,    134,    101,     83,     79,     85,     97,
       120,    145,    173,    205,    224,    255,    255,    255,
       255,    255
},
{
       255,    224,    192,    150,    120,    101,     92,     89,
        93,    102,    118,    134,    160,    182,    192,    224,
       224,    224
},
{
       255,    224,    224,    182,    155,    134,    118,    109,
       104,    102,    106,    111,    118,    131,    145,    160,
       173,    131
}
};

const opus_uint8 silk_rate_levels_iCDF[ 2 ][ 9 ] =
{
{
       241,    190,    178,    132,     87,     74,     41,     14,
         0
},
{
       223,    193,    157,    140,    106,     57,     39,     18,
         0
}
};

const opus_uint8 silk_rate_levels_BITS_Q5[ 2 ][ 9 ] =
{
{
       131,     74,    141,     79,     80,    138,     95,    104,
       134
},
{
        95,     99,     91,    125,     93,     76,    123,    115,
       123
}
};

const opus_uint8 silk_shell_code_table0[ 152 ] = {
       128,      0,    214,     42,      0,    235,    128,     21,
         0,    244,    184,     72,     11,      0,    248,    214,
       128,     42,      7,      0,    248,    225,    170,     80,
        25,      5,      0,    251,    236,    198,    126,     54,
        18,      3,      0,    250,    238,    211,    159,     82,
        35,     15,      5,      0,    250,    231,    203,    168,
       128,     88,     53,     25,      6,      0,    252,    238,
       216,    185,    148,    108,     71,     40,     18,      4,
         0,    253,    243,    225,    199,    166,    128,     90,
        57,     31,     13,      3,      0,    254,    246,    233,
       212,    183,    147,    109,     73,     44,     23,     10,
         2,      0,    255,    250,    240,    223,    198,    166,
       128,     90,     58,     33,     16,      6,      1,      0,
       255,    251,    244,    231,    210,    181,    146,    110,
        75,     46,     25,     12,      5,      1,      0,    255,
       253,    248,    238,    221,    196,    164,    128,     92,
        60,     35,     18,      8,      3,      1,      0,    255,
       253,    249,    242,    229,    208,    180,    146,    110,
        76,     48,     27,     14,      7,      3,      1,      0
};

const opus_uint8 silk_shell_code_table1[ 152 ] = {
       129,      0,    207,     50,      0,    236,    129,     20,
         0,    245,    185,     72,     10,      0,    249,    213,
       129,     42,      6,      0,    250,    226,    169,     87,
        27,      4,      0,    251,    233,    194,    130,     62,
        20,      4,      0,    250,    236,    207,    160,     99,
        47,     17,      3,      0,    255,    240,    217,    182,
       131,     81,     41,     11,      1,      0,    255,    254,
       233,    201,    159,    107,     61,     20,      2,      1,
         0,    255,    249,    233,    206,    170,    128,     86,
        50,     23,      7,      1,      0,    255,    250,    238,
       217,    186,    148,    108,     70,     39,     18,      6,
         1,      0,    255,    252,    243,    226,    200,    166,
       128,     90,     56,     30,     13,      4,      1,      0,
       255,    252,    245,    231,    209,    180,    146,    110,
        76,     47,     25,     11,      4,      1,      0,    255,
       253,    248,    237,    219,    194,    163,    128,     93,
        62,     37,     19,      8,      3,      1,      0,    255,
       254,    250,    241,    226,    205,    177,    145,    111,
        79,     51,     30,     15,      6,      2,      1,      0
};

const opus_uint8 silk_shell_code_table2[ 152 ] = {
       129,      0,    203,     54,      0,    234,    129,     23,
         0,    245,    184,     73,     10,      0,    250,    215,
       129,     41,      5,      0,    252,    232,    173,     86,
        24,      3,      0,    253,    240,    200,    129,     56,
        15,      2,      0,    253,    244,    217,    164,     94,
        38,     10,      1,      0,    253,    245,    226,    189,
       132,     71,     27,      7,      1,      0,    253,    246,
       231,    203,    159,    105,     56,     23,      6,      1,
         0,    255,    248,    235,    213,    179,    133,     85,
        47,     19,      5,      1,      0,    255,    254,    243,
       221,    194,    159,    117,     70,     37,     12,      2,
         1,      0,    255,    254,    248,    234,    208,    171,
       128,     85,     48,     22,      8,      2,      1,      0,
       255,    254,    250,    240,    220,    189,    149,    107,
        67,     36,     16,      6,      2,      1,      0,    255,
       254,    251,    243,    227,    201,    166,    128,     90,
        55,     29,     13,      5,      2,      1,      0,    255,
       254,    252,    246,    234,    213,    183,    147,    109,
        73,     43,     22,     10,      4,      2,      1,      0
};

const opus_uint8 silk_shell_code_table3[ 152 ] = {
       130,      0,    200,     58,      0,    231,    130,     26,
         0,    244,    184,     76,     12,      0,    249,    214,
       130,     43,      6,      0,    252,    232,    173,     87,
        24,      3,      0,    253,    241,    203,    131,     56,
        14,      2,      0,    254,    246,    221,    167,     94,
        35,      8,      1,      0,    254,    249,    232,    193,
       130,     65,     23,      5,      1,      0,    255,    251,
       239,    211,    162,     99,     45,     15,      4,      1,
         0,    255,    251,    243,    223,    186,    131,     74,
        33,     11,      3,      1,      0,    255,    252,    245,
       230,    202,    158,    105,     57,     24,      8,      2,
         1,      0,    255,    253,    247,    235,    214,    179,
       132,     84,     44,     19,      7,      2,      1,      0,
       255,    254,    250,    240,    223,    196,    159,    112,
        69,     36,     15,      6,      2,      1,      0,    255,
       254,    253,    245,    231,    209,    176,    136,     93,
        55,     27,     11,      3,      2,      1,      0,    255,
       254,    253,    252,    239,    221,    194,    158,    117,
        76,     42,     18,      4,      3,      2,      1,      0
};

const opus_uint8 silk_shell_code_table_offsets[ 17 ] = {
         0,      0,      2,      5,      9,     14,     20,     27,
        35,     44,     54,     65,     77,     90,    104,    119,
       135
};

const opus_uint8 silk_sign_iCDF[ 42 ] = {
       254,     49,     67,     77,     82,     93,     99,
       198,     11,     18,     24,     31,     36,     45,
       255,     46,     66,     78,     87,     94,    104,
       208,     14,     21,     32,     42,     51,     66,
       255,     94,    104,    109,    112,    115,    118,
       248,     53,     69,     80,     88,     95,    102
};
