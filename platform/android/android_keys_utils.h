/*************************************************************************/
/*  android_keys_utils.h                                                 */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2020 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2020 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#ifndef ANDROID_KEYS_UTILS_H
#define ANDROID_KEYS_UTILS_H

#include <core/os/keyboard.h>

/*
 * Android Key codes.
 */
enum {
	AKEYCODE_UNKNOWN = 0,
	AKEYCODE_SOFT_LEFT = 1,
	AKEYCODE_SOFT_RIGHT = 2,
	AKEYCODE_HOME = 3,
	AKEYCODE_BACK = 4,
	AKEYCODE_CALL = 5,
	AKEYCODE_ENDCALL = 6,
	AKEYCODE_0 = 7,
	AKEYCODE_1 = 8,
	AKEYCODE_2 = 9,
	AKEYCODE_3 = 10,
	AKEYCODE_4 = 11,
	AKEYCODE_5 = 12,
	AKEYCODE_6 = 13,
	AKEYCODE_7 = 14,
	AKEYCODE_8 = 15,
	AKEYCODE_9 = 16,
	AKEYCODE_STAR = 17,
	AKEYCODE_POUND = 18,
	AKEYCODE_DPAD_UP = 19,
	AKEYCODE_DPAD_DOWN = 20,
	AKEYCODE_DPAD_LEFT = 21,
	AKEYCODE_DPAD_RIGHT = 22,
	AKEYCODE_DPAD_CENTER = 23,
	AKEYCODE_VOLUME_UP = 24,
	AKEYCODE_VOLUME_DOWN = 25,
	AKEYCODE_POWER = 26,
	AKEYCODE_CAMERA = 27,
	AKEYCODE_CLEAR = 28,
	AKEYCODE_A = 29,
	AKEYCODE_B = 30,
	AKEYCODE_C = 31,
	AKEYCODE_D = 32,
	AKEYCODE_E = 33,
	AKEYCODE_F = 34,
	AKEYCODE_G = 35,
	AKEYCODE_H = 36,
	AKEYCODE_I = 37,
	AKEYCODE_J = 38,
	AKEYCODE_K = 39,
	AKEYCODE_L = 40,
	AKEYCODE_M = 41,
	AKEYCODE_N = 42,
	AKEYCODE_O = 43,
	AKEYCODE_P = 44,
	AKEYCODE_Q = 45,
	AKEYCODE_R = 46,
	AKEYCODE_S = 47,
	AKEYCODE_T = 48,
	AKEYCODE_U = 49,
	AKEYCODE_V = 50,
	AKEYCODE_W = 51,
	AKEYCODE_X = 52,
	AKEYCODE_Y = 53,
	AKEYCODE_Z = 54,
	AKEYCODE_COMMA = 55,
	AKEYCODE_PERIOD = 56,
	AKEYCODE_ALT_LEFT = 57,
	AKEYCODE_ALT_RIGHT = 58,
	AKEYCODE_SHIFT_LEFT = 59,
	AKEYCODE_SHIFT_RIGHT = 60,
	AKEYCODE_TAB = 61,
	AKEYCODE_SPACE = 62,
	AKEYCODE_SYM = 63,
	AKEYCODE_EXPLORER = 64,
	AKEYCODE_ENVELOPE = 65,
	AKEYCODE_ENTER = 66,
	AKEYCODE_DEL = 67,
	AKEYCODE_GRAVE = 68,
	AKEYCODE_MINUS = 69,
	AKEYCODE_EQUALS = 70,
	AKEYCODE_LEFT_BRACKET = 71,
	AKEYCODE_RIGHT_BRACKET = 72,
	AKEYCODE_BACKSLASH = 73,
	AKEYCODE_SEMICOLON = 74,
	AKEYCODE_APOSTROPHE = 75,
	AKEYCODE_SLASH = 76,
	AKEYCODE_AT = 77,
	AKEYCODE_NUM = 78,
	AKEYCODE_HEADSETHOOK = 79,
	AKEYCODE_FOCUS = 80, // *Camera* focus
	AKEYCODE_PLUS = 81,
	AKEYCODE_MENU = 82,
	AKEYCODE_NOTIFICATION = 83,
	AKEYCODE_SEARCH = 84,
	AKEYCODE_MEDIA_PLAY_PAUSE = 85,
	AKEYCODE_MEDIA_STOP = 86,
	AKEYCODE_MEDIA_NEXT = 87,
	AKEYCODE_MEDIA_PREVIOUS = 88,
	AKEYCODE_MEDIA_REWIND = 89,
	AKEYCODE_MEDIA_FAST_FORWARD = 90,
	AKEYCODE_MUTE = 91,
	AKEYCODE_PAGE_UP = 92,
	AKEYCODE_PAGE_DOWN = 93,
	AKEYCODE_PICTSYMBOLS = 94,
	AKEYCODE_SWITCH_CHARSET = 95,
	AKEYCODE_BUTTON_A = 96,
	AKEYCODE_BUTTON_B = 97,
	AKEYCODE_BUTTON_C = 98,
	AKEYCODE_BUTTON_X = 99,
	AKEYCODE_BUTTON_Y = 100,
	AKEYCODE_BUTTON_Z = 101,
	AKEYCODE_BUTTON_L1 = 102,
	AKEYCODE_BUTTON_R1 = 103,
	AKEYCODE_BUTTON_L2 = 104,
	AKEYCODE_BUTTON_R2 = 105,
	AKEYCODE_BUTTON_THUMBL = 106,
	AKEYCODE_BUTTON_THUMBR = 107,
	AKEYCODE_BUTTON_START = 108,
	AKEYCODE_BUTTON_SELECT = 109,
	AKEYCODE_BUTTON_MODE = 110,

	// NOTE: If you add a new keycode here you must also add it to several other files.
	//       Refer to frameworks/base/core/java/android/view/KeyEvent.java for the full list.
};

struct _WinTranslatePair {

	unsigned int keysym;
	unsigned int keycode;
};

static _WinTranslatePair _ak_to_keycode[] = {
	{ KEY_TAB, AKEYCODE_TAB },
	{ KEY_ENTER, AKEYCODE_ENTER },
	{ KEY_SHIFT, AKEYCODE_SHIFT_LEFT },
	{ KEY_SHIFT, AKEYCODE_SHIFT_RIGHT },
	{ KEY_ALT, AKEYCODE_ALT_LEFT },
	{ KEY_ALT, AKEYCODE_ALT_RIGHT },
	{ KEY_MENU, AKEYCODE_MENU },
	{ KEY_PAUSE, AKEYCODE_MEDIA_PLAY_PAUSE },
	{ KEY_ESCAPE, AKEYCODE_BACK },
	{ KEY_SPACE, AKEYCODE_SPACE },
	{ KEY_PAGEUP, AKEYCODE_PAGE_UP },
	{ KEY_PAGEDOWN, AKEYCODE_PAGE_DOWN },
	{ KEY_HOME, AKEYCODE_HOME }, //(0x24)
	{ KEY_LEFT, AKEYCODE_DPAD_LEFT },
	{ KEY_UP, AKEYCODE_DPAD_UP },
	{ KEY_RIGHT, AKEYCODE_DPAD_RIGHT },
	{ KEY_DOWN, AKEYCODE_DPAD_DOWN },
	{ KEY_PERIODCENTERED, AKEYCODE_DPAD_CENTER },
	{ KEY_BACKSPACE, AKEYCODE_DEL },
	{ KEY_0, AKEYCODE_0 }, ////0 key
	{ KEY_1, AKEYCODE_1 }, ////1 key
	{ KEY_2, AKEYCODE_2 }, ////2 key
	{ KEY_3, AKEYCODE_3 }, ////3 key
	{ KEY_4, AKEYCODE_4 }, ////4 key
	{ KEY_5, AKEYCODE_5 }, ////5 key
	{ KEY_6, AKEYCODE_6 }, ////6 key
	{ KEY_7, AKEYCODE_7 }, ////7 key
	{ KEY_8, AKEYCODE_8 }, ////8 key
	{ KEY_9, AKEYCODE_9 }, ////9 key
	{ KEY_A, AKEYCODE_A }, ////A key
	{ KEY_B, AKEYCODE_B }, ////B key
	{ KEY_C, AKEYCODE_C }, ////C key
	{ KEY_D, AKEYCODE_D }, ////D key
	{ KEY_E, AKEYCODE_E }, ////E key
	{ KEY_F, AKEYCODE_F }, ////F key
	{ KEY_G, AKEYCODE_G }, ////G key
	{ KEY_H, AKEYCODE_H }, ////H key
	{ KEY_I, AKEYCODE_I }, ////I key
	{ KEY_J, AKEYCODE_J }, ////J key
	{ KEY_K, AKEYCODE_K }, ////K key
	{ KEY_L, AKEYCODE_L }, ////L key
	{ KEY_M, AKEYCODE_M }, ////M key
	{ KEY_N, AKEYCODE_N }, ////N key
	{ KEY_O, AKEYCODE_O }, ////O key
	{ KEY_P, AKEYCODE_P }, ////P key
	{ KEY_Q, AKEYCODE_Q }, ////Q key
	{ KEY_R, AKEYCODE_R }, ////R key
	{ KEY_S, AKEYCODE_S }, ////S key
	{ KEY_T, AKEYCODE_T }, ////T key
	{ KEY_U, AKEYCODE_U }, ////U key
	{ KEY_V, AKEYCODE_V }, ////V key
	{ KEY_W, AKEYCODE_W }, ////W key
	{ KEY_X, AKEYCODE_X }, ////X key
	{ KEY_Y, AKEYCODE_Y }, ////Y key
	{ KEY_Z, AKEYCODE_Z }, ////Z key
	{ KEY_HOMEPAGE, AKEYCODE_EXPLORER },
	{ KEY_LAUNCH0, AKEYCODE_BUTTON_A },
	{ KEY_LAUNCH1, AKEYCODE_BUTTON_B },
	{ KEY_LAUNCH2, AKEYCODE_BUTTON_C },
	{ KEY_LAUNCH3, AKEYCODE_BUTTON_X },
	{ KEY_LAUNCH4, AKEYCODE_BUTTON_Y },
	{ KEY_LAUNCH5, AKEYCODE_BUTTON_Z },
	{ KEY_LAUNCH6, AKEYCODE_BUTTON_L1 },
	{ KEY_LAUNCH7, AKEYCODE_BUTTON_R1 },
	{ KEY_LAUNCH8, AKEYCODE_BUTTON_L2 },
	{ KEY_LAUNCH9, AKEYCODE_BUTTON_R2 },
	{ KEY_LAUNCHA, AKEYCODE_BUTTON_THUMBL },
	{ KEY_LAUNCHB, AKEYCODE_BUTTON_THUMBR },
	{ KEY_LAUNCHC, AKEYCODE_BUTTON_START },
	{ KEY_LAUNCHD, AKEYCODE_BUTTON_SELECT },
	{ KEY_LAUNCHE, AKEYCODE_BUTTON_MODE },
	{ KEY_VOLUMEMUTE, AKEYCODE_MUTE },
	{ KEY_VOLUMEDOWN, AKEYCODE_VOLUME_DOWN },
	{ KEY_VOLUMEUP, AKEYCODE_VOLUME_UP },
	{ KEY_BACK, AKEYCODE_MEDIA_REWIND },
	{ KEY_FORWARD, AKEYCODE_MEDIA_FAST_FORWARD },
	{ KEY_MEDIANEXT, AKEYCODE_MEDIA_NEXT },
	{ KEY_MEDIAPREVIOUS, AKEYCODE_MEDIA_PREVIOUS },
	{ KEY_MEDIASTOP, AKEYCODE_MEDIA_STOP },
	{ KEY_PLUS, AKEYCODE_PLUS },
	{ KEY_EQUAL, AKEYCODE_EQUALS }, // the '+' key
	{ KEY_COMMA, AKEYCODE_COMMA }, // the ',' key
	{ KEY_MINUS, AKEYCODE_MINUS }, // the '-' key
	{ KEY_SLASH, AKEYCODE_SLASH }, // the '/?' key
	{ KEY_BACKSLASH, AKEYCODE_BACKSLASH },
	{ KEY_BRACKETLEFT, AKEYCODE_LEFT_BRACKET },
	{ KEY_BRACKETRIGHT, AKEYCODE_RIGHT_BRACKET },
	{ KEY_UNKNOWN, 0 }
};
/*
TODO: map these android key:
    AKEYCODE_SOFT_LEFT       = 1,
    AKEYCODE_SOFT_RIGHT      = 2,
    AKEYCODE_CALL            = 5,
    AKEYCODE_ENDCALL         = 6,
    AKEYCODE_STAR            = 17,
    AKEYCODE_POUND           = 18,
    AKEYCODE_POWER           = 26,
    AKEYCODE_CAMERA          = 27,
    AKEYCODE_CLEAR           = 28,
    AKEYCODE_SYM             = 63,
    AKEYCODE_ENVELOPE        = 65,
    AKEYCODE_GRAVE           = 68,
    AKEYCODE_SEMICOLON       = 74,
    AKEYCODE_APOSTROPHE      = 75,
    AKEYCODE_AT              = 77,
    AKEYCODE_NUM             = 78,
    AKEYCODE_HEADSETHOOK     = 79,
    AKEYCODE_FOCUS           = 80,   // *Camera* focus
    AKEYCODE_NOTIFICATION    = 83,
    AKEYCODE_SEARCH          = 84,
    AKEYCODE_PICTSYMBOLS     = 94,
    AKEYCODE_SWITCH_CHARSET  = 95,
*/

unsigned int android_get_keysym(unsigned int p_code);

#endif // ANDROID_KEYS_UTILS_H
