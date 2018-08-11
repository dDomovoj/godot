/*************************************************************************/
/*  os_javascript.cpp                                                    */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2018 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2018 Godot Engine contributors (cf. AUTHORS.md)    */
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

#include "os_javascript.h"

#include "gles2/rasterizer_gles2.h"
#include "gles3/rasterizer_gles3.h"
#include "io/file_access_buffered_fa.h"
#include "main/main.h"
#include "servers/visual/visual_server_raster.h"
#include "unix/dir_access_unix.h"
#include "unix/file_access_unix.h"

#include <emscripten.h>
#include <stdlib.h>

#include "dom_keys.inc"

#define DOM_BUTTON_LEFT 0
#define DOM_BUTTON_MIDDLE 1
#define DOM_BUTTON_RIGHT 2
#define DOM_BUTTON_XBUTTON1 3
#define DOM_BUTTON_XBUTTON2 4

// Window (canvas)

static void focus_canvas() {

	/* clang-format off */
	EM_ASM(
		Module.canvas.focus();
	);
	/* clang-format on */
}

static bool is_canvas_focused() {

	/* clang-format off */
	return EM_ASM_INT_V(
		return document.activeElement == Module.canvas;
	);
	/* clang-format on */
}

static bool cursor_inside_canvas = true;

EM_BOOL OS_JavaScript::browser_resize_callback(int p_event_type, const EmscriptenUiEvent *p_event, void *p_user_data) {

	// The order of the fullscreen change event and the window size change
	// event varies, even within just one browser, so defer handling.
	get_singleton()->canvas_size_adjustment_requested = true;
	return false;
}

EM_BOOL OS_JavaScript::fullscreen_change_callback(int p_event_type, const EmscriptenFullscreenChangeEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();
	// Empty ID is canvas.
	String target_id = String::utf8(p_event->id);
	if (target_id.empty() || target_id == "canvas") {
		// This event property is the only reliable data on
		// browser fullscreen state.
		os->video_mode.fullscreen = p_event->isFullscreen;
		os->canvas_size_adjustment_requested = true;
	}
	return false;
}

void OS_JavaScript::set_video_mode(const VideoMode &p_video_mode, int p_screen) {

	video_mode = p_video_mode;
}

OS::VideoMode OS_JavaScript::get_video_mode(int p_screen) const {

	return video_mode;
}

Size2 OS_JavaScript::get_screen_size(int p_screen) const {

	EmscriptenFullscreenChangeEvent ev;
	EMSCRIPTEN_RESULT result = emscripten_get_fullscreen_status(&ev);
	ERR_FAIL_COND_V(result != EMSCRIPTEN_RESULT_SUCCESS, Size2());
	return Size2(ev.screenWidth, ev.screenHeight);
}

void OS_JavaScript::set_window_size(const Size2 p_size) {

	windowed_size = p_size;
	if (is_window_fullscreen()) {
		window_maximized = false;
		set_window_fullscreen(false);
	} else if (is_window_maximized()) {
		set_window_maximized(false);
	} else {
		video_mode.width = p_size.x;
		video_mode.height = p_size.y;
		emscripten_set_canvas_size(p_size.x, p_size.y);
	}
}

Size2 OS_JavaScript::get_window_size() const {

	int canvas[3];
	emscripten_get_canvas_size(canvas, canvas + 1, canvas + 2);
	return Size2(canvas[0], canvas[1]);
}

void OS_JavaScript::set_window_maximized(bool p_enabled) {

	window_maximized = p_enabled;
	if (is_window_fullscreen()) {
		set_window_fullscreen(false);
		return;
	}
	// Calling emscripten_enter_soft_fullscreen mutltiple times hides all
	// page elements except the canvas permanently, so track state.
	if (p_enabled && !soft_fullscreen_enabled) {

		EmscriptenFullscreenStrategy strategy;
		strategy.scaleMode = EMSCRIPTEN_FULLSCREEN_SCALE_STRETCH;
		strategy.canvasResolutionScaleMode = EMSCRIPTEN_FULLSCREEN_CANVAS_SCALE_STDDEF;
		strategy.filteringMode = EMSCRIPTEN_FULLSCREEN_FILTERING_DEFAULT;
		strategy.canvasResizedCallback = NULL;
		emscripten_enter_soft_fullscreen(NULL, &strategy);
		soft_fullscreen_enabled = true;
		video_mode.width = get_window_size().width;
		video_mode.height = get_window_size().height;
	} else if (!p_enabled) {

		emscripten_exit_soft_fullscreen();
		soft_fullscreen_enabled = false;
		video_mode.width = windowed_size.width;
		video_mode.height = windowed_size.height;
		emscripten_set_canvas_size(video_mode.width, video_mode.height);
	}
}

bool OS_JavaScript::is_window_maximized() const {

	return window_maximized;
}

void OS_JavaScript::set_window_fullscreen(bool p_enabled) {

	if (p_enabled == is_window_fullscreen()) {
		return;
	}

	// Just request changes here, if successful, canvas is resized in
	// _browser_resize_callback or _fullscreen_change_callback.
	EMSCRIPTEN_RESULT result;
	if (p_enabled) {
		if (window_maximized) {
			// Soft fullsreen during real fulllscreen can cause issues.
			set_window_maximized(false);
			window_maximized = true;
		}
		EmscriptenFullscreenStrategy strategy;
		strategy.scaleMode = EMSCRIPTEN_FULLSCREEN_SCALE_STRETCH;
		strategy.canvasResolutionScaleMode = EMSCRIPTEN_FULLSCREEN_CANVAS_SCALE_STDDEF;
		strategy.filteringMode = EMSCRIPTEN_FULLSCREEN_FILTERING_DEFAULT;
		strategy.canvasResizedCallback = NULL;
		emscripten_request_fullscreen_strategy(NULL, false, &strategy);
	} else {
		result = emscripten_exit_fullscreen();
		if (result != EMSCRIPTEN_RESULT_SUCCESS) {
			ERR_PRINTS("Failed to exit fullscreen: Code " + itos(result));
		}
	}
}

bool OS_JavaScript::is_window_fullscreen() const {

	return video_mode.fullscreen;
}

void OS_JavaScript::get_fullscreen_mode_list(List<VideoMode> *p_list, int p_screen) const {

	Size2 screen = get_screen_size();
	p_list->push_back(OS::VideoMode(screen.width, screen.height, true));
}

// Keys

template <typename T>
static void dom2godot_mod(T *emscripten_event_ptr, Ref<InputEventWithModifiers> godot_event) {

	godot_event->set_shift(emscripten_event_ptr->shiftKey);
	godot_event->set_alt(emscripten_event_ptr->altKey);
	godot_event->set_control(emscripten_event_ptr->ctrlKey);
	godot_event->set_metakey(emscripten_event_ptr->metaKey);
}

static Ref<InputEventKey> setup_key_event(const EmscriptenKeyboardEvent *emscripten_event) {

	Ref<InputEventKey> ev;
	ev.instance();
	ev->set_echo(emscripten_event->repeat);
	dom2godot_mod(emscripten_event, ev);
	ev->set_scancode(dom2godot_scancode(emscripten_event->keyCode));

	String unicode = String::utf8(emscripten_event->key);
	// Check if empty or multi-character (e.g. `CapsLock`).
	if (unicode.length() != 1) {
		// Might be empty as well, but better than nonsense.
		unicode = String::utf8(emscripten_event->charValue);
	}
	if (unicode.length() == 1) {
		ev->set_unicode(unicode[0]);
	}

	return ev;
}

EM_BOOL OS_JavaScript::keydown_callback(int p_event_type, const EmscriptenKeyboardEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();
	Ref<InputEventKey> ev = setup_key_event(p_event);
	ev->set_pressed(true);
	if (ev->get_unicode() == 0 && keycode_has_unicode(ev->get_scancode())) {
		// Defer to keypress event for legacy unicode retrieval.
		os->deferred_key_event = ev;
		// Do not suppress keypress event.
		return false;
	}
	os->input->parse_input_event(ev);
	return true;
}

EM_BOOL OS_JavaScript::keypress_callback(int p_event_type, const EmscriptenKeyboardEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();
	os->deferred_key_event->set_unicode(p_event->charCode);
	os->input->parse_input_event(os->deferred_key_event);
	return true;
}

EM_BOOL OS_JavaScript::keyup_callback(int p_event_type, const EmscriptenKeyboardEvent *p_event, void *p_user_data) {

	Ref<InputEventKey> ev = setup_key_event(p_event);
	ev->set_pressed(false);
	get_singleton()->input->parse_input_event(ev);
	return ev->get_scancode() != KEY_UNKNOWN && ev->get_scancode() != 0;
}

// Mouse

Point2 OS_JavaScript::get_mouse_position() const {

	return input->get_mouse_position();
}

int OS_JavaScript::get_mouse_button_state() const {

	return input->get_mouse_button_mask();
}

EM_BOOL OS_JavaScript::mouse_button_callback(int p_event_type, const EmscriptenMouseEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();

	Ref<InputEventMouseButton> ev;
	ev.instance();
	ev->set_pressed(p_event_type == EMSCRIPTEN_EVENT_MOUSEDOWN);
	ev->set_position(Point2(p_event->canvasX, p_event->canvasY));
	ev->set_global_position(ev->get_position());
	dom2godot_mod(p_event, ev);
	switch (p_event->button) {
		case DOM_BUTTON_LEFT: ev->set_button_index(BUTTON_LEFT); break;
		case DOM_BUTTON_MIDDLE: ev->set_button_index(BUTTON_MIDDLE); break;
		case DOM_BUTTON_RIGHT: ev->set_button_index(BUTTON_RIGHT); break;
		case DOM_BUTTON_XBUTTON1: ev->set_button_index(BUTTON_XBUTTON1); break;
		case DOM_BUTTON_XBUTTON2: ev->set_button_index(BUTTON_XBUTTON2); break;
		default: return false;
	}

	int mask = os->input->get_mouse_button_mask();
	int button_flag = 1 << (ev->get_button_index() - 1);
	if (ev->is_pressed()) {
		// Since the event is consumed, focus manually. The containing iframe,
		// if exists, may not have focus yet, so focus even if already focused.
		focus_canvas();
		mask |= button_flag;
	} else if (mask & button_flag) {
		mask &= ~button_flag;
	} else {
		// Received release event, but press was outside the canvas, so ignore.
		return false;
	}
	ev->set_button_mask(mask);

	os->input->parse_input_event(ev);
	// Prevent multi-click text selection and wheel-click scrolling anchor.
	// Context menu is prevented through contextmenu event.
	return true;
}

EM_BOOL OS_JavaScript::mousemove_callback(int p_event_type, const EmscriptenMouseEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();

	int input_mask = os->input->get_mouse_button_mask();
	Point2 pos = Point2(p_event->canvasX, p_event->canvasY);
	// For motion outside the canvas, only read mouse movement if dragging
	// started inside the canvas; imitating desktop app behaviour.
	if (!cursor_inside_canvas && !input_mask)
		return false;

	Ref<InputEventMouseMotion> ev;
	ev.instance();
	dom2godot_mod(p_event, ev);
	ev->set_button_mask(input_mask);

	ev->set_position(pos);
	ev->set_global_position(ev->get_position());

	ev->set_relative(Vector2(p_event->movementX, p_event->movementY));
	os->input->set_mouse_position(ev->get_position());
	ev->set_speed(os->input->get_last_mouse_speed());

	os->input->parse_input_event(ev);
	// Don't suppress mouseover/-leave events.
	return false;
}

static const char *godot2dom_cursor(OS::CursorShape p_shape) {

	switch (p_shape) {
		case OS::CURSOR_ARROW:
		default:
			return "auto";
		case OS::CURSOR_IBEAM: return "text";
		case OS::CURSOR_POINTING_HAND: return "pointer";
		case OS::CURSOR_CROSS: return "crosshair";
		case OS::CURSOR_WAIT: return "progress";
		case OS::CURSOR_BUSY: return "wait";
		case OS::CURSOR_DRAG: return "grab";
		case OS::CURSOR_CAN_DROP: return "grabbing";
		case OS::CURSOR_FORBIDDEN: return "no-drop";
		case OS::CURSOR_VSIZE: return "ns-resize";
		case OS::CURSOR_HSIZE: return "ew-resize";
		case OS::CURSOR_BDIAGSIZE: return "nesw-resize";
		case OS::CURSOR_FDIAGSIZE: return "nwse-resize";
		case OS::CURSOR_MOVE: return "move";
		case OS::CURSOR_VSPLIT: return "row-resize";
		case OS::CURSOR_HSPLIT: return "col-resize";
		case OS::CURSOR_HELP: return "help";
	}
}

static void set_css_cursor(const char *p_cursor) {

	/* clang-format off */
	EM_ASM_({
		Module.canvas.style.cursor = UTF8ToString($0);
	}, p_cursor);
	/* clang-format on */
}

static const char *get_css_cursor() {

	char cursor[16];
	/* clang-format off */
	EM_ASM_INT({
		stringToUTF8(Module.canvas.style.cursor ? Module.canvas.style.cursor : 'auto', $0, 16);
	}, cursor);
	/* clang-format on */
	return cursor;
}

void OS_JavaScript::set_cursor_shape(CursorShape p_shape) {

	ERR_FAIL_INDEX(p_shape, CURSOR_MAX);

	cursor_shape = p_shape;
	if (get_mouse_mode() != MOUSE_MODE_HIDDEN)
		set_css_cursor(godot2dom_cursor(cursor_shape));
}

void OS_JavaScript::set_custom_mouse_cursor(const RES &p_cursor, CursorShape p_shape, const Vector2 &p_hotspot) {
}

void OS_JavaScript::set_mouse_mode(OS::MouseMode p_mode) {

	ERR_EXPLAIN("MOUSE_MODE_CONFINED is not supported for the HTML5 platform");
	ERR_FAIL_COND(p_mode == MOUSE_MODE_CONFINED);
	if (p_mode == get_mouse_mode())
		return;

	if (p_mode == MOUSE_MODE_VISIBLE) {

		set_css_cursor(godot2dom_cursor(cursor_shape));
		emscripten_exit_pointerlock();

	} else if (p_mode == MOUSE_MODE_HIDDEN) {

		set_css_cursor("none");
		emscripten_exit_pointerlock();

	} else if (p_mode == MOUSE_MODE_CAPTURED) {

		EMSCRIPTEN_RESULT result = emscripten_request_pointerlock("canvas", false);
		ERR_EXPLAIN("MOUSE_MODE_CAPTURED can only be entered from within an appropriate input callback");
		ERR_FAIL_COND(result == EMSCRIPTEN_RESULT_FAILED_NOT_DEFERRED);
		ERR_FAIL_COND(result != EMSCRIPTEN_RESULT_SUCCESS);
		set_css_cursor(godot2dom_cursor(cursor_shape));
	}
}

OS::MouseMode OS_JavaScript::get_mouse_mode() const {

	if (String::utf8(get_css_cursor()) == "none")
		return MOUSE_MODE_HIDDEN;

	EmscriptenPointerlockChangeEvent ev;
	emscripten_get_pointerlock_status(&ev);
	return (ev.isActive && String::utf8(ev.id) == "canvas") ? MOUSE_MODE_CAPTURED : MOUSE_MODE_VISIBLE;
}

// Wheel

EM_BOOL OS_JavaScript::wheel_callback(int p_event_type, const EmscriptenWheelEvent *p_event, void *p_user_data) {

	ERR_FAIL_COND_V(p_event_type != EMSCRIPTEN_EVENT_WHEEL, false);
	if (!is_canvas_focused()) {
		if (cursor_inside_canvas) {
			focus_canvas();
		} else {
			return false;
		}
	}

	InputDefault *input = get_singleton()->input;
	Ref<InputEventMouseButton> ev;
	ev.instance();
	ev->set_button_mask(input->get_mouse_button_mask());
	ev->set_position(input->get_mouse_position());
	ev->set_global_position(ev->get_position());

	ev->set_shift(input->is_key_pressed(KEY_SHIFT));
	ev->set_alt(input->is_key_pressed(KEY_ALT));
	ev->set_control(input->is_key_pressed(KEY_CONTROL));
	ev->set_metakey(input->is_key_pressed(KEY_META));

	if (p_event->deltaY < 0)
		ev->set_button_index(BUTTON_WHEEL_UP);
	else if (p_event->deltaY > 0)
		ev->set_button_index(BUTTON_WHEEL_DOWN);
	else if (p_event->deltaX > 0)
		ev->set_button_index(BUTTON_WHEEL_LEFT);
	else if (p_event->deltaX < 0)
		ev->set_button_index(BUTTON_WHEEL_RIGHT);
	else
		return false;

	// Different browsers give wildly different delta values, and we can't
	// interpret deltaMode, so use default value for wheel events' factor.

	ev->set_pressed(true);
	input->parse_input_event(ev);

	ev->set_pressed(false);
	input->parse_input_event(ev);

	return true;
}

// Touch

bool OS_JavaScript::has_touchscreen_ui_hint() const {

	/* clang-format off */
	return EM_ASM_INT_V(
		return 'ontouchstart' in window;
	);
	/* clang-format on */
}

EM_BOOL OS_JavaScript::touch_press_callback(int p_event_type, const EmscriptenTouchEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();
	Ref<InputEventScreenTouch> ev;
	ev.instance();
	int lowest_id_index = -1;
	for (int i = 0; i < p_event->numTouches; ++i) {

		const EmscriptenTouchPoint &touch = p_event->touches[i];
		if (lowest_id_index == -1 || touch.identifier < p_event->touches[lowest_id_index].identifier)
			lowest_id_index = i;
		if (!touch.isChanged)
			continue;
		ev->set_index(touch.identifier);
		ev->set_position(Point2(touch.canvasX, touch.canvasY));
		os->touches[i] = ev->get_position();
		ev->set_pressed(p_event_type == EMSCRIPTEN_EVENT_TOUCHSTART);

		os->input->parse_input_event(ev);
	}
	return true;
}

EM_BOOL OS_JavaScript::touchmove_callback(int p_event_type, const EmscriptenTouchEvent *p_event, void *p_user_data) {

	OS_JavaScript *os = get_singleton();
	Ref<InputEventScreenDrag> ev;
	ev.instance();
	int lowest_id_index = -1;
	for (int i = 0; i < p_event->numTouches; ++i) {

		const EmscriptenTouchPoint &touch = p_event->touches[i];
		if (lowest_id_index == -1 || touch.identifier < p_event->touches[lowest_id_index].identifier)
			lowest_id_index = i;
		if (!touch.isChanged)
			continue;
		ev->set_index(touch.identifier);
		ev->set_position(Point2(touch.canvasX, touch.canvasY));
		Point2 &prev = os->touches[i];
		ev->set_relative(ev->get_position() - prev);
		prev = ev->get_position();

		os->input->parse_input_event(ev);
	}
	return true;
}

// Gamepad

EM_BOOL OS_JavaScript::gamepad_change_callback(int p_event_type, const EmscriptenGamepadEvent *p_event, void *p_user_data) {

	InputDefault *input = get_singleton()->input;
	if (p_event_type == EMSCRIPTEN_EVENT_GAMEPADCONNECTED) {

		String guid = "";
		if (String::utf8(p_event->mapping) == "standard")
			guid = "Default HTML5 Gamepad";
		input->joy_connection_changed(p_event->index, true, String::utf8(p_event->id), guid);
	} else {
		input->joy_connection_changed(p_event->index, false, "");
	}
	return true;
}

void OS_JavaScript::process_joypads() {

	int joypad_count = emscripten_get_num_gamepads();
	for (int joypad = 0; joypad < joypad_count; joypad++) {
		EmscriptenGamepadEvent state;
		emscripten_get_gamepad_status(joypad, &state);
		if (state.connected) {

			int button_count = MIN(state.numButtons, 18);
			int axis_count = MIN(state.numAxes, 8);
			for (int button = 0; button < button_count; button++) {

				float value = state.analogButton[button];
				if (String::utf8(state.mapping) == "standard" && (button == JOY_ANALOG_L2 || button == JOY_ANALOG_R2)) {
					InputDefault::JoyAxis joy_axis;
					joy_axis.min = 0;
					joy_axis.value = value;
					input->joy_axis(joypad, button, joy_axis);
				} else {
					input->joy_button(joypad, button, value);
				}
			}
			for (int axis = 0; axis < axis_count; axis++) {

				InputDefault::JoyAxis joy_axis;
				joy_axis.min = -1;
				joy_axis.value = state.axis[axis];
				input->joy_axis(joypad, axis, joy_axis);
			}
		}
	}
}

bool OS_JavaScript::is_joy_known(int p_device) {

	return input->is_joy_mapped(p_device);
}

String OS_JavaScript::get_joy_guid(int p_device) const {

	return input->get_joy_guid_remapped(p_device);
}

// Video

int OS_JavaScript::get_video_driver_count() const {

	return VIDEO_DRIVER_MAX;
}

const char *OS_JavaScript::get_video_driver_name(int p_driver) const {

	switch (p_driver) {
		case VIDEO_DRIVER_GLES3:
			return "GLES3";
		case VIDEO_DRIVER_GLES2:
			return "GLES2";
	}
	ERR_EXPLAIN("Invalid video driver index " + itos(p_driver));
	ERR_FAIL_V(NULL);
}

// Audio

int OS_JavaScript::get_audio_driver_count() const {

	return 1;
}

const char *OS_JavaScript::get_audio_driver_name(int p_driver) const {

	return "JavaScript";
}

// Lifecycle
int OS_JavaScript::get_current_video_driver() const {
	return video_driver_index;
}

void OS_JavaScript::initialize_core() {

	OS_Unix::initialize_core();
	FileAccess::make_default<FileAccessBufferedFA<FileAccessUnix> >(FileAccess::ACCESS_RESOURCES);
}

Error OS_JavaScript::initialize(const VideoMode &p_desired, int p_video_driver, int p_audio_driver) {

	EmscriptenWebGLContextAttributes attributes;
	emscripten_webgl_init_context_attributes(&attributes);
	attributes.alpha = false;
	attributes.antialias = false;
	ERR_FAIL_INDEX_V(p_video_driver, VIDEO_DRIVER_MAX, ERR_INVALID_PARAMETER);
	switch (p_video_driver) {
		case VIDEO_DRIVER_GLES3:
			attributes.majorVersion = 2;
			RasterizerGLES3::register_config();
			RasterizerGLES3::make_current();
			break;
		case VIDEO_DRIVER_GLES2:
			attributes.majorVersion = 1;
			RasterizerGLES2::register_config();
			RasterizerGLES2::make_current();
			break;
	}

	video_driver_index = p_video_driver;
	EMSCRIPTEN_WEBGL_CONTEXT_HANDLE ctx = emscripten_webgl_create_context(NULL, &attributes);
	ERR_EXPLAIN("WebGL " + itos(attributes.majorVersion) + ".0 not available");
	ERR_FAIL_COND_V(emscripten_webgl_make_context_current(ctx) != EMSCRIPTEN_RESULT_SUCCESS, ERR_UNAVAILABLE);

	video_mode = p_desired;
	// Can't fulfil fullscreen request during start-up due to browser security.
	video_mode.fullscreen = false;
	/* clang-format off */
	if (EM_ASM_INT_V({ return Module.resizeCanvasOnStart })) {
		/* clang-format on */
		set_window_size(Size2(video_mode.width, video_mode.height));
	} else {
		set_window_size(get_window_size());
	}

	char locale_ptr[16];
	/* clang-format off */
	EM_ASM_ARGS({
		stringToUTF8(Module.locale, $0, 16);
	}, locale_ptr);
	/* clang-format on */
	setenv("LANG", locale_ptr, true);

	AudioDriverManager::initialize(p_audio_driver);
	VisualServer *visual_server = memnew(VisualServerRaster());
	input = memnew(InputDefault);

	EMSCRIPTEN_RESULT result;
#define EM_CHECK(ev)                         \
	if (result != EMSCRIPTEN_RESULT_SUCCESS) \
	ERR_PRINTS("Error while setting " #ev " callback: Code " + itos(result))
#define SET_EM_CALLBACK(target, ev, cb)                               \
	result = emscripten_set_##ev##_callback(target, NULL, true, &cb); \
	EM_CHECK(ev)
#define SET_EM_CALLBACK_NOTARGET(ev, cb)                      \
	result = emscripten_set_##ev##_callback(NULL, true, &cb); \
	EM_CHECK(ev)
	// These callbacks from Emscripten's html5.h suffice to access most
	// JavaScript APIs. For APIs that are not (sufficiently) exposed, EM_ASM
	// is used below.
	SET_EM_CALLBACK("#window", mousemove, mousemove_callback)
	SET_EM_CALLBACK("#canvas", mousedown, mouse_button_callback)
	SET_EM_CALLBACK("#window", mouseup, mouse_button_callback)
	SET_EM_CALLBACK("#window", wheel, wheel_callback)
	SET_EM_CALLBACK("#window", touchstart, touch_press_callback)
	SET_EM_CALLBACK("#window", touchmove, touchmove_callback)
	SET_EM_CALLBACK("#window", touchend, touch_press_callback)
	SET_EM_CALLBACK("#window", touchcancel, touch_press_callback)
	SET_EM_CALLBACK("#canvas", keydown, keydown_callback)
	SET_EM_CALLBACK("#canvas", keypress, keypress_callback)
	SET_EM_CALLBACK("#canvas", keyup, keyup_callback)
	SET_EM_CALLBACK(NULL, resize, browser_resize_callback)
	SET_EM_CALLBACK(NULL, fullscreenchange, fullscreen_change_callback)
	SET_EM_CALLBACK_NOTARGET(gamepadconnected, gamepad_change_callback)
	SET_EM_CALLBACK_NOTARGET(gamepaddisconnected, gamepad_change_callback)
#undef SET_EM_CALLBACK_NODATA
#undef SET_EM_CALLBACK
#undef EM_CHECK

	/* clang-format off */
	EM_ASM_ARGS({
		const send_notification = cwrap('send_notification', null, ['number']);
		const notifications = arguments;
		(['mouseover', 'mouseleave', 'focus', 'blur']).forEach(function(event, index) {
			Module.canvas.addEventListener(event, send_notification.bind(null, notifications[index]));
		});
	},
		MainLoop::NOTIFICATION_WM_MOUSE_ENTER,
		MainLoop::NOTIFICATION_WM_MOUSE_EXIT,
		MainLoop::NOTIFICATION_WM_FOCUS_IN,
		MainLoop::NOTIFICATION_WM_FOCUS_OUT
	);
	/* clang-format on */

	visual_server->init();

	return OK;
}

void OS_JavaScript::set_main_loop(MainLoop *p_main_loop) {

	main_loop = p_main_loop;
	input->set_main_loop(p_main_loop);
}

MainLoop *OS_JavaScript::get_main_loop() const {

	return main_loop;
}

void OS_JavaScript::run_async() {

	main_loop->init();
	emscripten_set_main_loop(main_loop_callback, -1, false);
}

void OS_JavaScript::main_loop_callback() {

	get_singleton()->main_loop_iterate();
}

bool OS_JavaScript::main_loop_iterate() {

	if (is_userfs_persistent() && sync_wait_time >= 0) {
		int64_t current_time = get_ticks_msec();
		int64_t elapsed_time = current_time - last_sync_check_time;
		last_sync_check_time = current_time;

		sync_wait_time -= elapsed_time;

		if (sync_wait_time < 0) {
			/* clang-format off */
			EM_ASM(
				FS.syncfs(function(err) {
					if (err) { console.warn('Failed to save IDB file system: ' + err.message); }
				});
			);
			/* clang-format on */
		}
	}
	process_joypads();
	if (canvas_size_adjustment_requested) {

		if (video_mode.fullscreen || window_maximized) {
			video_mode.width = get_window_size().width;
			video_mode.height = get_window_size().height;
		}
		if (!video_mode.fullscreen) {
			set_window_maximized(window_maximized);
		}
		canvas_size_adjustment_requested = false;
	}
	return Main::iteration();
}

void OS_JavaScript::delete_main_loop() {

	memdelete(main_loop);
}

void OS_JavaScript::finalize() {

	memdelete(input);
}

// Miscellaneous

extern "C" EMSCRIPTEN_KEEPALIVE void send_notification(int p_notification) {

	if (p_notification == MainLoop::NOTIFICATION_WM_MOUSE_ENTER || p_notification == MainLoop::NOTIFICATION_WM_MOUSE_EXIT) {
		cursor_inside_canvas = p_notification == MainLoop::NOTIFICATION_WM_MOUSE_ENTER;
	}
	OS_JavaScript::get_singleton()->get_main_loop()->notification(p_notification);
}

bool OS_JavaScript::_check_internal_feature_support(const String &p_feature) {

	if (p_feature == "HTML5" || p_feature == "web")
		return true;

#ifdef JAVASCRIPT_EVAL_ENABLED
	if (p_feature == "JavaScript")
		return true;
#endif

	EMSCRIPTEN_WEBGL_CONTEXT_HANDLE ctx = emscripten_webgl_get_current_context();
	// All extensions are already automatically enabled, this function allows
	// checking WebGL extension support without inline JavaScript
	if (p_feature == "s3tc")
		return emscripten_webgl_enable_extension(ctx, "WEBGL_compressed_texture_s3tc_srgb");
	if (p_feature == "etc")
		return emscripten_webgl_enable_extension(ctx, "WEBGL_compressed_texture_etc1");
	if (p_feature == "etc2")
		return emscripten_webgl_enable_extension(ctx, "WEBGL_compressed_texture_etc");

	return false;
}

void OS_JavaScript::alert(const String &p_alert, const String &p_title) {

	/* clang-format off */
	EM_ASM_({
		window.alert(UTF8ToString($0));
	}, p_alert.utf8().get_data());
	/* clang-format on */
}

void OS_JavaScript::set_window_title(const String &p_title) {

	/* clang-format off */
	EM_ASM_({
		document.title = UTF8ToString($0);
	}, p_title.utf8().get_data());
	/* clang-format on */
}

String OS_JavaScript::get_executable_path() const {

	return OS::get_executable_path();
}

Error OS_JavaScript::shell_open(String p_uri) {

	// Open URI in a new tab, browser will deal with it by protocol.
	/* clang-format off */
	EM_ASM_({
		window.open(UTF8ToString($0), '_blank');
	}, p_uri.utf8().get_data());
	/* clang-format on */
	return OK;
}

String OS_JavaScript::get_name() {

	return "HTML5";
}

bool OS_JavaScript::can_draw() const {

	return true; // Always?
}

String OS_JavaScript::get_user_data_dir() const {

	return "/userfs";
};

String OS_JavaScript::get_resource_dir() const {

	return "/";
}

OS::PowerState OS_JavaScript::get_power_state() {

	WARN_PRINT("Power management is not supported for the HTML5 platform, defaulting to POWERSTATE_UNKNOWN");
	return OS::POWERSTATE_UNKNOWN;
}

int OS_JavaScript::get_power_seconds_left() {

	WARN_PRINT("Power management is not supported for the HTML5 platform, defaulting to -1");
	return -1;
}

int OS_JavaScript::get_power_percent_left() {

	WARN_PRINT("Power management is not supported for the HTML5 platform, defaulting to -1");
	return -1;
}

void OS_JavaScript::file_access_close_callback(const String &p_file, int p_flags) {

	OS_JavaScript *os = get_singleton();
	if (os->is_userfs_persistent() && p_file.begins_with("/userfs") && p_flags & FileAccess::WRITE) {
		os->last_sync_check_time = OS::get_singleton()->get_ticks_msec();
		// Wait five seconds in case more files are about to be closed.
		os->sync_wait_time = 5000;
	}
}

void OS_JavaScript::set_idb_available(bool p_idb_available) {

	idb_available = p_idb_available;
}

bool OS_JavaScript::is_userfs_persistent() const {

	return idb_available;
}

OS_JavaScript *OS_JavaScript::get_singleton() {

	return static_cast<OS_JavaScript *>(OS::get_singleton());
}

OS_JavaScript::OS_JavaScript(int p_argc, char *p_argv[]) {

	List<String> arguments;
	for (int i = 1; i < p_argc; i++) {
		arguments.push_back(String::utf8(p_argv[i]));
	}
	set_cmdline(p_argv[0], arguments);

	window_maximized = false;
	soft_fullscreen_enabled = false;
	canvas_size_adjustment_requested = false;

	main_loop = NULL;

	idb_available = false;
	sync_wait_time = -1;

	AudioDriverManager::add_driver(&audio_driver_javascript);

	Vector<Logger *> loggers;
	loggers.push_back(memnew(StdLogger));
	_set_logger(memnew(CompositeLogger(loggers)));

	FileAccessUnix::close_notification_func = file_access_close_callback;
}
