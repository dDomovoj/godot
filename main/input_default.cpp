/*************************************************************************/
/*  input_default.cpp                                                    */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2017 Juan Linietsky, Ariel Manzur.                 */
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
#include "input_default.h"
#include "servers/visual_server.h"
#include "os/os.h"
#include "input_map.h"
#include "scene/resources/texture.h"

void InputDefault::SpeedTrack::update(const Vector2& p_delta_p) {

	uint64_t tick = OS::get_singleton()->get_ticks_usec();
	uint32_t tdiff = tick-last_tick;
	float delta_t = tdiff / 1000000.0;
	last_tick=tick;


	accum+=p_delta_p;
	accum_t+=delta_t;

	if (accum_t>max_ref_frame*10)
		accum_t=max_ref_frame*10;

	while( accum_t>=min_ref_frame ) {

		float slice_t = min_ref_frame / accum_t;
		Vector2 slice = accum*slice_t;
		accum=accum-slice;
		accum_t-=min_ref_frame;

		speed=(slice/min_ref_frame).linear_interpolate(speed,min_ref_frame/max_ref_frame);
	}

}

void InputDefault::SpeedTrack::reset() {
	last_tick = OS::get_singleton()->get_ticks_usec();
	speed=Vector2();
	accum_t=0;
}

InputDefault::SpeedTrack::SpeedTrack() {

	 min_ref_frame=0.1;
	 max_ref_frame=0.3;
	 reset();
}

bool InputDefault::is_key_pressed(int p_scancode) const {

	_THREAD_SAFE_METHOD_
	return keys_pressed.has(p_scancode);
}

bool InputDefault::is_mouse_button_pressed(int p_button) const {

	_THREAD_SAFE_METHOD_
	return (mouse_button_mask&(1<<p_button))!=0;
}


static int _combine_device(int p_value,int p_device) {

	return p_value|(p_device<<20);
}

bool InputDefault::is_joy_button_pressed(int p_device, int p_button) const{

	_THREAD_SAFE_METHOD_
	return joy_buttons_pressed.has(_combine_device(p_button,p_device));
}

bool InputDefault::is_action_pressed(const StringName& p_action) const{

	return action_state.has(p_action) && action_state[p_action].pressed;
#if 0
	if (custom_action_press.has(p_action))
		return true; //simpler

	const List<InputEvent> *alist = InputMap::get_singleton()->get_action_list(p_action);
	if (!alist)
		return false;


	for (const List<InputEvent>::Element *E=alist->front();E;E=E->next()) {


		int device=E->get().device;

		switch(E->get().type) {

			case InputEvent::KEY: {

				const InputEventKey &iek=E->get().key;
				if ((keys_pressed.has(iek.scancode)))
					return true;
			} break;
			case InputEvent::MOUSE_BUTTON: {

				const InputEventMouseButton &iemb=E->get().mouse_button;
				 if(mouse_button_mask&(1<<iemb.button_index))
					 return true;
			} break;
			case InputEvent::JOYPAD_BUTTON: {

				const InputEventJoypadButton &iejb=E->get().joy_button;
				int c = _combine_device(iejb.button_index,device);
				if (joy_buttons_pressed.has(c))
					return true;
			} break;
			case InputEvent::JOYPAD_MOTION: {

				const InputEventJoypadMotion &iejm=E->get().joy_motion;
				int c = _combine_device(iejm.axis,device);
				if (_joy_axis.has(c)) {
					if (iejm.axis_value < 0) {
						if (_joy_axis[c] < -0.5f) return true;
					}
					else
						if (_joy_axis[c] > 0.5f) return true;
				}
			} break;
		}
	}

	return false;
#endif
}

bool InputDefault::is_action_just_pressed(const StringName& p_action) const {

	const Map<StringName,Action>::Element *E=action_state.find(p_action);
	if (!E)
		return false;

	if (Engine::get_singleton()->is_in_fixed_frame()) {
		return E->get().pressed && E->get().fixed_frame==Engine::get_singleton()->get_fixed_frames();
	} else {
		return E->get().pressed && E->get().idle_frame==Engine::get_singleton()->get_idle_frames();
	}
}

bool InputDefault::is_action_just_released(const StringName& p_action) const{

	const Map<StringName,Action>::Element *E=action_state.find(p_action);
	if (!E)
		return false;

	if (Engine::get_singleton()->is_in_fixed_frame()) {
		return !E->get().pressed && E->get().fixed_frame==Engine::get_singleton()->get_fixed_frames();
	} else {
		return !E->get().pressed && E->get().idle_frame==Engine::get_singleton()->get_idle_frames();
	}
}


float InputDefault::get_joy_axis(int p_device,int p_axis) const{

	_THREAD_SAFE_METHOD_
	int c = _combine_device(p_axis,p_device);
	if (_joy_axis.has(c)) {
		return _joy_axis[c];
	} else {
		return 0;
	}
}

String InputDefault::get_joy_name(int p_idx) {

	_THREAD_SAFE_METHOD_
	return joy_names[p_idx].name;
};

Vector2 InputDefault::get_joy_vibration_strength(int p_device) {
	if (joy_vibration.has(p_device)) {
		return Vector2(joy_vibration[p_device].weak_magnitude, joy_vibration[p_device].strong_magnitude);
	} else {
		return Vector2(0, 0);
	}
}

uint64_t InputDefault::get_joy_vibration_timestamp(int p_device) {
	if (joy_vibration.has(p_device)) {
		return joy_vibration[p_device].timestamp;
	} else {
		return 0;
	}
}

float InputDefault::get_joy_vibration_duration(int p_device) {
	if (joy_vibration.has(p_device)) {
		return joy_vibration[p_device].duration;
	} else {
		return 0.f;
	}
}

static String _hex_str(uint8_t p_byte) {

	static const char* dict = "0123456789abcdef";
	char ret[3];
	ret[2] = 0;

	ret[0] = dict[p_byte>>4];
	ret[1] = dict[p_byte & 0xf];

	return ret;
};

void InputDefault::joy_connection_changed(int p_idx, bool p_connected, String p_name, String p_guid) {

	_THREAD_SAFE_METHOD_
	Joypad js;
	js.name = p_connected ? p_name : "";
	js.uid = p_connected ? p_guid : "";
	js.mapping = -1;
	js.hat_current = 0;

	if (p_connected) {

		String uidname = p_guid;
		if (p_guid == "") {
			int uidlen = MIN(p_name.length(), 16);
			for (int i=0; i<uidlen; i++) {
				uidname = uidname + _hex_str(p_name[i]);
			};
		};
		js.uid = uidname;
		js.connected = true;
		int mapping = fallback_mapping;
		for (int i=0; i < map_db.size(); i++) {
			if (js.uid == map_db[i].uid) {
				mapping = i;
				js.name = map_db[i].name;
			};
		};
		js.mapping = mapping;
	}
	else {
		js.connected = false;
		for (int i = 0; i < JOY_BUTTON_MAX; i++) {

			if (i < JOY_AXIS_MAX)
				set_joy_axis(p_idx, i, 0.0f);

			int c = _combine_device(i, p_idx);
			joy_buttons_pressed.erase(c);
		};
	};
	joy_names[p_idx] = js;

	emit_signal("joy_connection_changed", p_idx, p_connected);
};

Vector3 InputDefault::get_gravity() const{

	_THREAD_SAFE_METHOD_
	return gravity;
}

Vector3 InputDefault::get_accelerometer() const{

	_THREAD_SAFE_METHOD_
	return accelerometer;
}

Vector3 InputDefault::get_magnetometer() const{

	_THREAD_SAFE_METHOD_
	return magnetometer;
}

Vector3 InputDefault::get_gyroscope() const {

	_THREAD_SAFE_METHOD_
	return gyroscope;
}

void InputDefault::parse_input_event(const InputEvent& p_event) {

	_THREAD_SAFE_METHOD_
	switch(p_event.type) {

		case InputEvent::KEY: {

			if (p_event.key.echo)
				break;
			if (p_event.key.scancode==0)
				break;

			//print_line(p_event);

			if (p_event.key.pressed)
				keys_pressed.insert(p_event.key.scancode);
			else
				keys_pressed.erase(p_event.key.scancode);
		} break;
		case InputEvent::MOUSE_BUTTON: {

			if (p_event.mouse_button.doubleclick)
				break;

			if (p_event.mouse_button.pressed)
				mouse_button_mask|=(1<<p_event.mouse_button.button_index);
			else
				mouse_button_mask&=~(1<<p_event.mouse_button.button_index);

			if (main_loop && emulate_touch && p_event.mouse_button.button_index==1) {
				InputEventScreenTouch touch_event;
				touch_event.index=0;
				touch_event.pressed=p_event.mouse_button.pressed;
				touch_event.x=p_event.mouse_button.x;
				touch_event.y=p_event.mouse_button.y;
				InputEvent ev;
				ev.type=InputEvent::SCREEN_TOUCH;
				ev.screen_touch=touch_event;
				main_loop->input_event(ev);
			}
		} break;
		case InputEvent::MOUSE_MOTION: {

			if (main_loop && emulate_touch && p_event.mouse_motion.button_mask&1) {
				InputEventScreenDrag drag_event;
				drag_event.index=0;
				drag_event.x=p_event.mouse_motion.x;
				drag_event.y=p_event.mouse_motion.y;
				drag_event.relative_x=p_event.mouse_motion.relative_x;
				drag_event.relative_y=p_event.mouse_motion.relative_y;
				drag_event.speed_x=p_event.mouse_motion.speed_x;
				drag_event.speed_y=p_event.mouse_motion.speed_y;

				InputEvent ev;
				ev.type=InputEvent::SCREEN_DRAG;
				ev.screen_drag=drag_event;

				main_loop->input_event(ev);
			}

		} break;
		case InputEvent::JOYPAD_BUTTON: {

			int c = _combine_device(p_event.joy_button.button_index,p_event.device);

			if (p_event.joy_button.pressed)
				joy_buttons_pressed.insert(c);
			else
				joy_buttons_pressed.erase(c);
		} break;
		case InputEvent::JOYPAD_MOTION: {
			set_joy_axis(p_event.device, p_event.joy_motion.axis, p_event.joy_motion.axis_value);
		} break;

	}


	if (!p_event.is_echo()) {
		for (const Map<StringName,InputMap::Action>::Element *E=InputMap::get_singleton()->get_action_map().front();E;E=E->next()) {

			if (InputMap::get_singleton()->event_is_action(p_event,E->key()) && is_action_pressed(E->key()) != p_event.is_pressed()) {
				Action action;
				action.fixed_frame=Engine::get_singleton()->get_fixed_frames();
				action.idle_frame=Engine::get_singleton()->get_idle_frames();
				action.pressed=p_event.is_pressed();
				action_state[E->key()]=action;
			}
		}
	}

	if (main_loop)
		main_loop->input_event(p_event);

}

void InputDefault::set_joy_axis(int p_device,int p_axis,float p_value) {

	_THREAD_SAFE_METHOD_
	int c = _combine_device(p_axis,p_device);
	_joy_axis[c]=p_value;
}

void InputDefault::start_joy_vibration(int p_device, float p_weak_magnitude, float p_strong_magnitude, float p_duration) {
	_THREAD_SAFE_METHOD_
	if (p_weak_magnitude < 0.f || p_weak_magnitude > 1.f || p_strong_magnitude < 0.f || p_strong_magnitude > 1.f) {
		return;
	}
	VibrationInfo vibration;
	vibration.weak_magnitude = p_weak_magnitude;
	vibration.strong_magnitude = p_strong_magnitude;
	vibration.duration = p_duration;
	vibration.timestamp = OS::get_singleton()->get_ticks_usec();
	joy_vibration[p_device] = vibration;
}

void InputDefault::stop_joy_vibration(int p_device) {
	_THREAD_SAFE_METHOD_
	VibrationInfo vibration;
	vibration.weak_magnitude = 0;
	vibration.strong_magnitude = 0;
	vibration.duration = 0;
	vibration.timestamp = OS::get_singleton()->get_ticks_usec();
	joy_vibration[p_device] = vibration;
}

void InputDefault::set_gravity(const Vector3& p_gravity) {

	_THREAD_SAFE_METHOD_

	gravity=p_gravity;

}

void InputDefault::set_accelerometer(const Vector3& p_accel) {

	_THREAD_SAFE_METHOD_

	accelerometer=p_accel;

}

void InputDefault::set_magnetometer(const Vector3& p_magnetometer) {

	_THREAD_SAFE_METHOD_

	magnetometer=p_magnetometer;

}

void InputDefault::set_gyroscope(const Vector3& p_gyroscope) {

	_THREAD_SAFE_METHOD_

	gyroscope=p_gyroscope;

}

void InputDefault::set_main_loop(MainLoop *p_main_loop) {
	main_loop=p_main_loop;

}

void InputDefault::set_mouse_pos(const Point2& p_posf) {

	mouse_speed_track.update(p_posf-mouse_pos);
	mouse_pos=p_posf;
	if (custom_cursor.is_valid()) {
		VisualServer::get_singleton()->cursor_set_pos(get_mouse_pos());
	}
}

Point2 InputDefault::get_mouse_pos() const {

	return mouse_pos;
}
Point2 InputDefault::get_last_mouse_speed() const {

	return mouse_speed_track.speed;
}

int InputDefault::get_mouse_button_mask() const {

	return mouse_button_mask;// do not trust OS implementaiton, should remove it - OS::get_singleton()->get_mouse_button_state();
}

void InputDefault::warp_mouse_pos(const Vector2& p_to) {

	OS::get_singleton()->warp_mouse_pos(p_to);
}


void InputDefault::iteration(float p_step) {


}

void InputDefault::action_press(const StringName& p_action) {

	Action action;

	action.fixed_frame=Engine::get_singleton()->get_fixed_frames();
	action.idle_frame=Engine::get_singleton()->get_idle_frames();
	action.pressed=true;

	action_state[p_action]=action;

}

void InputDefault::action_release(const StringName& p_action){

	Action action;

	action.fixed_frame=Engine::get_singleton()->get_fixed_frames();
	action.idle_frame=Engine::get_singleton()->get_idle_frames();
	action.pressed=false;

	action_state[p_action]=action;
}

void InputDefault::set_emulate_touch(bool p_emulate) {

	emulate_touch=p_emulate;
}

bool InputDefault::is_emulating_touchscreen() const {

	return emulate_touch;
}

void InputDefault::set_custom_mouse_cursor(const RES& p_cursor,const Vector2& p_hotspot) {
	if (custom_cursor==p_cursor)
		return;

	custom_cursor=p_cursor;

	if (p_cursor.is_null()) {
		set_mouse_mode(MOUSE_MODE_VISIBLE);
		VisualServer::get_singleton()->cursor_set_visible(false);
	} else {
		Ref<AtlasTexture> atex = custom_cursor;
		Rect2 region = atex.is_valid() ? atex->get_region() : Rect2();
		set_mouse_mode(MOUSE_MODE_HIDDEN);
		VisualServer::get_singleton()->cursor_set_visible(true);
		VisualServer::get_singleton()->cursor_set_texture(custom_cursor->get_rid(),p_hotspot, 0, region);
		VisualServer::get_singleton()->cursor_set_pos(get_mouse_pos());
	}
}

void InputDefault::set_mouse_in_window(bool p_in_window) {

	if (custom_cursor.is_valid()) {

		if (p_in_window) {
			set_mouse_mode(MOUSE_MODE_HIDDEN);
			VisualServer::get_singleton()->cursor_set_visible(true);
		} else {
			set_mouse_mode(MOUSE_MODE_VISIBLE);
			VisualServer::get_singleton()->cursor_set_visible(false);
		}

	}
}

// from github.com/gabomdq/SDL_GameControllerDB
static const char *s_ControllerMappings [] = {
	#ifdef WINDOWS_ENABLED
	"00f00300000000000000504944564944,RetroUSB.com RetroPad,a:b1,b:b5,x:b0,y:b4,back:b2,start:b3,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"00f0f100000000000000504944564944,RetroUSB.com Super RetroPort,a:b1,b:b5,x:b0,y:b4,back:b2,start:b3,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"02200090000000000000504944564944,8Bitdo NES30 PRO USB,a:b0,b:b1,x:b3,y:b4,leftshoulder:b6,rightshoulder:b7,lefttrigger:b8,righttrigger:b9,back:b10,start:b11,leftstick:b13,rightstick:b14,leftx:a0,lefty:a1,rightx:a3,righty:a4,dpup:h0.1,dpright:h0.2,dpdown:h0.4,dpleft:h0.8,",
	"0d0f4900000000000000504944564944,Hatsune Miku Sho Controller,a:b1,b:b2,x:b0,y:b3,back:b8,guide:b12,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"0d0f6e00000000000000504944564944,HORIPAD 4,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"10080100000000000000504944564944,PS1 USB,a:b2,b:b1,x:b3,y:b0,back:b8,start:b9,leftshoulder:b6,rightshoulder:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a2,lefttrigger:b4,righttrigger:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,",
	"10080300000000000000504944564944,PS2 USB,a:b2,b:b1,y:b0,x:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a4,righty:a2,lefttrigger:b4,righttrigger:b5,",
	"10280900000000000000504944564944,8Bitdo SFC30 GamePad,a:b1,b:b0,y:b3,x:b4,start:b11,back:b10,leftshoulder:b6,leftx:a0,lefty:a1,rightshoulder:b7,",
	"20380900000000000000504944564944,8Bitdo NES30 PRO Wireless,a:b0,b:b1,x:b3,y:b4,leftshoulder:b6,rightshoulder:b7,lefttrigger:b8,righttrigger:b9,back:b10,start:b11,leftstick:b13,rightstick:b14,leftx:a0,lefty:a1,rightx:a3,righty:a4,dpup:h0.1,dpright:h0.2,dpdown:h0.4,dpleft:h0.8,",
	"25090500000000000000504944564944,PS3 DualShock,a:b2,b:b1,back:b9,dpdown:h0.8,dpleft:h0.4,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b6,leftstick:b10,lefttrigger:b4,leftx:a0,lefty:a1,rightshoulder:b7,rightstick:b11,righttrigger:b5,rightx:a2,righty:a3,start:b8,x:b0,y:b3,",
	"2509e803000000000000504944564944,Mayflash Wii Classic Controller,a:b1,b:b0,x:b3,y:b2,back:b8,guide:b10,start:b9,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:b11,dpdown:b13,dpleft:b12,dpright:b14,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"28040140000000000000504944564944,GamePad Pro USB,a:b1,b:b2,x:b0,y:b3,back:b8,start:b9,leftshoulder:b4,rightshoulder:b5,leftx:a0,lefty:a1,lefttrigger:b6,righttrigger:b7,",
	"300f1001000000000000504944564944,Saitek P480 Rumble Pad,a:b2,b:b3,x:b0,y:b1,back:b8,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a2,lefttrigger:b5,righttrigger:b7,",
	"341a0108000000000000504944564944,EXEQ RF USB Gamepad 8206,a:b0,b:b1,x:b2,y:b3,leftshoulder:b4,rightshoulder:b5,leftstick:b8,rightstick:b7,back:b8,start:b9,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"341a3608000000000000504944564944,Afterglow PS3 Controller,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b12,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"36280100000000000000504944564944,OUYA Controller,a:b0,b:b3,y:b2,x:b1,start:b14,guide:b15,leftstick:b6,rightstick:b7,leftshoulder:b4,rightshoulder:b5,dpup:b8,dpleft:b10,dpdown:b9,dpright:b11,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:b12,righttrigger:b13,",
	"49190204000000000000504944564944,Ipega PG-9023,a:b0,b:b1,x:b3,y:b4,back:b10,start:b11,leftstick:b13,rightstick:b14,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:b8,righttrigger:b9",
	"4b12014d000000000000504944564944,NYKO AIRFLO,a:b0,b:b1,x:b2,y:b3,back:b8,guide:b10,start:b9,leftstick:a0,rightstick:a2,leftshoulder:a3,rightshoulder:b5,dpup:h0.1,dpdown:h0.0,dpleft:h0.8,dpright:h0.2,leftx:h0.6,lefty:h0.12,rightx:h0.9,righty:h0.4,lefttrigger:b6,righttrigger:b7,",
	"4c056802000000000000504944564944,PS3 Controller,a:b14,b:b13,back:b0,dpdown:b6,dpleft:b7,dpright:b5,dpup:b4,guide:b16,leftshoulder:b10,leftstick:b1,lefttrigger:b8,leftx:a0,lefty:a1,rightshoulder:b11,rightstick:b2,righttrigger:b9,rightx:a2,righty:a3,start:b3,x:b15,y:b12,",
	"4c05c405000000000000504944564944,PS4 Controller,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b12,leftshoulder:b4,leftstick:b10,lefttrigger:a3,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:a4,rightx:a2,righty:a5,start:b9,x:b0,y:b3,",
	"4f0400b3000000000000504944564944,Thrustmaster Firestorm Dual Power,a:b0,b:b2,y:b3,x:b1,start:b10,guide:b8,back:b9,leftstick:b11,rightstick:b12,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b5,righttrigger:b7,",
	"4f0415b3000000000000504944564944,Thrustmaster Dual Analog 3.2,x:b1,a:b0,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b5,rightshoulder:b6,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"4f0423b3000000000000504944564944,Dual Trigger 3-in-1,a:b1,b:b2,x:b0,y:b3,back:b8,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a5,lefttrigger:b6,righttrigger:b7",
	"63252305000000000000504944564944,USB Vibration Joystick (BM),x:b3,a:b2,b:b1,y:b0,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"6d0416c2000000000000504944564944,Generic DirectInput Controller,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"6d0418c2000000000000504944564944,Logitech RumblePad 2 USB,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"6d0419c2000000000000504944564944,Logitech F710 Gamepad,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"6f0e1e01000000000000504944564944,Rock Candy Gamepad for PS3,a:b1,b:b2,x:b0,y:b3,back:b8,start:b9,guide:b12,leftshoulder:b4,rightshoulder:b5,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,",
	"79000018000000000000504944564944,Mayflash WiiU Pro Game Controller Adapter (DInput),a:b1,b:b2,x:b0,y:b3,back:b8,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"79000600000000000000504944564944,G-Shark GS-GP702,a:b2,b:b1,x:b3,y:b0,back:b8,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a4,lefttrigger:b6,righttrigger:b7",
	"79004318000000000000504944564944,Mayflash GameCube Controller Adapter,a:b1,b:b2,x:b0,y:b3,back:b0,start:b9,guide:b0,leftshoulder:b4,rightshoulder:b7,leftstick:b0,rightstick:b0,leftx:a0,lefty:a1,rightx:a5,righty:a2,lefttrigger:a3,righttrigger:a4,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,",
	"79000600000000000000504944564944,Generic Speedlink,a:b2,b:b1,y:b0,x:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a4,lefttrigger:b6,righttrigger:b7,",
	"83056020000000000000504944564944,iBuffalo USB 2-axis 8-button Gamepad,a:b1,b:b0,y:b2,x:b3,start:b7,back:b6,leftshoulder:b4,rightshoulder:b5,leftx:a0,lefty:a1,",
	"88880803000000000000504944564944,PS3 Controller,a:b2,b:b1,back:b8,dpdown:h0.8,dpleft:h0.4,dpright:h0.2,dpup:h0.1,guide:b12,leftshoulder:b4,leftstick:b9,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:b7,rightx:a3,righty:a4,start:b11,x:b0,y:b3,",
	"8f0e0300000000000000504944564944,Trust GXT 28,a:b2,b:b1,y:b0,x:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"8f0e0d31000000000000504944564944,Multilaser JS071 USB,a:b1,b:b2,y:b3,x:b0,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"8f0e1200000000000000504944564944,Acme,x:b2,a:b0,b:b1,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b5,rightshoulder:b6,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a2,",
	"9000318000000000000504944564944,Mayflash Wiimote PC Adapter,a:b2,b:h0.4,x:b0,y:b1,back:b4,start:b5,guide:b11,leftshoulder:b6,rightshoulder:b3,leftx:a0,lefty:a1,",
	"a3060cff000000000000504944564944,Saitek P2500,a:b2,b:b3,y:b1,x:b0,start:b4,guide:b10,back:b5,leftstick:b8,rightstick:b9,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"c0111352000000000000504944564944,Battalife Joystick,x:b4,a:b6,b:b7,y:b5,back:b2,start:b3,leftshoulder:b0,rightshoulder:b1,leftx:a0,lefty:a1,",
	"c911f055000000000000504944564944,GAMEPAD,a:b0,b:b1,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b2,y:b3,",
	"d6206dca000000000000504944564944,PowerA Pro Ex,a:b1,b:b2,x:b0,y:b3,back:b8,guide:b12,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.0,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"ff113133000000000000504944564944,Gembird JPD-DualForce,a:b2,b:b3,x:b0,y:b1,start:b9,back:b8,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a4,lefttrigger:b6,righttrigger:b7,leftstick:b10,rightstick:b11,",
	"ff113133000000000000504944564944,SVEN X-PAD,a:b2,b:b3,y:b1,x:b0,start:b5,back:b4,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a4,lefttrigger:b8,righttrigger:b9,",
	"ffff0000000000000000504944564944,GameStop Gamepad,a:b0,b:b1,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b2,y:b3,",
	"__XINPUT_DEVICE__,XInput Gamepad,a:b12,b:b13,x:b14,y:b15,start:b4,back:b5,leftstick:b6,rightstick:b7,leftshoulder:b8,rightshoulder:b9,dpup:b0,dpdown:b1,dpleft:b2,dpright:b3,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:a4,righttrigger:a5,",
	#endif

	#ifdef OSX_ENABLED
	"0500000047532047616d657061640000,GameStop Gamepad,a:b0,b:b1,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b2,y:b3,",
	"050000005769696d6f74652028303000,Wii Remote,a:b4,b:b5,y:b9,x:b10,start:b6,guide:b8,back:b7,dpup:b2,dpleft:b0,dpdown:b3,dpright:b1,leftx:a0,lefty:a1,lefttrigger:b12,righttrigger:,leftshoulder:b11,",
	"050000005769696d6f74652028313800,Wii U Pro Controller,a:b16,b:b15,x:b18,y:b17,back:b7,guide:b8,start:b6,leftstick:b23,rightstick:b24,leftshoulder:b19,rightshoulder:b20,dpup:b11,dpdown:b12,dpleft:b13,dpright:b14,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b21,righttrigger:b22,",
	"0d0f0000000000004d00000000000000,HORI Gem Pad 3,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"0d0f0000000000006600000000000000,HORIPAD FPS PLUS 4,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a5,lefttrigger:b6,righttrigger:a4,",
	"10280000000000000900000000000000,8Bitdo SFC30 GamePad,a:b1,b:b0,x:b4,y:b3,back:b10,start:b11,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"2509000000000000e803000000000000,Mayflash Wii Classic Controller,a:b1,b:b0,x:b3,y:b2,back:b8,guide:b10,start:b9,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:b11,dpdown:b13,dpleft:b12,dpright:b14,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"351200000000000021ab000000000000,SFC30 Joystick,a:b1,b:b0,x:b4,y:b3,back:b10,start:b11,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"4c050000000000006802000000000000,PS3 Controller,a:b14,b:b13,back:b0,dpdown:b6,dpleft:b7,dpright:b5,dpup:b4,guide:b16,leftshoulder:b10,leftstick:b1,lefttrigger:b8,leftx:a0,lefty:a1,rightshoulder:b11,rightstick:b2,righttrigger:b9,rightx:a2,righty:a3,start:b3,x:b15,y:b12,",
	"4c05000000000000c405000000000000,PS4 Controller,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b12,leftshoulder:b4,leftstick:b10,lefttrigger:a3,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:a4,rightx:a2,righty:a5,start:b9,x:b0,y:b3,",
	"4f0400000000000000b3000000000000,Thrustmaster Firestorm Dual Power,a:b0,b:b2,y:b3,x:b1,start:b10,guide:b8,back:b9,leftstick:b11,rightstick:,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b5,righttrigger:b7,",
	"4f0400000000000015b3000000000000,Thrustmaster Dual Analog 3.2,x:b1,a:b0,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b5,rightshoulder:b6,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"5e040000000000008e02000000000000,X360 Controller,a:b0,b:b1,back:b9,dpdown:b12,dpleft:b13,dpright:b14,dpup:b11,guide:b10,leftshoulder:b4,leftstick:b6,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b7,righttrigger:a5,rightx:a3,righty:a4,start:b8,x:b2,y:b3,",
	"5e04000000000000dd02000000000000,Xbox One Wired Controller,x:b2,a:b0,b:b1,y:b3,back:b9,guide:b10,start:b8,dpleft:b13,dpdown:b12,dpright:b14,dpup:b11,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b6,rightstick:b7,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"6d0400000000000016c2000000000000,Logitech F310 Gamepad (DInput),a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,", /* Guide button doesn't seem to be sent in DInput mode. */
	"6d0400000000000018c2000000000000,Logitech F510 Gamepad (DInput),a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"6d0400000000000019c2000000000000,Logitech Wireless Gamepad (DInput),a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,", /* This includes F710 in DInput mode and the "Logitech Cordless RumblePad 2", at the very least. */
	"6d040000000000001fc2000000000000,Logitech F710 Gamepad (XInput),a:b0,b:b1,back:b9,dpdown:b12,dpleft:b13,dpright:b14,dpup:b11,guide:b10,leftshoulder:b4,leftstick:b6,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b7,righttrigger:a5,rightx:a3,righty:a4,start:b8,x:b2,y:b3,",
	"79000000000000000018000000000000,Mayflash WiiU Pro Game Controller Adapter (DInput),a:b4,b:b8,x:b0,y:b12,back:b32,start:b36,leftstick:b40,rightstick:b44,leftshoulder:b16,rightshoulder:b20,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a4,rightx:a8,righty:a12,lefttrigger:b24,righttrigger:b28,",
	"79000000000000000600000000000000,G-Shark GP-702,a:b2,b:b1,x:b3,y:b0,back:b8,start:b9,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:b6,righttrigger:b7,",
	"83050000000000006020000000000000,iBuffalo USB 2-axis 8-button Gamepad,a:b1,b:b0,x:b3,y:b2,back:b6,start:b7,leftshoulder:b4,rightshoulder:b5,leftx:a0,lefty:a1,",
	"891600000000000000fd000000000000,Razer Onza Tournament,a:b0,b:b1,y:b3,x:b2,start:b8,guide:b10,back:b9,leftstick:b6,rightstick:b7,leftshoulder:b4,rightshoulder:b5,dpup:b11,dpleft:b13,dpdown:b12,dpright:b14,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,",
	"8f0e0000000000000300000000000000,Piranha xtreme,x:b3,a:b2,b:b1,y:b0,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b6,lefttrigger:b4,rightshoulder:b7,righttrigger:b5,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a2,",
	"ad1b00000000000001f9000000000000,Gamestop BB-070 X360 Controller,a:b0,b:b1,back:b9,dpdown:b12,dpleft:b13,dpright:b14,dpup:b11,guide:b10,leftshoulder:b4,leftstick:b6,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b7,righttrigger:a5,rightx:a3,righty:a4,start:b8,x:b2,y:b3,",
	"b4040000000000000a01000000000000,Sega Saturn USB Gamepad,a:b0,b:b1,x:b3,y:b4,back:b5,guide:b2,start:b8,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"d814000000000000cecf000000000000,MC Cthulhu,leftx:,lefty:,rightx:,righty:,lefttrigger:b6,a:b1,b:b2,y:b3,x:b0,start:b9,back:b8,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,righttrigger:b7,",
	#endif

	#if X11_ENABLED
	"0000000058626f782047616d65706100,Xbox Gamepad (userspace driver),a:b0,b:b1,x:b2,y:b3,start:b7,back:b6,guide:b8,dpup:h0.1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,leftshoulder:b4,rightshoulder:b5,lefttrigger:a5,righttrigger:a4,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"0300000000f000000300000000010000,RetroUSB.com RetroPad,a:b1,b:b5,x:b0,y:b4,back:b2,start:b3,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"0300000000f00000f100000000010000,RetroUSB.com Super RetroPort,a:b1,b:b5,x:b0,y:b4,back:b2,start:b3,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"030000000d0f00001000000011010000,HORI CO.,LTD. FIGHTING STICK 3,x:b0,a:b1,b:b2,y:b3,back:b8,guide:b12,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7",
	"030000000d0f00002200000011010000,HORI CO.,LTD. REAL ARCADE Pro.V3,x:b0,a:b1,b:b2,y:b3,back:b8,guide:b12,start:b9,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,",
	"030000000d0f00004d00000011010000,HORI Gem Pad 3,x:b0,a:b1,b:b2,y:b3,back:b8,guide:b12,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"03000000100800000100000010010000,Twin USB PS2 Adapter,a:b2,b:b1,y:b0,x:b3,start:b9,guide:,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a2,lefttrigger:b4,righttrigger:b5,",
	"03000000250900000500000000010000,Sony PS2 pad with SmartJoy adapter,a:b2,b:b1,y:b0,x:b3,start:b8,back:b9,leftstick:b10,rightstick:b11,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b4,righttrigger:b5,",
	"03000000260900008888000000010000,GameCube {WiseGroup USB box},a:b0,b:b2,y:b3,x:b1,start:b7,leftshoulder:,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,rightstick:,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:a4,righttrigger:a5,",
	"03000000280400000140000000010000,Gravis GamePad Pro USB ,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftx:a0,lefty:a1,",
	"03000000341a000005f7000010010000,GameCube {HuiJia USB box},a:b1,b:b2,y:b3,x:b0,start:b9,guide:,back:,leftstick:,rightstick:,leftshoulder:,dpleft:b15,dpdown:b14,dpright:b13,leftx:a0,lefty:a1,rightx:a5,righty:a2,lefttrigger:a3,righttrigger:a4,rightshoulder:b7,dpup:b12,",
	"03000000380700001647000010040000,Mad Catz Wired Xbox 360 Controller,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000004c0500006802000011010000,PS3 Controller,a:b14,b:b13,back:b0,dpdown:b6,dpleft:b7,dpright:b5,dpup:b4,guide:b16,leftshoulder:b10,leftstick:b1,lefttrigger:b8,leftx:a0,lefty:a1,rightshoulder:b11,rightstick:b2,righttrigger:b9,rightx:a2,righty:a3,start:b3,x:b15,y:b12,",
	"030000004c050000c405000011010000,Sony DualShock 4,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a5,lefttrigger:a3,righttrigger:a4,",
	"030000004f04000000b3000010010000,Thrustmaster Firestorm Dual Power,a:b0,b:b2,y:b3,x:b1,start:b10,guide:b8,back:b9,leftstick:b11,rightstick:b12,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b5,righttrigger:b7,",
	"030000004f04000008d0000000010000,Thrustmaster Run N Drive  Wireless,a:b1,b:b2,x:b0,y:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a5,lefttrigger:b6,righttrigger:b7,",
	"030000004f04000009d0000000010000,Thrustmaster Run N Drive Wireless PS3,a:b1,b:b2,x:b0,y:b3,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"030000004f04000015b3000010010000,Thrustmaster Dual Analog 4,a:b0,b:b2,x:b1,y:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b5,righttrigger:b7,",
	"030000004f04000020b3000010010000,Thrustmaster 2 in 1 DT,a:b0,b:b2,y:b3,x:b1,start:b9,guide:,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b6,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b5,righttrigger:b7,",
	"030000004f04000023b3000000010000,Thrustmaster Dual Trigger 3-in-1,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a5,",
	"030000005e0400001907000000010000,X360 Wireless Controller,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b10,rightshoulder:b5,rightx:a3,start:b7,righty:a4,dpleft:h0.8,lefttrigger:a2,x:b2,dpup:h0.1,back:b6,leftstick:b9,leftshoulder:b4,y:b3,a:b0,dpright:h0.2,righttrigger:a5,b:b1,",
	"030000005e0400008902000021010000,Microsoft X-Box pad v2 (US),x:b3,a:b0,b:b1,y:b4,back:b6,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b5,lefttrigger:a2,rightshoulder:b2,righttrigger:a5,leftstick:b8,rightstick:b9,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000005e0400008e02000001000000,Microsoft X-Box 360 pad,leftstick:b9,leftx:a0,lefty:a1,dpdown:h0.1,rightstick:b10,rightshoulder:b5,rightx:a3,start:b7,righty:a4,dpleft:h0.2,lefttrigger:a2,x:b2,dpup:h0.4,back:b6,leftshoulder:b4,y:b3,a:b0,dpright:h0.8,righttrigger:a5,b:b1,",
	"030000005e0400008e02000004010000,Microsoft X-Box 360 pad,a:b0,b:b1,x:b2,y:b3,back:b6,start:b7,guide:b8,leftshoulder:b4,rightshoulder:b5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,",
	"030000005e0400008e02000010010000,X360 Controller,a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"030000005e0400008e02000014010000,X360 Controller,a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"030000005e0400008e02000020200000,SpeedLink XEOX Pro Analog Gamepad pad,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000005e0400008e02000062230000,Microsoft X-Box 360 pad,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000005e0400009102000007010000,X360 Wireless Controller,a:b0,b:b1,y:b3,x:b2,start:b7,guide:b8,back:b6,leftstick:b9,rightstick:b10,leftshoulder:b4,rightshoulder:b5,dpup:b13,dpleft:b11,dpdown:b14,dpright:b12,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,",
	"030000005e040000d102000001010000,Microsoft X-Box One pad,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000005e040000dd02000003020000,Microsoft X-Box One pad v2,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"03000000666600000488000000010000,Super Joy Box 5 Pro,a:b2,b:b1,x:b3,y:b0,back:b9,start:b8,leftshoulder:b6,rightshoulder:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b4,righttrigger:b5,dpup:b12,dpleft:b15,dpdown:b14,dpright:b13,",
	"030000006d04000011c2000010010000,Logitech WingMan Cordless RumblePad,a:b0,b:b1,y:b4,x:b3,start:b8,guide:b5,back:b2,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:b9,righttrigger:b10,",
	"030000006d04000016c2000010010000,Logitech Logitech Dual Action,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"030000006d04000018c2000010010000,Logitech Logitech RumblePad 2 USB,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"030000006d04000019c2000010010000,Logitech Cordless RumblePad 2,a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"030000006d04000019c2000011010000,Logitech F710 Gamepad (DInput),a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b0,y:b3,",
	"030000006d0400001dc2000014400000,Logitech F310 Gamepad (XInput),a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"030000006d0400001ec2000020200000,Logitech F510 Gamepad (XInput),a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"030000006d0400001fc2000005030000,Logitech F710 Gamepad (XInput),a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"030000006e0500000320000010010000,JC-U3613M - DirectInput Mode,x:b0,a:b2,b:b3,y:b1,back:b10,guide:b12,start:b11,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b8,rightstick:b9,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"030000006f0e00000103000000020000,Logic3 Controller,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000006f0e00001304000000010000,Generic X-Box pad,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000006f0e00001e01000011010000,Rock Candy Gamepad for PS3,a:b1,b:b2,x:b0,y:b3,back:b8,start:b9,guide:b12,leftshoulder:b4,rightshoulder:b5,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,",
	"030000006f0e00001f01000000010000,Generic X-Box pad,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000006f0e00002801000011010000,PDP Rock Candy Wireless Controller for PS3,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b11,rightshoulder:b5,rightx:a2,start:b9,righty:a3,dpleft:h0.8,lefttrigger:b6,x:b0,dpup:h0.1,back:b8,leftstick:b10,leftshoulder:b4,y:b3,a:b1,dpright:h0.2,righttrigger:b7,b:b2,",
	"030000006f0e00003001000001010000,EA Sports PS3 Controller,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"030000006f0e00003901000020060000,Afterglow Wired Controller for Xbox One,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000006f0e00004601000001010000,Rock Candy Wired Controller for Xbox One,a:b0,b:b1,x:b2,y:b3,leftshoulder:b4,rightshoulder:b5,back:b6,start:b7,guide:b8,leftstick:b9,rightstick:b10,lefttrigger:a2,righttrigger:a5,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"03000000790000000600000010010000,DragonRise Inc. Generic USB Joystick,x:b3,a:b2,b:b1,y:b0,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"03000000790000001100000010010000,RetroLink Saturn Classic Controller,x:b3,a:b0,b:b1,y:b4,back:b5,guide:b2,start:b8,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"03000000830500006020000010010000,iBuffalo USB 2-axis 8-button Gamepad,a:b1,b:b0,x:b3,y:b2,back:b6,start:b7,leftshoulder:b4,rightshoulder:b5,leftx:a0,lefty:a1,",
	"030000008916000000fd000024010000,Razer Onza Tournament,a:b0,b:b1,y:b3,x:b2,start:b7,guide:b8,back:b6,leftstick:b9,rightstick:b10,leftshoulder:b4,rightshoulder:b5,dpup:b13,dpleft:b11,dpdown:b14,dpright:b12,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,",
	"030000008916000001fd000024010000,Razer Onza Classic Edition,x:b2,a:b0,b:b1,y:b3,back:b6,guide:b8,start:b7,dpleft:b11,dpdown:b14,dpright:b12,dpup:b13,leftshoulder:b4,lefttrigger:a2,rightshoulder:b5,righttrigger:a5,leftstick:b9,rightstick:b10,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"030000008f0e00000300000010010000,GreenAsia Inc. USB Joystick,x:b3,a:b2,b:b1,y:b0,back:b8,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b6,lefttrigger:b4,rightshoulder:b7,righttrigger:b5,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a2,",
	"030000008f0e00001200000010010000,GreenAsia Inc. USB Joystick,x:b2,a:b0,b:b1,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b5,rightshoulder:b6,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a2,",
	"03000000a30600000901000000010000,Saitek P880,a:b2,b:b3,y:b1,x:b0,leftstick:b8,rightstick:b9,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a2,lefttrigger:b6,righttrigger:b7,",
	"03000000a30600000c04000011010000,Saitek P2900 Wireless Pad,a:b1,b:b2,y:b3,x:b0,start:b12,guide:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b6,rightshoulder:b7,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a2,lefttrigger:b4,righttrigger:b5,",
	"03000000a306000018f5000010010000,Saitek PLC Saitek P3200 Rumble Pad,x:b0,a:b1,b:b2,y:b3,back:b8,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:a2,rightshoulder:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a3,righty:a4,",
	"03000000a306000023f6000011010000,Saitek Cyborg V.1 Game Pad,a:b1,b:b2,y:b3,x:b0,start:b9,guide:b12,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a4,lefttrigger:b6,righttrigger:b7,",
	"03000000ad1b000001f5000033050000,Hori Pad EX Turbo 2,a:b0,b:b1,y:b3,x:b2,start:b7,guide:b8,back:b6,leftstick:b9,rightstick:b10,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,",
	"03000000ad1b000016f0000090040000,Mad Catz Xbox 360 Controller,a:b0,b:b1,y:b3,x:b2,start:b7,guide:b8,back:b6,leftstick:b9,rightstick:b10,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a3,righty:a4,lefttrigger:a2,righttrigger:a5,",
	"03000000ba2200002010000001010000,Jess Technology USB Game Controller,a:b2,b:b1,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b4,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,righttrigger:b7,rightx:a3,righty:a2,start:b9,x:b3,y:b0,",
	"03000000c9110000f055000011010000,HJC Game GAMEPAD,platform:Linux,x:b2,a:b0,b:b1,y:b3,back:b4,back:b8,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"03000000d814000007cd000011010000,Toodles 2008 Chimp PC/PS3,a:b0,b:b1,y:b2,x:b3,start:b9,back:b8,leftshoulder:b4,rightshoulder:b5,leftx:a0,lefty:a1,lefttrigger:b6,righttrigger:b7,",
	"03000000d81400000862000011010000,HitBox (PS3/PC) Analog Mode,a:b1,b:b2,y:b3,x:b0,start:b12,guide:b9,back:b8,leftshoulder:b4,rightshoulder:b5,lefttrigger:b6,righttrigger:b7,leftx:a0,lefty:a1,",
	"03000000de280000ff11000001000000,Valve Streaming Gamepad,a:b0,b:b1,back:b6,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b8,leftshoulder:b4,leftstick:b9,lefttrigger:a2,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b10,righttrigger:a5,rightx:a3,righty:a4,start:b7,x:b2,y:b3,",
	"03000000f0250000c183000010010000,Goodbetterbest Ltd USB Controller,x:b0,a:b1,b:b2,y:b3,back:b8,guide:b12,start:b9,dpleft:h0.8,dpdown:h0.0,dpdown:h0.4,dpright:h0.0,dpright:h0.2,dpup:h0.0,dpup:h0.1,leftshoulder:h0.0,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"03000000fd0500002a26000000010000,3dfx InterAct HammerHead FX,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b5,rightshoulder:b7,rightx:a2,start:b11,righty:a3,dpleft:h0.8,lefttrigger:b8,x:b0,dpup:h0.1,back:b10,leftstick:b2,leftshoulder:b6,y:b1,a:b3,dpright:h0.2,righttrigger:b9,b:b4,",
	"03000000ff1100003133000010010000,PC Game Controller,a:b2,b:b1,y:b0,x:b3,start:b9,back:b8,leftstick:b10,rightstick:b11,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"05000000010000000100000003000000,Nintendo Wiimote,a:b0,b:b1,y:b3,x:b2,start:b9,guide:b10,back:b8,leftstick:b11,rightstick:b12,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,",
	"05000000102800000900000000010000,8Bitdo SFC30 GamePad,x:b4,a:b1,b:b0,y:b3,back:b10,start:b11,leftshoulder:b6,rightshoulder:b7,leftx:a0,lefty:a1,",
	"05000000362800000100000002010000,OUYA Game Controller,leftx:a0,lefty:a1,dpdown:b9,rightstick:b7,rightshoulder:b5,rightx:a3,start:b16,righty:a4,dpleft:b10,lefttrigger:b12,x:b1,dpup:b8,back:b14,leftstick:b6,leftshoulder:b4,y:b2,a:b0,dpright:b11,righttrigger:b13,b:b3,",
	"05000000362800000100000003010000,OUYA Game Controller,leftx:a0,lefty:a1,dpdown:b9,rightstick:b7,rightshoulder:b5,rightx:a3,start:b16,righty:a4,dpleft:b10,lefttrigger:b12,x:b1,dpup:b8,back:b14,leftstick:b6,leftshoulder:b4,y:b2,a:b0,dpright:b11,righttrigger:b13,b:b3,",
	"05000000362800000100000004010000,OUYA Game Controller,leftx:a0,lefty:a1,dpdown:b9,rightstick:b7,rightshoulder:b5,rightx:a3,start:b16,righty:a4,dpleft:b10,lefttrigger:b12,x:b1,dpup:b8,back:b14,leftstick:b6,leftshoulder:b4,y:b2,a:b0,dpright:b11,righttrigger:b13,b:b3,",
	"05000000380700006652000025010000,Mad Catz C.T.R.L.R,x:b0,a:b1,b:b2,y:b3,back:b8,guide:b12,start:b9,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,dpup:h0.1,leftshoulder:b4,lefttrigger:b6,rightshoulder:b5,righttrigger:b7,leftstick:b10,rightstick:b11,leftx:a0,lefty:a1,rightx:a2,righty:a3,",
	"0500000047532047616d657061640000,GameStop Gamepad,a:b0,b:b1,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:,leftshoulder:b4,leftstick:b10,lefttrigger:b6,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:b7,rightx:a2,righty:a3,start:b9,x:b2,y:b3,",
	"050000004c0500006802000000010000,PS3 Controller (Bluetooth),a:b14,b:b13,y:b12,x:b15,start:b3,guide:b16,back:b0,leftstick:b1,rightstick:b2,leftshoulder:b10,rightshoulder:b11,dpup:b4,dpleft:b7,dpdown:b6,dpright:b5,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b8,righttrigger:b9,",
	"050000004c050000c405000000010000,PS4 Controller (Bluetooth),a:b1,b:b2,back:b8,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,guide:b12,leftshoulder:b4,leftstick:b10,lefttrigger:a3,leftx:a0,lefty:a1,rightshoulder:b5,rightstick:b11,righttrigger:a4,rightx:a2,righty:a5,start:b9,x:b0,y:b3,",
	"050000007e0500003003000001000000,Nintendo Wii U Pro Controller,a:b0,b:b1,x:b3,y:b2,back:b8,start:b9,guide:b10,leftshoulder:b4,rightshoulder:b5,leftstick:b11,rightstick:b12,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b6,righttrigger:b7,dpup:b13,dpleft:b15,dpdown:b14,dpright:b16,",
	"05000000d6200000ad0d000001000000,Moga Pro,a:b0,b:b1,y:b3,x:b2,start:b6,leftstick:b7,rightstick:b8,leftshoulder:b4,rightshoulder:b5,dpup:h0.1,dpleft:h0.8,dpdown:h0.4,dpright:h0.2,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:a5,righttrigger:a4,",
	"060000004c0500006802000000010000,PS3 Controller (Bluetooth),a:b14,b:b13,y:b12,x:b15,start:b3,guide:b16,back:b0,leftstick:b1,rightstick:b2,leftshoulder:b10,rightshoulder:b11,dpup:b4,dpleft:b7,dpdown:b6,dpright:b5,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:b8,righttrigger:b9,",
	#endif

	#if defined(__ANDROID__)
	"Default Android Gamepad,Default Controller,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b8,rightshoulder:b10,rightx:a2,start:b6,righty:a3,dpleft:h0.8,lefttrigger:a4,x:b2,dpup:h0.1,back:b4,leftstick:b7,leftshoulder:b9,y:b3,a:b0,dpright:h0.2,righttrigger:a5,b:b1,",
	"47656e6572696320582d426f78207061,Logitech F-310,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b8,rightshoulder:b10,rightx:a2,start:b6,righty:a3,dpleft:h0.8,lefttrigger:a5,x:b2,dpup:h0.1,leftstick:b7,leftshoulder:b9,y:b3,a:b0,dpright:h0.2,righttrigger:a4,b:b1,",
	"484f524920434f2e2c4c544420205041,Hori Gem Pad 3,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b6,rightshoulder:b18,rightx:a2,start:b16,righty:a3,dpleft:h0.8,lefttrigger:b9,x:b0,dpup:h0.1,back:b15,leftstick:b4,leftshoulder:b3,y:b2,a:b1,dpright:h0.2,righttrigger:b10,b:b17,",
	"4d6963726f736f667420582d426f7820,Microsoft X-Box 360 pad,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b8,rightshoulder:b10,rightx:a2,start:b6,righty:a3,dpleft:h0.8,lefttrigger:a4,x:b2,dpup:h0.1,leftstick:b7,leftshoulder:b9,y:b3,a:b0,dpright:h0.2,righttrigger:a5,b:b1,",
	"4e564944494120436f72706f72617469,NVIDIA Controller,a:b0,b:b1,dpdown:h0.4,dpleft:h0.8,dpright:h0.2,dpup:h0.1,leftshoulder:b9,leftstick:b7,lefttrigger:a4,leftx:a0,lefty:a1,rightshoulder:b10,rightstick:b8,righttrigger:a5,rightx:a2,righty:a3,start:b6,x:b2,y:b3,",
	"532e542e442e20496e74657261637420,3dfx InterAct HammerHead FX,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b25,rightshoulder:b27,rightx:a2,start:b31,righty:a3,dpleft:h0.8,lefttrigger:b28,x:b20,dpup:h0.1,back:b30,leftstick:b22,leftshoulder:b26,y:b21,a:b23,dpright:h0.2,righttrigger:b29,b:b24,",
	"506572666f726d616e63652044657369,PDP Rock Candy Wireless Controller for PS3,leftx:a0,lefty:a1,dpdown:h0.4,rightstick:b6,rightshoulder:b18,rightx:a2,start:b16,righty:a3,dpleft:h0.8,lefttrigger:b9,x:b0,dpup:h0.1,back:h0.2,leftstick:b4,leftshoulder:b3,y:b2,a:b1,dpright:h0.2,righttrigger:b10,b:b17,",
	#endif

	#ifdef JAVASCRIPT_ENABLED
	"Default HTML5 Gamepad, Default Mapping,leftx:a0,lefty:a1,dpdown:b13,rightstick:b11,rightshoulder:b5,rightx:a2,start:b9,righty:a3,dpleft:b14,lefttrigger:a6,x:b2,dpup:b12,back:b8,leftstick:b10,leftshoulder:b4,y:b3,a:b0,dpright:b15,righttrigger:a7,b:b1,",
	"303435652d303238652d4d6963726f73,Wired X360 Controller,leftx:a0,lefty:a1,dpdown:a7,rightstick:b10,rightshoulder:b5,rightx:a3,start:b7,righty:a4,dpleft:a6,lefttrigger:a2,x:b2,dpup:a7,back:b6,leftstick:b9,leftshoulder:b4,y:b3,a:b0,dpright:a6,righttrigger:a5,b:b1,",
	"303435652d303731392d58626f782033,Wireless X360 Controller,leftx:a0,lefty:a1,dpdown:b14,rightstick:b10,rightshoulder:b5,rightx:a3,start:b7,righty:a4,dpleft:b11,lefttrigger:a2,x:b2,dpup:b13,back:b6,leftstick:b9,leftshoulder:b4,y:b3,a:b0,dpright:b12,righttrigger:a5,b:b1,",
	"303534632d303236382d536f6e792050,PS3 Controller USB/Linux,leftx:a0,lefty:a1,dpdown:b6,rightstick:b2,rightshoulder:b11,rightx:a2,start:b3,righty:a3,dpleft:b7,lefttrigger:b8,x:b15,dpup:b4,back:b0,leftstick:b1,leftshoulder:b10,y:b12,a:b14,dpright:b5,righttrigger:b9,b:b13,",
	"303534632d303563342d536f6e792043,PS4 Controller USB/Linux,leftx:a0,lefty:a1,dpdown:a7,rightstick:b11,rightshoulder:b5,rightx:a2,start:b9,righty:a5,dpleft:a6,lefttrigger:a3,x:b0,dpup:a7,back:b8,leftstick:b10,leftshoulder:b4,y:b3,a:b1,dpright:a6,righttrigger:a4,b:b2,",
	"303534632d303563342d576972656c65,PS4 Controller USB/Win,leftx:a0,lefty:a1,dpdown:b15,rightstick:b11,rightshoulder:b5,rightx:a2,start:b9,righty:a5,lefttrigger:a3,x:b0,dpup:b14,dpleft:b16,dpright:b17,back:b8,leftstick:b10,leftshoulder:b4,y:b3,a:b1,righttrigger:b7,b:b2,",
	"c2a94d6963726f736f66742058626f78,Wireless X360 Controller,leftx:a0,lefty:a1,dpdown:b14,rightstick:b10,rightshoulder:b5,rightx:a3,start:b7,righty:a4,dpleft:b11,lefttrigger:a2,x:b2,dpup:b13,back:b6,leftstick:b9,leftshoulder:b4,y:b3,a:b0,dpright:b12,righttrigger:a5,b:b1,",
	#endif

	#ifdef UWP_ENABLED
	"__UWP_GAMEPAD__,Xbox Controller,a:b2,b:b3,x:b4,y:b5,start:b0,back:b1,leftstick:b12,rightstick:b13,leftshoulder:b10,rightshoulder:b11,dpup:b6,dpdown:b7,dpleft:b8,dpright:b9,leftx:a0,lefty:a1,rightx:a2,righty:a3,lefttrigger:a4,righttrigger:a5,",
	#endif
	NULL
};

InputDefault::InputDefault() {

	mouse_button_mask=0;
	emulate_touch=false;
	main_loop=NULL;

	hat_map_default[HAT_UP].type = TYPE_BUTTON;
	hat_map_default[HAT_UP].index = JOY_DPAD_UP;
	hat_map_default[HAT_UP].value = 0;

	hat_map_default[HAT_RIGHT].type = TYPE_BUTTON;
	hat_map_default[HAT_RIGHT].index = JOY_DPAD_RIGHT;
	hat_map_default[HAT_RIGHT].value = 0;

	hat_map_default[HAT_DOWN].type = TYPE_BUTTON;
	hat_map_default[HAT_DOWN].index = JOY_DPAD_DOWN;
	hat_map_default[HAT_DOWN].value = 0;

	hat_map_default[HAT_LEFT].type = TYPE_BUTTON;
	hat_map_default[HAT_LEFT].index = JOY_DPAD_LEFT;
	hat_map_default[HAT_LEFT].value = 0;

	fallback_mapping = -1;

	String env_mapping = OS::get_singleton()->get_environment("SDL_GAMECONTROLLERCONFIG");
	if (env_mapping != "") {

		Vector<String> entries = env_mapping.split("\n");
		for (int i=0; i < entries.size(); i++) {
			if (entries[i] == "")
				continue;
			parse_mapping(entries[i]);
		};
	};

	int i = 0;
	while (s_ControllerMappings[i]) {

		parse_mapping(s_ControllerMappings[i++]);
	};
}


uint32_t InputDefault::joy_button(uint32_t p_last_id, int p_device, int p_button, bool p_pressed) {

	_THREAD_SAFE_METHOD_;
	Joypad& joy = joy_names[p_device];
	//printf("got button %i, mapping is %i\n", p_button, joy.mapping);
	if (joy.last_buttons[p_button] == p_pressed) {
		return p_last_id;
		//printf("same button value\n");
	}
	joy.last_buttons[p_button] = p_pressed;
	if (joy.mapping == -1) {
		return _button_event(p_last_id, p_device, p_button, p_pressed);
	};

	Map<int,JoyEvent>::Element* el = map_db[joy.mapping].buttons.find(p_button);
	if (!el) {
		//don't process un-mapped events for now, it could mess things up badly for devices with additional buttons/axis
		//return _button_event(p_last_id, p_device, p_button, p_pressed);
		return p_last_id;
	};

	JoyEvent map = el->get();
	if (map.type == TYPE_BUTTON) {
		//fake additional axis event for triggers
		if (map.index == JOY_L2 || map.index == JOY_R2) {
			float value = p_pressed ? 1.0f : 0.0f;
			int axis = map.index == JOY_L2 ? JOY_ANALOG_L2 : JOY_ANALOG_R2;
			p_last_id = _axis_event(p_last_id, p_device, axis, value);
		}
		return _button_event(p_last_id, p_device, map.index, p_pressed);
	};

	if (map.type == TYPE_AXIS) {
		return _axis_event(p_last_id, p_device, map.index, p_pressed ? 1.0 : 0.0);
	};

	return p_last_id; // no event?
};

uint32_t InputDefault::joy_axis(uint32_t p_last_id, int p_device, int p_axis, const JoyAxis& p_value) {

	_THREAD_SAFE_METHOD_;

	Joypad& joy = joy_names[p_device];

	if (joy.last_axis[p_axis] == p_value.value) {
		return p_last_id;
	}

	if (p_value.value > joy.last_axis[p_axis]) {

		if (p_value.value < joy.last_axis[p_axis] + joy.filter ) {

			return p_last_id;
		}
	}
	else if (p_value.value > joy.last_axis[p_axis] - joy.filter) {

		return p_last_id;
	}

	if (ABS(joy.last_axis[p_axis]) > 0.5 && joy.last_axis[p_axis] * p_value.value < 0) {
		//changed direction quickly, insert fake event to release pending inputmap actions
		JoyAxis jx;
		jx.min = p_value.min;
		jx.value = p_value.value < 0 ? 0.1 : -0.1;
		p_last_id = joy_axis(p_last_id, p_device, p_axis, jx);
	}

	joy.last_axis[p_axis] = p_value.value;
	float val = p_value.min == 0 ? -1.0f + 2.0f * p_value.value : p_value.value;

	if (joy.mapping == -1) {
		return _axis_event(p_last_id, p_device, p_axis, val);
	};

	Map<int,JoyEvent>::Element* el = map_db[joy.mapping].axis.find(p_axis);
	if (!el) {
		//return _axis_event(p_last_id, p_device, p_axis, p_value);
		return p_last_id;
	};


	JoyEvent map = el->get();

	if (map.type == TYPE_BUTTON) {
		//send axis event for triggers
		if (map.index == JOY_L2 || map.index == JOY_R2) {
			float value = p_value.min == 0 ? p_value.value : 0.5f + p_value.value / 2.0f;
			int axis = map.index == JOY_L2 ? JOY_ANALOG_L2 : JOY_ANALOG_R2;
			p_last_id = _axis_event(p_last_id, p_device, axis, value);
		}

		if (map.index == JOY_DPAD_UP || map.index == JOY_DPAD_DOWN) {
			bool pressed = p_value.value != 0.0f;
			int button = p_value.value < 0 ? JOY_DPAD_UP : JOY_DPAD_DOWN;

			if (!pressed) {
				if (joy_buttons_pressed.has(_combine_device(JOY_DPAD_UP, p_device))) {
					p_last_id = _button_event(p_last_id, p_device, JOY_DPAD_UP, false);
				}
				if (joy_buttons_pressed.has(_combine_device(JOY_DPAD_DOWN, p_device))) {
					p_last_id = _button_event(p_last_id, p_device, JOY_DPAD_DOWN, false);
				}
			}
			if ( pressed == joy_buttons_pressed.has(_combine_device(button, p_device))) {
				return p_last_id;
			}
			return _button_event(p_last_id, p_device, button, true);
		}
		if (map.index == JOY_DPAD_LEFT || map.index == JOY_DPAD_RIGHT) {
			bool pressed = p_value.value != 0.0f;
			int button = p_value.value < 0 ? JOY_DPAD_LEFT : JOY_DPAD_RIGHT;

			if (!pressed) {
				if (joy_buttons_pressed.has(_combine_device(JOY_DPAD_LEFT, p_device))) {
					p_last_id = _button_event(p_last_id, p_device, JOY_DPAD_LEFT, false);
				}
				if (joy_buttons_pressed.has(_combine_device(JOY_DPAD_RIGHT, p_device))) {
					p_last_id = _button_event(p_last_id, p_device, JOY_DPAD_RIGHT, false);
				}
			}
			if ( pressed == joy_buttons_pressed.has(_combine_device(button, p_device))) {
				return p_last_id;
			}
			return _button_event(p_last_id, p_device, button, true);
		}
		float deadzone = p_value.min == 0 ? 0.5f : 0.0f;
		bool pressed = p_value.value > deadzone ? true : false;
		if (pressed == joy_buttons_pressed.has(_combine_device(map.index,p_device))) {
			// button already pressed or released, this is an axis bounce value
			return p_last_id;
		};
		return _button_event(p_last_id, p_device, map.index, pressed);
	};

	if (map.type == TYPE_AXIS) {

		return _axis_event(p_last_id, p_device, map.index, val );
	};
	//printf("invalid mapping\n");
	return p_last_id;
};

uint32_t InputDefault::joy_hat(uint32_t p_last_id, int p_device, int p_val) {

	_THREAD_SAFE_METHOD_;
	const Joypad& joy = joy_names[p_device];

	JoyEvent* map;

	if (joy.mapping == -1) {
		map = hat_map_default;
	} else {
		map = map_db[joy.mapping].hat;
	};

	int cur_val = joy_names[p_device].hat_current;

	if ( (p_val & HAT_MASK_UP) != (cur_val & HAT_MASK_UP) ) {
		p_last_id = _button_event(p_last_id, p_device, map[HAT_UP].index, p_val & HAT_MASK_UP);
	};

	if ( (p_val & HAT_MASK_RIGHT) != (cur_val & HAT_MASK_RIGHT) ) {
		p_last_id = _button_event(p_last_id, p_device, map[HAT_RIGHT].index, p_val & HAT_MASK_RIGHT);
	};
	if ( (p_val & HAT_MASK_DOWN) != (cur_val & HAT_MASK_DOWN) ) {
		p_last_id = _button_event(p_last_id, p_device, map[HAT_DOWN].index, p_val & HAT_MASK_DOWN);
	};
	if ( (p_val & HAT_MASK_LEFT) != (cur_val & HAT_MASK_LEFT) ) {
		p_last_id = _button_event(p_last_id, p_device, map[HAT_LEFT].index, p_val & HAT_MASK_LEFT);
	};

	joy_names[p_device].hat_current = p_val;

	return p_last_id;
};

uint32_t InputDefault::_button_event(uint32_t p_last_id, int p_device, int p_index, bool p_pressed) {

	InputEvent ievent;
	ievent.type = InputEvent::JOYPAD_BUTTON;
	ievent.device = p_device;
	ievent.ID = ++p_last_id;
	ievent.joy_button.button_index = p_index;
	ievent.joy_button.pressed = p_pressed;

	parse_input_event(ievent);

	return p_last_id;
};

uint32_t InputDefault::_axis_event(uint32_t p_last_id, int p_device, int p_axis, float p_value) {

	InputEvent ievent;
	ievent.type = InputEvent::JOYPAD_MOTION;
	ievent.device = p_device;
	ievent.ID = ++p_last_id;
	ievent.joy_motion.axis = p_axis;
	ievent.joy_motion.axis_value = p_value;

	parse_input_event( ievent );

	return p_last_id;
};

InputDefault::JoyEvent InputDefault::_find_to_event(String p_to) {

	// string names of the SDL buttons in the same order as input_event.h godot buttons
	static const char* buttons[] = {"a", "b", "x", "y", "leftshoulder", "rightshoulder", "lefttrigger", "righttrigger", "leftstick", "rightstick", "back", "start", "dpup", "dpdown", "dpleft", "dpright", "guide", NULL };

	static const char* axis[] = {"leftx", "lefty", "rightx", "righty", NULL };

	JoyEvent ret;
	ret.type = -1;

	int i=0;
	while (buttons[i]) {

		if (p_to == buttons[i]) {
			//printf("mapping button %s\n", buttons[i]);
			ret.type = TYPE_BUTTON;
			ret.index = i;
			ret.value = 0;
			return ret;
		};
		++i;
	};

	i = 0;
	while (axis[i]) {

		if (p_to == axis[i]) {
			ret.type = TYPE_AXIS;
			ret.index = i;
			ret.value = 0;
			return ret;
		};
		++i;
	};

	return ret;
};

void InputDefault::parse_mapping(String p_mapping) {

	_THREAD_SAFE_METHOD_;
	JoyDeviceMapping mapping;

	Vector<String> entry = p_mapping.split(",");
	CharString uid;
	uid.resize(17);

	mapping.uid = entry[0];
	mapping.name = entry[1];

	int idx = 1;
	while (++idx < entry.size()) {

		if (entry[idx] == "")
			continue;

		String from = entry[idx].get_slice(":", 1).replace(" ", "");
		String to = entry[idx].get_slice(":", 0).replace(" ", "");

		JoyEvent to_event = _find_to_event(to);
		if (to_event.type == -1)
			continue;

		String etype = from.substr(0, 1);
		if (etype == "a") {

			int aid = from.substr(1, from.length()-1).to_int();
			mapping.axis[aid] = to_event;

		} else if (etype == "b") {

			int bid = from.substr(1, from.length()-1).to_int();
			mapping.buttons[bid] = to_event;

		} else if (etype == "h") {

			int hat_value = from.get_slice(".", 1).to_int();
			switch (hat_value) {
			case 1:
				mapping.hat[HAT_UP] = to_event;
				break;
			case 2:
				mapping.hat[HAT_RIGHT] = to_event;
				break;
			case 4:
				mapping.hat[HAT_DOWN] = to_event;
				break;
			case 8:
				mapping.hat[HAT_LEFT] = to_event;
				break;
			};
		};
	};
	map_db.push_back(mapping);
	//printf("added mapping with uuid %ls\n", mapping.uid.c_str());
};

void InputDefault::add_joy_mapping(String p_mapping, bool p_update_existing) {
	parse_mapping(p_mapping);
	if (p_update_existing) {
		Vector<String> entry = p_mapping.split(",");
		String uid = entry[0];
		for (int i=0; i<joy_names.size(); i++) {
			if (uid == joy_names[i].uid) {
				joy_names[i].mapping = map_db.size() -1;
			}
		}
	}
}

void InputDefault::remove_joy_mapping(String p_guid) {
	for (int i=map_db.size()-1; i >= 0;i--) {
		if (p_guid == map_db[i].uid) {
			map_db.remove(i);
		}
	}
	for (int i=0; i<joy_names.size(); i++) {
		if (joy_names[i].uid == p_guid) {
			joy_names[i].mapping = -1;
		}
	}
}

void InputDefault::set_fallback_mapping(String p_guid) {

	for (int i = 0; i < map_db.size(); i++) {
		if (map_db[i].uid == p_guid) {
			fallback_mapping = i;
			return;
		}
	}
}

//Defaults to simple implementation for platforms with a fixed gamepad layout, like consoles.
bool InputDefault::is_joy_known(int p_device) {

	return OS::get_singleton()->is_joy_known(p_device);
}

String InputDefault::get_joy_guid(int p_device) const {
	return OS::get_singleton()->get_joy_guid(p_device);
}

//platforms that use the remapping system can override and call to these ones
bool InputDefault::is_joy_mapped(int p_device) {
	int mapping = joy_names[p_device].mapping;
	return mapping != -1 ? (mapping != fallback_mapping) : false;
}

String InputDefault::get_joy_guid_remapped(int p_device) const {
	ERR_FAIL_COND_V(!joy_names.has(p_device), "");
	return joy_names[p_device].uid;
}

Array InputDefault::get_connected_joypads() {
	Array ret;
	Map<int, Joypad>::Element *elem = joy_names.front();
	while (elem) {
		if (elem->get().connected) {
			ret.push_back(elem->key());
		}
		elem = elem->next();
	}
	return ret;
}

static const char* _buttons[] = {
	"Face Button Bottom",
	"Face Button Right",
	"Face Button Left",
	"Face Button Top",
	"L",
	"R",
	"L2",
	"R2",
	"L3",
	"R3",
	"Select",
	"Start",
	"DPAD Up",
	"DPAD Down",
	"DPAD Left",
	"DPAD Right"
};

static const char* _axes[] = {
	"Left Stick X",
	"Left Stick Y",
	"Right Stick X",
	"Right Stick Y",
	"",
	"",
	"L2",
	"R2"
};

String InputDefault::get_joy_button_string(int p_button) {
	ERR_FAIL_INDEX_V(p_button, JOY_BUTTON_MAX, "");
	return _buttons[p_button];
}

int InputDefault::get_joy_button_index_from_string(String p_button) {
	for (int i = 0; i < JOY_BUTTON_MAX; i++) {
		if (p_button == _buttons[i]) {
			return i;
		}
	}
	ERR_FAIL_V(-1);
}

String InputDefault::get_joy_axis_string(int p_axis) {
	ERR_FAIL_INDEX_V(p_axis, JOY_AXIS_MAX, "");
	return _axes[p_axis];
}

int InputDefault::get_joy_axis_index_from_string(String p_axis) {
	for (int i = 0; i < JOY_AXIS_MAX; i++) {
		if (p_axis == _axes[i]) {
			return i;
		}
	}
	ERR_FAIL_V(-1);
}
