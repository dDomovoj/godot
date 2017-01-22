#include "audio_effect_stereo_enhance.h"
#include "servers/audio_server.h"
void AudioEffectStereoEnhanceInstance::process(const AudioFrame *p_src_frames,AudioFrame *p_dst_frames,int p_frame_count) {


	float intensity=base->pan_pullout;
	bool surround_mode=base->surround>0;
	float surround_amount=base->surround;
	unsigned int delay_frames=(base->time_pullout/1000.0)*AudioServer::get_singleton()->get_mix_rate();

	for (int i=0;i<p_frame_count;i++) {

		float l=p_src_frames[i].l;
		float r=p_src_frames[i].r;

		float center=(l+r)/2.0f;

		l=( center+(l-center)*intensity );
		r=( center+(r-center)*intensity );

		if (surround_mode) {

			float val=(l+r)/2.0;

			delay_ringbuff[ringbuff_pos&ringbuff_mask]=val;

			float out=delay_ringbuff[(ringbuff_pos-delay_frames)&ringbuff_mask]*surround_amount;

			l+=out;
			r+=-out;
		} else {

			float val=r;

			delay_ringbuff[ringbuff_pos&ringbuff_mask]=val;

			//r is delayed
			r=delay_ringbuff[(ringbuff_pos-delay_frames)&ringbuff_mask];;


		}

		p_dst_frames[i].l=l;
		p_dst_frames[i].r=r;
		ringbuff_pos++;

	}

}


AudioEffectStereoEnhanceInstance::~AudioEffectStereoEnhanceInstance() {

	memdelete_arr(delay_ringbuff);
}

Ref<AudioEffectInstance> AudioEffectStereoEnhance::instance() {
	Ref<AudioEffectStereoEnhanceInstance> ins;
	ins.instance();

	ins->base=Ref<AudioEffectStereoEnhance>(this);


	float ring_buffer_max_size=AudioEffectStereoEnhanceInstance::MAX_DELAY_MS+2;
	ring_buffer_max_size/=1000.0;//convert to seconds
	ring_buffer_max_size*=AudioServer::get_singleton()->get_mix_rate();

	int ringbuff_size=(int)ring_buffer_max_size;

	int bits=0;

	while(ringbuff_size>0) {
		bits++;
		ringbuff_size/=2;
	}

	ringbuff_size=1<<bits;
	ins->ringbuff_mask=ringbuff_size-1;
	ins->ringbuff_pos=0;

	ins->delay_ringbuff = memnew_arr(float,ringbuff_size );

	return ins;
}

void  AudioEffectStereoEnhance::set_pan_pullout(float p_amount) {

	pan_pullout=p_amount;
}

float  AudioEffectStereoEnhance::get_pan_pullout() const {

	return pan_pullout;
}

void  AudioEffectStereoEnhance::set_time_pullout(float p_amount) {

	time_pullout=p_amount;
}
float  AudioEffectStereoEnhance::get_time_pullout() const {

	return time_pullout;
}

void  AudioEffectStereoEnhance::set_surround(float p_amount) {

	surround=p_amount;
}
float  AudioEffectStereoEnhance::get_surround() const {

	return surround;
}

void AudioEffectStereoEnhance::_bind_methods() {

	ClassDB::bind_method(_MD("set_pan_pullout","amount"),&AudioEffectStereoEnhance::set_pan_pullout);
	ClassDB::bind_method(_MD("get_pan_pullout"),&AudioEffectStereoEnhance::get_pan_pullout);

	ClassDB::bind_method(_MD("set_time_pullout","amount"),&AudioEffectStereoEnhance::set_time_pullout);
	ClassDB::bind_method(_MD("get_time_pullout"),&AudioEffectStereoEnhance::get_time_pullout);

	ClassDB::bind_method(_MD("set_surround","amount"),&AudioEffectStereoEnhance::set_surround);
	ClassDB::bind_method(_MD("get_surround"),&AudioEffectStereoEnhance::get_surround);

	ADD_PROPERTY(PropertyInfo(Variant::REAL,"pan_pullout",PROPERTY_HINT_RANGE,"0,4,0.01"),_SCS("set_pan_pullout"),_SCS("get_pan_pullout"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"time_pullout_ms",PROPERTY_HINT_RANGE,"0,50,0.01"),_SCS("set_time_pullout"),_SCS("get_time_pullout"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"surround",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_surround"),_SCS("get_surround"));
}

AudioEffectStereoEnhance::AudioEffectStereoEnhance()
{
	pan_pullout=1;
	time_pullout=0;
	surround=0;
}
