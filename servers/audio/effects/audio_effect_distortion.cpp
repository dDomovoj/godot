#include "audio_effect_distortion.h"
#include "servers/audio_server.h"



void AudioEffectDistortionInstance::process(const AudioFrame *p_src_frames,AudioFrame *p_dst_frames,int p_frame_count) {

	const float *src = (const float*)p_src_frames;
	float *dst = (float*)p_dst_frames;

	//float lpf_c=expf(-2.0*M_PI*keep_hf_hz.get()/(mix_rate*(float)OVERSAMPLE));
	float lpf_c=expf(-2.0*M_PI*base->keep_hf_hz/(AudioServer::get_singleton()->get_mix_rate()));
	float lpf_ic=1.0-lpf_c;

	float drive_f=base->drive;
	float pregain_f=Math::db2linear(base->pre_gain);
	float postgain_f=Math::db2linear(base->post_gain);

	float atan_mult=pow(10,drive_f*drive_f*3.0)-1.0+0.001;
	float atan_div=1.0/(atanf(atan_mult)*(1.0+drive_f*8));

	float lofi_mult=powf(2.0,2.0+(1.0-drive_f)*14); //goes from 16 to 2 bits

	for (int i=0;i<p_frame_count*2;i++) {

		float out=undenormalise(src[i]*lpf_ic+lpf_c*h[i&1]);
		h[i&1]=out;
		float a=out;
		float ha=src[i]-out; //high freqs
		a*=pregain_f;

		switch (base->mode) {

			case AudioEffectDistortion::MODE_CLIP: {

				a=powf(a,1.0001-drive_f);
				if (a>1.0)
					a=1.0;
				else if (a<(-1.0))
					a=-1.0;

			} break;
			case AudioEffectDistortion::MODE_ATAN: {


				a=atanf(a*atan_mult)*atan_div;

			} break;
			case AudioEffectDistortion::MODE_LOFI: {

				a = floorf(a*lofi_mult+0.5)/lofi_mult;

			} break;
			case AudioEffectDistortion::MODE_OVERDRIVE: {


				const double x = a * 0.686306;
				const double z = 1 + exp (sqrt (fabs (x)) * -0.75);
				a = (expf(x) - expf(-x * z)) / (expf(x) + expf(-x));
			} break;
			case AudioEffectDistortion::MODE_WAVESHAPE: {
				float x = a;
				float k= 2*drive_f/(1.00001-drive_f);

				a = (1.0+k)*x/(1.0+k*fabsf(x));


			} break;
		}

		dst[i]=a*postgain_f+ha;

	}


}


Ref<AudioEffectInstance> AudioEffectDistortion::instance() {
	Ref<AudioEffectDistortionInstance> ins;
	ins.instance();
	ins->base=Ref<AudioEffectDistortion>(this);
	ins->h[0]=0;
	ins->h[1]=0;

	return ins;
}


void AudioEffectDistortion::set_mode(Mode p_mode) {

	mode=p_mode;
}

AudioEffectDistortion::Mode AudioEffectDistortion::get_mode() const{

	return mode;
}

void AudioEffectDistortion::set_pre_gain(float p_pre_gain){

	pre_gain=p_pre_gain;
}
float AudioEffectDistortion::get_pre_gain() const{

	return pre_gain;
}

void AudioEffectDistortion::set_keep_hf_hz(float p_keep_hf_hz){

	keep_hf_hz=p_keep_hf_hz;
}
float AudioEffectDistortion::get_keep_hf_hz() const{

	return keep_hf_hz;
}

void AudioEffectDistortion::set_drive(float p_drive){

	drive=p_drive;
}
float AudioEffectDistortion::get_drive() const{

	return drive;
}

void AudioEffectDistortion::set_post_gain(float p_post_gain){

	post_gain=p_post_gain;
}
float AudioEffectDistortion::get_post_gain() const{

	return post_gain;
}


void AudioEffectDistortion::_bind_methods() {

	ClassDB::bind_method(_MD("set_mode","mode"),&AudioEffectDistortion::set_mode);
	ClassDB::bind_method(_MD("get_mode"),&AudioEffectDistortion::get_mode);

	ClassDB::bind_method(_MD("set_pre_gain","pre_gain"),&AudioEffectDistortion::set_pre_gain);
	ClassDB::bind_method(_MD("get_pre_gain"),&AudioEffectDistortion::get_pre_gain);

	ClassDB::bind_method(_MD("set_keep_hf_hz","keep_hf_hz"),&AudioEffectDistortion::set_keep_hf_hz);
	ClassDB::bind_method(_MD("get_keep_hf_hz"),&AudioEffectDistortion::get_keep_hf_hz);

	ClassDB::bind_method(_MD("set_drive","drive"),&AudioEffectDistortion::set_drive);
	ClassDB::bind_method(_MD("get_drive"),&AudioEffectDistortion::get_drive);


	ClassDB::bind_method(_MD("set_post_gain","post_gain"),&AudioEffectDistortion::set_post_gain);
	ClassDB::bind_method(_MD("get_post_gain"),&AudioEffectDistortion::get_post_gain);

	ADD_PROPERTY(PropertyInfo(Variant::INT,"mode",PROPERTY_HINT_ENUM,"Clip,ATan,LoFi,Overdrive,WaveShape"),_SCS("set_mode"),_SCS("get_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"pre_gain",PROPERTY_HINT_RANGE,"-60,60,0.01"),_SCS("set_pre_gain"),_SCS("get_pre_gain"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"keep_hf_hz",PROPERTY_HINT_RANGE,"1,20000,1"),_SCS("set_keep_hf_hz"),_SCS("get_keep_hf_hz"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"drive",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_drive"),_SCS("get_drive"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"post_gain",PROPERTY_HINT_RANGE,"-80,24,0.01"),_SCS("set_post_gain"),_SCS("get_post_gain"));
}

AudioEffectDistortion::AudioEffectDistortion()
{
	mode=MODE_CLIP;
	pre_gain=0;
	post_gain=0;
	keep_hf_hz=16000;
	drive=0;
}

