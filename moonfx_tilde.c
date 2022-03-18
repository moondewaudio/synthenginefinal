// moonfx~.cpp : Defines the exported functions for the DLL application.
//
#define PD_LONGINTTYPE long long
#include "m_pd.h"
//#include "stdafx.h"
#include <stdlib.h>
#include <math.h>

static t_class *moonfx_tilde_class;

typedef struct _moonfx_tilde {
	t_object x_obj;

	t_float oosamplerate, samplerate;
	t_float pi;

	// variables for folding
	t_float x_gain;
	t_float x_offset;

	// variables for BPF
	t_float x_cutoffBPF;
	double in0bp, in1bp, in2bp, out0bp, out1bp, out2bp;
	double gainbp, bwbp;

	// variables for delay
	t_float x_time;
	t_float x_fb;
	t_float x_mix;
	t_float peak;
	t_float *delayline;
	long delaysize;
	long wp;

	// variables for LPF
	t_float x_cutoffLPF;
	double in1, in2;

	// variables for phase shifter
	float pfreq, phase, pdepth;
	double pina, pinb, pinc, pind, pouta, poutb, poutc, poutd;

	t_float f;

	t_outlet *x_out;
} t_moonfx_tilde;


// peak detector for delay compression
float moonfx_tilde_peak_detect(t_moonfx_tilde *x, float in) {
	float rect;
	rect = in;
	if (rect < 0.0) rect *= -1.0;			// absolute value
	if (x->peak < rect)
		x->peak += (rect - x->peak) * 0.3f;	// smoothing on the peak value tracking to avoid large jumps and better show average level
	else									
		x->peak *= 0.999f;					// decay peak over time to track new ones
	return(x->peak);
}

void *moonfx_tilde_new(t_floatarg f) {
	t_moonfx_tilde *x = (t_moonfx_tilde *)pd_new(moonfx_tilde_class);
	outlet_new(&x->x_obj, gensym("signal"));

	// wavefolder inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_gain"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_offset"));

	// phasor inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("pfreq"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("pdepth"));

	// BPF inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_cutoffBPF"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("gainbp"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("bwbp"));

	// delay inlets
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_time"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_fb"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_mix"));

	// LPF inlet
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("x_cutoffLPF"));

	// defaults
	x->x_gain = 0.5f;
	x->x_offset = 0.0f;
	x->x_cutoffLPF = 20000.0f;
	x->x_cutoffBPF = 5000.0f;
	x->x_time = 1000.0f;
	x->x_fb = 0.0f;
	x->x_mix = 0.5f;
	x->pi = 4.0 * atan(1.0);
	x->peak = 0.1f;
	x->delaysize = 1048576;
	x->delayline = (t_float *)malloc(sizeof(t_float) * x->delaysize);
	x->wp = 0;
	x->in1 = x->in2 = x->in0bp = x->in1bp = x->in2bp = x->out0bp = x->out1bp = x->out2bp = 0.0;
	x->gainbp = 1;
	x->bwbp = 1;
	x->phase = 0.0f;
	x->pdepth = 1.0f;
	x->pfreq = 0.5f;
	x->pina = x->pinb = x->pinc = x->pind = x->pouta = x->poutb = x->poutc = x->poutd = 0.0;
	return (void *)x;
}

void moonfx_tilde_settime(t_moonfx_tilde *x, float f) {
	if (f > 15000) f = 15000.0f;	// maximum delay time in ms set to 15 seconds
	if (f < 1) f = 1.0f;			// delay time cant be below 1 ms
	x->x_time = f;
}

void moonfx_tilde_setmix(t_moonfx_tilde *x, float f) {
	if (f > 100.0) f = 100.0f;		// signal cannot be more than 100% wet or 0% wet
	if (f < 0.0) f = 0.0f;
	x->x_mix = f / 100.0f;			// scale from 0 to 1
}

void moonfx_tilde_setfb(t_moonfx_tilde *x, float f) {
	if (f > 120.0) f = 120.0f;		// maximum feedback above 100 to allow for self oscillation in the delay
	if (f < 0.0) f = 0.0f;			
	x->x_fb = f / 100.0f;			// scale from 0 to 1.2
}

void moonfx_tilde_setcutoffLPF(t_moonfx_tilde *x, float f) {
	if (f > 20000) f = 20000.0f;	// maximum cutoff for LPF is 20kHz
	if (f < 1) f = 1.0f;
	x->x_cutoffLPF = f;
}

void moonfx_tilde_setcutoffBPF(t_moonfx_tilde *x, float f) {
	if (f > 20000) f = 20000.0f;	// maximum cutoff for BPF is 20kHz
	if (f < 1) f = 1.0f;
	x->x_cutoffBPF = f;
}

void moonfx_tilde_setgain(t_moonfx_tilde *x, float f) {
	if (f < 0) f = 0.0f;			// signal cant have negative gain, but can scale infinitely
	x->x_gain = f / 10.0f;
}

void moonfx_tilde_setgainbp(t_moonfx_tilde *x, float f) {
	if (f < 0) f = 0.0f;			// BP gain must be between 0 and 1
	if (f > 100.0f) f = 100.0f;
	x->gainbp = f / 100.0f;			// scale knob for easier modification
}

void moonfx_tilde_setpdepth(t_moonfx_tilde *x, float f) {
	if (f < 0) f = 0.0f;			// phasor depth must be between 0 and 1
	if (f > 100) f = 100.0f;
	x->pdepth = f / 100.0f;			// scale knob for easier modification
}

void moonfx_tilde_setpfreq(t_moonfx_tilde *x, float f) {
	if (f > 100.0) f = 100.0f;		// max limit for phasor rate
	if (f < 0.0) f = 0.0f;			// rate cannot be below 0hz
	x->pfreq = f / 10.0f;			// scale the knob
}

void moonfx_tilde_setbwbp(t_moonfx_tilde *x, float f) {
	if (f < 0) f = 0.0f;			// BW must be above 0
	x->bwbp = f;			// scale knob for easier modification
}

void moonfx_tilde_setoffset(t_moonfx_tilde *x, float f) {
	x->x_offset = f / 100.0f;		// dc offset on signal can go as high or low as user wants for maximum waveshaping possibilities
}

// future plans: filter after clipping before wavefolder to give more intersting shape for folder to play with
// clip -> lpf -> folder -> BPF -> phase -> delay

static t_int *moonfx_tilde_perform(t_int *w)
{
	t_moonfx_tilde *x = (t_moonfx_tilde *)(w[1]);
	t_float *in = (t_float *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
	float cascade;
	float ingain;
	double tf, c; // tan freq, coefficient for LPF

	// vars for BP
	double d, tfbp, cbp;

	// vars for phasor
	float papfreq, ptf, pc, pinc, px0;

	// setup variables for the LPf
	tf = tan(x->pi * (x->x_cutoffLPF / x->samplerate));			// convert the cutoffLPF frequency from hz to radians
	c = (tf - 1.0) / (tf + 1.0);								// filter coefficient
	// and for the BPF
	d = -cos(2.0 * x->pi * (x->x_cutoffBPF / x->samplerate));	// filter coeff
	tfbp = tan(x->pi * (x->bwbp / x->samplerate));				// tan BW in radians
	cbp = (tfbp - 1.0) / (tfbp + 1.0);							// filter coeff

	// delay line
	float rpf, frac, x0, x1;
	long rp, rpp1;

	float fbinsum, input, rectif;

	float delaysamples = (x->x_time * 0.001 * x->samplerate);
	if (delaysamples < 1.0) delaysamples = 1.0;								// delay cant be less than 1 sample or longer than the delay line
	if (delaysamples >= x->delaysize) delaysamples = x->delaysize - 1.0;	

	int sample;
	for (sample = 0; sample < n; sample++) {
		// wavefolding based on a triangle waveshape, from http://synthnotes.ucsd.edu/wp4/index.php/2019/10/31/wavefolding/
		ingain = (x->x_gain * *(in + sample)) + x->x_offset;
		cascade = cos(0.5 * x->pi * ingain) 
			- 1.0 / 9.0 * cos(1.5 * x->pi * ingain)
			+ 1.0 / 25.0 * cos(2.5 * x->pi * ingain)
			- 1.0 / 49.0 * cos(3.5 * x->pi * ingain);

		// phaser post folding
		pinc = x->pfreq / x->samplerate;	// phase increment
		papfreq = (sin(x->phase * 2 * x->pi) + 1.0) * 2400.0f + 200.0f;	// all pass filter center freq, sweeps between 200 and 5k
		ptf = tan(x->pi * papfreq / x->samplerate);						// center freq in radians
		pc = (ptf - 1.0f) / (ptf + 1.0f);
		px0 = cascade;
		// cascade of APFs
		x->pouta = (pc * px0) + x->pina - pc * x->pouta;
		x->pina = px0;
		x->poutb = (pc * x->pouta) + x->pinb - pc * x->poutb;
		x->pinb = x->pouta;
		x->poutc = (pc * x->poutb) + x->pinc - pc * x->poutc;
		x->pinc = x->poutb;
		x->poutd = (pc * x->poutc) + x->pind - pc * x->poutc;
		x->pind = x->poutc;

		// sample output with depth
		cascade = cascade + x->pdepth * x->poutd;

		// increment the phase
		x->phase = x->phase + pinc;
		if (x->phase >= 1.0)
			x->phase = x->phase - 1.0;

		// bandpass the output post phase
		x->in0bp = cascade;
		x->out0bp = -cbp * x->in0bp + (d - d * cbp) * x->in1bp + x->in2bp - (d - d * cbp) * x->out1bp + cbp * x->out2bp;
		// move samples along
		x->in2bp = x->in1bp;
		x->in1bp = x->in0bp;
		x->out2bp = x->out1bp;
		x->out1bp = x->out0bp;
		cascade = (x->in0bp + x->out0bp * x->gainbp) * 0.5;

		// fractional delay line with mix control
		rpf = x->wp + delaysamples;
		rp = (long)rpf;
		frac = rpf - rp;

		rpp1 = rp + 1;
		if (rp >= x->delaysize) rp -= x->delaysize;		// wrap the read pointers to lie within the delay line
		if (rpp1 >= x->delaysize) rpp1 -= x->delaysize;
		x0 = *(x->delayline + rp);
		x1 = *(x->delayline + rpp1);
		
		input = cascade;
		fbinsum = x0 + (x1 - x0) * frac;				// output from the delay line, using linear interpolation to estimate values between samples
		
		// all pass based LPF to filter the repeats, from http://synthnotes.ucsd.edu/wp4/index.php/2019/10/14/first-order-low-pass-and-high-pass-filters/
		x->in2 = (c * fbinsum) + x->in1 - (c * x->in2);	// output
		x->in1 = fbinsum;								// current sample
		fbinsum = (x->in1 + x->in2)*0.5;
		*(out + sample) = input * (1 - x->x_mix) + fbinsum * x->x_mix;	// mix the input with the output. at 0, signal is fully dry. at 1, signal is fully wet.

		fbinsum = input + x->x_fb * fbinsum;							// new input to the delay line

		// compression
		moonfx_tilde_peak_detect(x, fbinsum);
		if (x->peak > 2.0) fbinsum *= 0.5;
		if (x->peak > 0.5)
			fbinsum *= (1.601539 + x->peak * (-1.605725 + x->peak * (0.8883899 - 0.180484 * x->peak))); // polynomial compression algorithm, from http://synthnotes.ucsd.edu/wp4/index.php/2019/11/16/delay-with-feedback/
		if (fbinsum > 1.0) fbinsum = 1.0; // set feedback gain to be within -1 and 1 to avoid any clipping in the delay
		if (fbinsum < -1.0) fbinsum = -1.0;
		*(x->delayline + x->wp) = fbinsum; // write to the delay line
		(x->wp)--;
		if (x->wp < 0) x->wp += x->delaysize; // wrap the write pointer back into the delay line if necessary
	}
	return(w + 5);
}

void moonfx_tilde_free(t_moonfx_tilde *x) {
	free(x->delayline);
}

void moonfx_tilde_dsp(t_moonfx_tilde *x, t_signal **sp) {
	x->oosamplerate = 1.0 / sp[0]->s_sr;
	x->samplerate = sp[0]->s_sr;
	dsp_add(moonfx_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void moonfx_tilde_setup(void) {
	moonfx_tilde_class = class_new(gensym("moonfx~"),
		(t_newmethod)moonfx_tilde_new,
		(t_method)moonfx_tilde_free,
		sizeof(t_moonfx_tilde),
		CLASS_DEFAULT,
		A_DEFFLOAT, 0);

	CLASS_MAINSIGNALIN(moonfx_tilde_class, t_moonfx_tilde, f);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_dsp, gensym("dsp"), A_CANT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setmix,
		gensym("x_mix"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_settime,
		gensym("x_time"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setfb,
		gensym("x_fb"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setcutoffLPF,
		gensym("x_cutoffLPF"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setgain,
		gensym("x_gain"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setoffset,
		gensym("x_offset"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setcutoffLPF,
		gensym("x_cutoffLPF"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setgainbp,
		gensym("gainbp"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setbwbp,
		gensym("bwbp"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setcutoffBPF,
		gensym("x_cutoffBPF"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setpdepth,
		gensym("pdepth"), A_DEFFLOAT, 0);

	class_addmethod(moonfx_tilde_class,
		(t_method)moonfx_tilde_setpfreq,
		gensym("pfreq"), A_DEFFLOAT, 0);
}


