/* Created by Language version: 5.5.1 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#define exp hoc_Exp
extern double hoc_Exp();
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define delta_t dt
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double dt;
 extern double t;
 /* declaration of user functions */
 static int _hoc_install_vector_fitness();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range("feature");
 _p = _prop->param; _ppvar = _prop->dparam;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_feature", _hoc_setdata,
 "install_vector_fitness", _hoc_install_vector_fitness,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static nrn_alloc(), nrn_init(), nrn_state();
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "5.5.1",
"feature",
 0,
 0,
 0,
 0};
 
static nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = (double *)ecalloc(0, sizeof(double));
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 0;
 
}
 static _initlists();
 _feature_reg_() {
	int _vectorized = 0;
  _initlists();
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 feature /Users/hines/neuron/nrndmg/src/nrnoc/../../../nrn/src/nrnoc/feature.mod\n");
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static install_vector_fitness();
 
/*VERBATIM*/
static double width(void* vv) {
	int i, nx, ny, i1, i2;
	double *x, *y;
	double h;
	h = *getarg(2);
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) return 0.;
	for (i1 = 0; i1 < nx; ++i1) {
		if (y[i1] >= h) { break;}
	}
	for (i2 = i1+1; i2 < nx; ++i2) {
		if (y[i2] <= h) { break; }
	}
	return x[i2] - x[i1];
}

/* xval, yval must be sorted in increasing order for xval */
/* xval is relative to the peak of the data */
/* xpeak is the peak location of the simulation */
/* xfitness is a measure of match in the x-dimension relative to
peak location for particular values of y */
/* yfitness is a measure of match in the y-dimension relative to
peak location for particular values of relative x */

static double xfitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	j = 0;
	sum = 0.;
	for (i = 0; i < nx; ++i) {
		if (y[i] >= yval[j]) {
			d = (x[i] - xpeak) - xval[j];
			sum += d*d;
			++j;
			if (j >= nxval) return sum;
			if (x[i] > xpeak) break;
		}
	}
	for (++i; i < nx; ++i) {
		if (y[i] <= yval[j]) {
			d = (x[i] - xpeak) - xval[j];
			sum += d*d;
			++j;
			if (j >= nxval) return sum;
		}
	}
	return 1e9;
}
static double yfitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	j = 0;
	sum = 0.;
	for (i = 0; i < nx; ++i) {
		if (x[i] - xpeak >= xval[j]) {
			d = y[i] - yval[j];
			sum += d*d;
			++j;
			if (j >= nxval) return sum;
		}
	}
	return 1e9;
}


static double firstpeak(void* vv) {
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = 0;
	while (i < ny) {
		if (y[i] >= -20) {
			if (y[i] > y[i+1]) {
				return (double) i;
			}
			i = i + 1;
		} else {
			i = i + 2;
		}
	}
	return 0.;
}
 
static int  install_vector_fitness (  )  {
   
/*VERBATIM*/
  {static int once; if (!once) { once = 1;
	install_vector_method("width", width);
	install_vector_method("xfitness", xfitness);
	install_vector_method("yfitness", yfitness);
	install_vector_method("firstpeak", firstpeak);
  }}
  return 0; }
 static int _hoc_install_vector_fitness() {
 double _r;
 _r = 1.;
 install_vector_fitness (  ) ;
 ret(_r);
}

static initmodel() {
  int _i; double _save;_ninits++;
{

}
}

static nrn_init(_nd, _pp, _ppd) Node *_nd; double *_pp; Datum* _ppd; {
 double _v;
 _p = _pp; _ppvar = _ppd;
 _v = _nd->_v;
 v = _v;
 initmodel();
}

static double _nrn_current(_v) double _v;{double _current=0.;v=_v;{
} return _current;
}

static nrn_state(_nd, _pp, _ppd) Node *_nd; double *_pp; Datum* _ppd; {
 double _break, _save;
 double _v;
 _p = _pp; _ppvar = _ppd;
 _v = _nd->_v;
 _break = t + .5*dt; _save = t; delta_t = dt;
 v=_v;
{

}
}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
