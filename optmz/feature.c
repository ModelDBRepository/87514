/*********************   feature.c   in src/nrnoc *****************************/

/* Created by Language version: 5.5.1 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#define CHRIS_DBG 0
#define DBG_NEGSHP 0
#define EPS 0.000001

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
   static double lin_interp();
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
	//FILE *chrisfp;
	//fprintf(stdout,"XFITNESS IN C FILE\n");
	//chrisfp=fopen("chris_xfitness.txt","w");
	//fprintf(chrisfp,"XFITNESS IN C FILE\n");
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	//fprintf(chrisfp,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
	j = 0;
	sum = 0.;
	for (i = 0; i < nx; ++i) {
		if (y[i] >= yval[j]) {
			d = (x[i] - xpeak) - xval[j];
			sum += d*d;
			//fprintf(chrisfp,"[1]\tx,y=(%g, %g)\txval,yval=(%g, %g); xpeak %g; d = %g; sum %g\n",x[i],y[i],xval[j],yval[j],xpeak,d,sum);
			++j;
			if (j >= nxval) 
			  //{ fprintf(chrisfp,"final sum %g\n",sum); fclose(chrisfp);
return sum;
			//}
			if (x[i] > xpeak) break;
		}
	}
	for (++i; i < nx; ++i) {
		if (y[i] <= yval[j]) {
			d = (x[i] - xpeak) - xval[j];
			sum += d*d;
			//fprintf(chrisfp,"[2]\tx,y=(%g, %g)\txval,yval=(%g, %g); xpeak %g; d = %g; sum %g\n",x[i],y[i],xval[j],yval[j],xpeak,d,sum);
			++j;
			if (j >= nxval) //{ fprintf(chrisfp,"final sum %g\n",sum); fclose(chrisfp);
return sum;
			//}
		}
	}
	//fclose(chrisfp);
	fprintf(stdout,"X FITNESS: End of model X values; returning sum = 1e9\n");
	return 1e9;
}

/***********************************************************

Same as yfitness below, except that this function reports an
error proportional to the size of the experimental window.  
This is most useful in APs near the start or the end of a 
model run, where data in the width of the entire experimental 
window may not be available.  

    This function is called by the Multiple Run Fitter error functions that use AP shape.  The function is called by, e.g.:
	    eval = $o1.ywnscl_fitness($o2, peak, ytmp, xtmp,ytmp.firstpeak,tmp_modind,mnind,mxind)

    so that the variables defined below correspond to the following:

    	y		the y-values of model output
        x		the t-values of model output
        xpeak		time of the specified model AP peak
	yval		the y-values of the target output
	xval		the t-values of the target output
	val_pkind	index of the AP peak of the target data
        pkind		index of peak of the model AP whose shape error will be calculated
        mnind		first index of model data within time bounds
        mxind		last  index of model data within time bounds

***********************************************************/
static double ywnscl_fitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind, mod_start, mod_end, mnind, mxind, n_pts;
	double ytmp, val_winsz, mod_winsz;

#if CHRIS_DBG
FILE *chrisfp;
// when debugging several instances at once, this will generate a HUGE file!
chrisfp=fopen("ywnscl_fit.m","w");
//chrisfp=fopen("ywnscl_fit.m","a");
fprintf(chrisfp,"\n\n%%YWINFIT IN C FILE\n");
fprintf(chrisfp,"%%[i,modelx-modelpeak,exp_peak] modelx,modely=(,) expx,expy=(,)\n");
#endif

#if DBG_NEGSHP
FILE *chrisfp;
chrisfp=fopen("dbg_negshp.out","a");
#endif

	//
	// get parameters from NEURON function
	ny = vector_instance_px(vv, &y);		// y is the vector
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
	mnind = *getarg(7);
	mxind = *getarg(8);

#if CHRIS_DBG
fprintf(chrisfp,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
 fprintf(chrisfp,"model pkind = %d [%g, %g], exp pk ind %d [%g, %g]\n",pkind,x[pkind],y[pkind],val_pkind,xval[val_pkind],yval[val_pkind]);
//fprintf(stdout,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
//fprintf(stdout,"model pkind = %d [%g, %g], exp pk ind %d [%g, %g]; time indices [%g, %g]\n",pkind,x[pkind],y[pkind],val_pkind,xval[val_pkind],yval[val_pkind],x[mnind],x[mxind]);
#endif

	j = 0;
	sum = 0.;
	n_pts = 0;

	//
	// Calculate the error directly at, before, and after the
	// peak, for each of the experimental data points included in (xval,yval).
	// Interpolate when data points do not occur at the same x-values.

#if CHRIS_DBG
fprintf(chrisfp,"[%d, %d] Peak:  M (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",pkind,val_pkind,x[pkind],xpeak,x[pkind]-xpeak,y[pkind],xval[val_pkind],yval[val_pkind],d,sum);
#endif

        // start at the beginning of the model window, and the beginning 
        // of the experimental window. 
 mod_start = mxind+1;
 mod_end = mnind-1;
        i = mnind;
        j = 0; 

#if DBG_NEGSHP
	if( i >= mxind || j >= nyval ) { 
fprintf(chrisfp,"\tAFTER loop never entered:  pkind=%d, mxind=%d, val_pkind+1=%d, nyval=%d.\n",i,mxind,j,nyval);
	} 
	if( i <= 0 ) { 
	  fprintf(chrisfp,"\ti = pkind = %d, j = val_pkind+1 %d vs. nyval %d\n",pkind,j,nyval);
	}
#endif

	while( i < mxind && j < nyval ) {

#if CHRIS_DBG
fprintf(chrisfp,"\tLooking to match [%d] E (%g, %g)\n",j,xval[j],yval[j]);
#endif
	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is greater than the current experiment value
	    if( xval[j] > (x[i] - xpeak) ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak && (fabs(xval[j] - (x[i] - xpeak)) > EPS) ) {
		  #if CHRIS_DBG
		  fprintf(chrisfp,"\t[%d, %d] compare %g > %g = %g Skipping M (%g - %g = %g, %g) for E (%g, %g)\n",\
                          i,j,xval[j],x[i]-xpeak,xval[j]-(x[i]-xpeak),x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j]);
		  #endif
	            i++;
		}
	    }


	    #if CHRIS_DBG
	    fprintf(chrisfp,"\tFound [%d] (%g - %g = %g, %g) to match [%d] E (%g, %g)\n",\
                    i,x[i],xpeak,x[i]-xpeak,y[i],j,xval[j],yval[j]);
	    #endif

	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
#if CHRIS_DBG
fprintf(chrisfp,"I\t[%d, %d] X equal (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",\
        i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],d,sum);
#endif
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		// end of both experiment & model; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M and E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",\
        i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",\
        ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
fprintf(stdout,"II.  \t[%d, %d] end M and E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",\
        i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(stdout,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",\
        ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
#endif
	    } else if((i==ny-1) && (j < nyval) ) {
		// j is not at end of list; interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M not E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
fprintf(stdout,"III.\t[%d, %d] end M not E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(stdout,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
#endif
	    } else {
	        //interpolate model backward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] middle M and E; IntMod Bck  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
fprintf(stdout,"IV.\t[%d, %d] middle M and E; IntMod Bck  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
fprintf(stdout,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
#endif
	    }

	    j++;
	}

	// end of model AP window
	mod_end = i;
	if( mod_end == nx ) { 
#if DBG_NEGSHP
	fprintf(chrisfp,"\tdecrementing mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
	fprintf(stdout,"\tdecrementing mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
#endif
#if CHRIS_DBG
	fprintf(chrisfp,"\tdecrementing mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
	//fprintf(stdout,"\tdecrementing mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
#endif
	     mod_end--; 
	}

#if CHRIS_DBG
    if( n_pts != nyval ) {
	fprintf(stdout,"SHAPE:  End of model window i = %d/%d; j = %d/%d.  Found mod_start %d, mod_end %d, Sum = %g from %d points\n",i,mxind,j,nyval,mod_start,mod_end,sum,n_pts);
fprintf(stdout,"\tmodel pkind = %d [%g, %g], exp pk ind %d [%g, %g]; time indices [%g, %g]\n",pkind,x[pkind],y[pkind],val_pkind,xval[val_pkind],yval[val_pkind],x[mnind],x[mxind]);
    }
#endif

#if DBG_NEGSHP
	fprintf(chrisfp,"\tSET mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
	if( mod_end <= 0 ) { 
	  fprintf(chrisfp,"\tmod_end = %d\tj = %d\n",mod_end,j);
	} 
#endif
#if CHRIS_DBG
	fprintf(chrisfp,"\tSET mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
	if( mod_end <= 0 ) { 
	  fprintf(chrisfp,"\tmod_end = %d\tj = %d\n",mod_end,j);
	} 
	//fprintf(stdout,"\tSET mod_end = %d (%g), max size %d\n",mod_end,x[mod_end],nx);
	if( mod_end <= 0 ) { 
	  fprintf(stdout,"\tmod_end = %d\tj = %d\n",mod_end,j);
	} 
        //fprintf(stdout,"Error sum = %g, from %d points.\n",sum,n_pts);
#endif
	
	// root mean squared error
        if( n_pts == 0 ) { 
	    sum = -1;
        } else {
  	    sum = sqrt(sum/(double)n_pts);
	}

/****
#if CHRIS_DBG
	fprintf(stdout,"\n\n\nBefore window scale, sum = %g\n",sum);
#endif
****/

#if DBG_NEGSHP
	fprintf(stdout,"\n\n\nBefore window scale, sum = %g\n",sum);
#endif

        // total size (in ms) of the experimental & model windows;
	// scale the calculated error proportional to the experimental error size
	val_winsz = xval[nxval-1] - xval[0];
        if( n_pts == 0 ) { 
	    mod_winsz = 0;
	} else {
	    mod_winsz = x[mod_end] - x[mod_start];
	}
#if CHRIS_DBG
	fprintf(chrisfp,"\n\nsize_full = %g - %g = %g;\n",xval[nxval-1], xval[0], val_winsz);
	fprintf(chrisfp,"size_trunc = x[%d] - x[%d] = %g;\n",mod_end, mod_start, mod_winsz);
	fprintf(chrisfp,"before scale, sum = %g\nafter scale = %g;\n",sum,mod_winsz*sum/val_winsz);
#endif
#if DBG_NEGSHP
	fprintf(chrisfp,"\nDONE mod_start=%d (x %g), mod_end=%d (x %g)\n",mod_start,x[mod_start],mod_end,x[mod_end]);
	fprintf(chrisfp,"\nsize_full = %g - %g = %g;\t",xval[nxval-1], xval[0], val_winsz);
	fprintf(chrisfp,"size_trunc = %g - %g = %g; [n_pts = %d]\n",x[mod_end], x[mod_start], mod_winsz,n_pts);
	fprintf(chrisfp,"before scale, sum = %g\tafter scale = %g;\n",sum,mod_winsz*sum/val_winsz);
#endif
        if( mod_winsz > 0 ) {
	    sum = mod_winsz*sum/val_winsz; 
	} else {
	    sum = -1.0;
	}
/****
#if CHRIS_DBG
	fprintf(stdout,"now after scale = %g;\n",sum);
#endif
****/

#if CHRIS_DBG
fclose(chrisfp);
	//fprintf(stdout,"exp window = x[%d] %g to x[%d] %g \t size %g\n",0,xval[0],nxval-1,xval[nxval-1],val_winsz);
	//fprintf(stdout,"mod window = x[%d] %g to x[%d] %g \t size %g\n",mod_start,x[mod_start],mod_end,x[mod_end],mod_winsz);
	//fprintf(stdout,"returning sum = %g\n",sum);
#endif
#if DBG_NEGSHP
fprintf(stdout,"now after scale = %g;\n",sum);
fclose(chrisfp);
#endif
	return sum;

}


/***********************************************************

Same as yfitness below, except that this function reports an
error proportional to the size of the experimental window.  
This is most useful in APs near the start or the end of a 
model run, where data in the width of the entire experimental 
window may not be available.  

***********************************************************/
static double ywinfitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind, mod_start, mod_end;
	double ytmp, val_winsz, mod_winsz;
#if CHRIS_DBG
FILE *chrisfp;
chrisfp=fopen("chris_ywinfit.txt","a");
fprintf(chrisfp,"\n\nYWINFIT IN C FILE\n");
fprintf(chrisfp,"%%[i,modelx-modelpeak,exp_peak] modelx,modely=(,) expx,expy=(,)\n");
#endif
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
#if CHRIS_DBG
fprintf(chrisfp,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
 fprintf(chrisfp,"model pkind = %d [%g, %g], exp pk ind %d [%g, %g]\n",pkind,x[pkind],y[pkind],val_pkind,xval[val_pkind],yval[val_pkind]);
#endif

	j = 0;
	sum = 0.;

	//
	// Calculate the error directly at, before, and after the
	// peak, for each of the experimental data points included in (xval,yval).
	// Interpolate when data points do not occur at the same x-values.

	// first, error at the peak
	d = y[pkind] - yval[val_pkind];
	sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"[%d, %d] Peak:  M (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",pkind,val_pkind,x[pkind],xpeak,x[pkind]-xpeak,y[pkind],xval[val_pkind],yval[val_pkind],d,sum);
fprintf(chrisfp,"\n\nBEFORE PEAK:\n");
#endif

	// now, error before the peak
        i = pkind; 
	j = val_pkind -1;
	while( i > 0 && j >= 0 ) {

#if CHRIS_DBG
fprintf(chrisfp,"\tLooking to match [%d] E (%g, %g)\n",j,xval[j],yval[j]);
#endif
	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is less than the current experiment value
	    if( xval[j] < x[i] - xpeak ) {
	        while( (i>0) && xval[j] < x[i] - xpeak ) {
		  //#if CHRIS_DBG
		  //fprintf(chrisfp,"\t[%d, %d] Skipping M (%g - %g = %g, %g) for E (%g, %g)\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j]);
		  //#endif

	            i--;
		}
	    }

	    //#if CHRIS_DBG
	    //fprintf(chrisfp,"\tFound [%d] (%g - %g = %g, %g) to match [%d] E (%g, %g)\n",i,x[i],xpeak,x[i]-xpeak,y[i],j,xval[j],yval[j]);
	    //#endif

	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] X equal (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],d,sum);
#endif
	    } else if( j==0 && i==0 && (xval[j] < x[i]-xpeak) ) {
		// interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] start M and E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
#endif
	    } else if( i==0 && j>0 ) {
		// j is nonzero; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] start M not E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
#endif
	    } else {
	        // interpolate model forward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
	        d = ytmp - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] middle M and E; IntMod Fwd  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
#endif
	    }

	    j--;
	}

	// start of model AP window
	mod_start = i;

#if CHRIS_DBG
fprintf(chrisfp,"\n\nAFTER PEAK:\n");
#endif
	// error after the peak
        i = pkind; // + 1;
	j = val_pkind + 1;
	while( i < ny-1 && j < nyval ) {

#if CHRIS_DBG
fprintf(chrisfp,"\tLooking to match [%d] E (%g, %g)\n",j,xval[j],yval[j]);
#endif
	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is greater than the current experiment value
	    if( xval[j] > x[i] - xpeak ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak ) {
		  //#if CHRIS_DBG
		  //fprintf(chrisfp,"\t[%d, %d] Skipping M (%g - %g = %g, %g) for E (%g, %g)\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j]);
		  //#endif
	            i++;
		}
	    }


	    //#if CHRIS_DBG
	    //fprintf(chrisfp,"\tFound [%d] (%g - %g = %g, %g) to match [%d] E (%g, %g)\n",i,x[i],xpeak,x[i]-xpeak,y[i],j,xval[j],yval[j]);
	    //#endif

	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] X equal (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],d,sum);
#endif
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		// end of both experiment & model; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M and E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
#endif
	    } else if((i==ny-1) && (j < nyval) ) {
		// j is not at end of list; interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M not E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
#endif
	    } else {
	        //interpolate model backward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] middle M and E; IntMod Bck  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
#endif
	    }

	    j++;
	}

	// end of model AP window
	mod_end = i;

#if CHRIS_DBG
fclose(chrisfp);
	fprintf(stdout,"\n\n\nBefore window scale, sum = %g\n",sum);
#endif

        // total size (in ms) of the experimental & model windows;
	// scale the calculated error proportional to the experimental error size
	val_winsz = xval[nxval-1] - xval[0];
	mod_winsz = x[mod_end] - x[mod_start];
	sum = (val_winsz / mod_winsz) * sum;

#if CHRIS_DBG
fclose(chrisfp);
	fprintf(stdout,"Y FITNESS: End of model X values; \n");
	fprintf(stdout,"exp window = x[%d] %g to x[%d] %g \t size %g\n",0,xval[0],nxval-1,xval[nxval-1],val_winsz);
	fprintf(stdout,"mod window = x[%d] %g to x[%d] %g \t size %g\n",mod_start,x[mod_start],mod_end,x[mod_end],mod_winsz);
	fprintf(stdout,"returning sum = %g\n",sum);
#endif
	return sum;

}




// Modified, Aug 2004 by Christina Weaver
static double yfitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind;
	double ytmp;
#if CHRIS_DBG
FILE *chrisfp;
chrisfp=fopen("chris_yfitness.txt","a");
fprintf(chrisfp,"\n\nYFITNESS IN C FILE\n");
fprintf(chrisfp,"%%[i,modelx-modelpeak,exp_peak] modelx,modely=(,) expx,expy=(,)\n");
#endif
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
#if CHRIS_DBG
fprintf(chrisfp,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
 fprintf(chrisfp,"model pkind = %d [%g, %g], exp pk ind %d [%g, %g]\n",pkind,x[pkind],y[pkind],val_pkind,xval[val_pkind],yval[val_pkind]);
#endif
	j = 0;
	sum = 0.;

	// DON'T PENALIZE if the model AP window is not complete at start of the run,
	// or at end of the run.
	//
	// Calculate the error directly at, before, and after the
	// peak, for each of the experimental data points included in (xval,yval).
	// Interpolate when data points do not occur at the same x-values.

	// first, error at the peak
	d = y[pkind] - yval[val_pkind];
	sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"[%d, %d] Peak:  M (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",pkind,val_pkind,x[pkind],xpeak,x[pkind]-xpeak,y[pkind],xval[val_pkind],yval[val_pkind],d,sum);
fprintf(chrisfp,"\n\nBEFORE PEAK:\n");
#endif

	// now, error before the peak
        i = pkind; // -1;
	j = val_pkind -1;
	while( i > 0 && j >= 0 ) {

#if CHRIS_DBG
fprintf(chrisfp,"\tLooking to match [%d] E (%g, %g)\n",j,xval[j],yval[j]);
#endif
	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is less than the current experiment value
	    if( xval[j] < x[i] - xpeak ) {
	        while( (i>0) && xval[j] < x[i] - xpeak ) {
		  //#if CHRIS_DBG
		  //fprintf(chrisfp,"\t[%d, %d] Skipping M (%g - %g = %g, %g) for E (%g, %g)\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j]);
		  //#endif

	            i--;
		}
	    }

	    //#if CHRIS_DBG
	    //fprintf(chrisfp,"\tFound [%d] (%g - %g = %g, %g) to match [%d] E (%g, %g)\n",i,x[i],xpeak,x[i]-xpeak,y[i],j,xval[j],yval[j]);
	    //#endif

	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] X equal (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],d,sum);
#endif
	    } else if( j==0 && i==0 && (xval[j] < x[i]-xpeak) ) {
		// interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] start M and E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
#endif
	    } else if( i==0 && j>0 ) {
		// j is nonzero; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] start M not E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
#endif
	    } else {
	        // interpolate model forward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
	        d = ytmp - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] middle M and E; IntMod Fwd  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
#endif
	    }

	    j--;
	}

#if CHRIS_DBG
fprintf(chrisfp,"\n\nAFTER PEAK:\n");
#endif
	// error after the peak
        i = pkind; // + 1;
	j = val_pkind + 1;
	while( i < ny-1 && j < nyval ) {

#if CHRIS_DBG
fprintf(chrisfp,"\tLooking to match [%d] E (%g, %g)\n",j,xval[j],yval[j]);
#endif
	    // model & experiment may have different (even variable) values of dt;
	    // find the appropriate x value, and interpolate accordingly.

	    // find next model value which is greater than the current experiment value
	    if( xval[j] > x[i] - xpeak ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak ) {
		  //#if CHRIS_DBG
		  //fprintf(chrisfp,"\t[%d, %d] Skipping M (%g - %g = %g, %g) for E (%g, %g)\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j]);
		  //#endif
	            i++;
		}
	    }


	    //#if CHRIS_DBG
	    //fprintf(chrisfp,"\tFound [%d] (%g - %g = %g, %g) to match [%d] E (%g, %g)\n",i,x[i],xpeak,x[i]-xpeak,y[i],j,xval[j],yval[j]);
	    //#endif

	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        // both model & experiment were evaluated at the same t value
	        d = y[i] - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] X equal (%g - %g = %g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],d,sum);
#endif
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		// end of both experiment & model; interpolate experiment backward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M and E; IntExp Bck  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
#endif
	    } else if((i==ny-1) && (j < nyval) ) {
		// j is not at end of list; interpolate experiment forward
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] end M not E; IntExp Fwd  (%g - %g = %g, %g), E (%g, %g), I (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],yval[j],x[i],ytmp,d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
#endif
	    } else {
	        //interpolate model backward
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"\t[%d, %d] middle M and E; IntMod Bck  (%g - %g = %g, %g), I (%g, %g), E (%g, %g)\td = %g, err = %g\n",i,j,x[i],xpeak,x[i]-xpeak,y[i],xval[j],ytmp,xval[j],yval[j],d,sum);
//fprintf(chrisfp,"\tInt %g = y for x = %g, from (%g, %g) to (%g, %g)\n",ytmp,xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
#endif
	    }

	    j++;
	}

#if CHRIS_DBG
fclose(chrisfp);
	fprintf(stdout,"Y FITNESS: End of model X values; returning sum = %g\n",sum);
#endif
	return sum;

	/****
	// this assumes that there are sufficient points at the beginning of the model 
	// trace to calculate the errors with respect to the experimental AP.
	// HOW SHOULD WE HANDLE IT OTHERWISE?
	for (i = 0; i < nx; ++i) {
	  if (x[i] - xpeak >= xval[j]) {
			d = y[i] - yval[j];
			sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"[%d, %g, %g]\tx,y=(%g, %g)\txval,yval=(%g, %g); xpeak %g; d = %g; sum %g\n",i,x[i]-xpeak,xval[j],x[i],y[i],xval[j],yval[j],xpeak,d,sum);
#endif
			++j;
			if (j >= nxval) { 
#if CHRIS_DBG
fprintf(chrisfp,"final sum %g\n",sum); fclose(chrisfp);
#endif
return sum;
}
		}
	}
#if CHRIS_DBG
fclose(chrisfp);
	fprintf(stdout,"Y FITNESS: End of model X values; returning sum = 1e9\n");
#endif
	return 1e9;
	****/
}


// Original version of yfitness coded by Michael Hines
static double yfitness_hines(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
#if CHRIS_DBG
FILE *chrisfp;
fprintf(stdout,"\n\nYFITNESS IN C FILE\n");
chrisfp=fopen("chris_yfitness.txt","a");
fprintf(chrisfp,"\n\nYFITNESS IN C FILE\n");
fprintf(chrisfp,"%%[i,modelx-modelpeak,exp_peak] modelx,modely=(,) expx,expy=(,)\n");
#endif
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
#if CHRIS_DBG
fprintf(chrisfp,"nx %d ny %d nyval %d nxval %d\n",nx,ny,nyval,nxval);
#endif
	j = 0;
	sum = 0.;

	// this assumes that there are sufficient points at the beginning of the model 
	// trace to calculate the errors with respect to the experimental AP.
	// HOW SHOULD WE HANDLE IT OTHERWISE?
	for (i = 0; i < nx; ++i) {
	  if (x[i] - xpeak >= xval[j]) {
			d = y[i] - yval[j];
			sum += d*d;
#if CHRIS_DBG
fprintf(chrisfp,"[%d, %g, %g]\tx,y=(%g, %g)\txval,yval=(%g, %g); xpeak %g; d = %g; sum %g\n",i,x[i]-xpeak,xval[j],x[i],y[i],xval[j],yval[j],xpeak,d,sum);
#endif
			++j;
			if (j >= nxval) { 
#if CHRIS_DBG
fprintf(chrisfp,"final sum %g\n",sum); fclose(chrisfp);
#endif
return sum;
}
		}
	}
#if CHRIS_DBG
fclose(chrisfp);
	fprintf(stdout,"Y FITNESS: End of model X values; returning sum = 1e9\n");
#endif
	return 1e9;
}


static double lin_interp(xstar, x1, y1, x2, y2) double xstar, x1, y1, x2, y2; {
        double ystar;

	if( fabs(x2-x1) < EPS ) return 0.5*(y1+y2);

	ystar = y1 + ((y2-y1)/(x2-x1))*(xstar-x1);
	return ystar;
}

/************************************

cmw NOTE:  the firstpeak() function only recognizes 
peaks with voltage >= -20 mV.  Otherwise 0 is returned.

************************************/
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

/************************************

    firstmax

    this function finds the first local max of the vector y.

************************************/
static double firstmax(void* vv) { 
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = 0;
	while (i < ny) {
		if (y[i] > y[i+1]) {
			return (double) i;
		}
		i = i + 1;
	}
	return 0.;
}


static double nextpeak(void* vv) {
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = *getarg(1);
	//i = *getarg(2);
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
 

static	float	sqrarg;

#define	SQR(a) (sqrarg=(a),sqrarg*sqrarg)

/*
*       From NUMERICAL RECIPES IN C.
*   
*	Given a set of points x[strt ... end], y[strt ... end], with standard
*	deviations sig[strt ... end], fit them to
*
*			y = a + bx		(!!!!!!!!!)
*
*	by minimizing chi_sq.
*	Returned are a,b, siga, sigb, chi_sq (chi2), and the goodness of fit
*	probability q (that the fit would have chi_sq this large or larger).
*	If mwt = 0 *	on input, then the standard deviations are assumed
*	unavailable, q is returned as 1.0, and the normalization of chi2
*	is to unit standard deviation on all points.
*/

void	linfit(void* vv) {

  /** mwt = 0, no 'sig' array is given , no q returned.  
      references to these values have been deleted. **/
	float	*x, *y;
	int	strt, end;
	float	*a, *b, *siga, *sigb, *chi2;

	int	i,nx,ny;
	float	wt, t, sxoss, ss, sigdat;
	float	sx  = 0.0;
	float	sy  = 0.0;
	float	st2 = 0.0;

	ny = vector_instance_px(vv,&y);
	nx = vector_instance_px(1,&x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	strt = *getarg(2);
	end  = *getarg(3);
	*a    = *getarg(4);
	*b    = *getarg(5);
	*siga = *getarg(6);
	*sigb = *getarg(7);
	*chi2 = *getarg(8);

#if CHRIS_DBG
	fprintf(stdout,"x,y size = %d; start = %g end = %d\n",ny,strt,end);
	fprintf(stdout,"a=%g; b=%g; siga=%g; sigb=%g; chi2=%g\n",*a,*b,*siga,*sigb,*chi2);
#endif

	*b = 0.0;
	for( i = strt;  i <= end;  i++ )
	{
		sx += x[i];	sy += y[i];
	}
	ss = end+1-strt;

	sxoss = sx/ss;

	for( i = strt;  i <= end;  i++ )
	{
		t = x[i] - sxoss;
		st2 += t*t;
		*b += t*y[i];
	}

	*b /= st2;
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0 + sx*sx/(ss*st2))/ss);
	*sigb = sqrt(1.0/st2);
	*chi2 = 0.0;

	for( i = strt;  i <= end;  i++ ) {
		*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
	}
	sigdat = sqrt((*chi2)/(end-2));
	*siga *= sigdat;
	*sigb *= sigdat;
}


static int  install_vector_fitness (  )  {
   
/*VERBATIM*/
  {static int once; if (!once) { once = 1;
	install_vector_method("width", width);
	install_vector_method("xfitness", xfitness);
	install_vector_method("yfitness", yfitness);
	install_vector_method("ywinfitness", ywinfitness);
	install_vector_method("ywnscl_fitness", ywnscl_fitness);
	install_vector_method("firstpeak", firstpeak);
	install_vector_method("firstmax", firstmax);
	install_vector_method("nextpeak", nextpeak);
	install_vector_method("linfit", linfit);
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
