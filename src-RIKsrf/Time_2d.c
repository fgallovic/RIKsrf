/*---------------------------------------------------------------------------*/
/* TIME_2D:                                                                  */
/*        FINITE-DIFFERENCES COMPUTATION OF FIRST TRAVEL TIMES IN 2D.        */
/*        P.Podvin, Geophysique, Ecole des Mines de Paris,                   */
/*        Fontainebleau, France. October 1990.                               */
/*        Last revision : 6 May 1993                                         */
/*        (E-mail: PODVIN@GEOPHY.ENSMP.FR)                                   */
/*                                                                           */
/* USAGE :  ( From a program written in C as well as in FORTRAN )            */
/*                                                                           */
/*     (int) time_2d(HS,T,NX,NY,XS,YS,EPS_INIT,MESSAGES)                     */
/*                                                                           */
/* Note: Fortran interface is obtained through C-function time_2d_().        */
/*       Sun Fortran compiler appends the missing '_' to Fortran             */
/*       function name 'time_2d'. This convention is not standardized,       */
/*       so that full name 'time_2d_' may be required on other systems,      */
/*       according to C/FORTRAN interface conventions.                       */
/*       FORTRAN interface is designed in accordance with standard           */
/*       FORTRAN conventions. In particular, array indexes start at 1,       */
/*       and multidimensional arrays appear transposed with respect to their */
/*       equivalent in C (tab(i,j) in Fortran "means" tab[j][i] in C).       */
/*       Apart from values T(*), all other arguments are left unchanged.     */
/*       Type equivalences: int=INTEGER*4, float=REAL*4                      */
/*       In case of problems :                                               */
/*       Use compile-option -DDEBUG_ARGS to see what function you actually   */
/*       called from Fortran, (it may be NOT the one you think !) and what   */
/*       arguments were effectively passed to it !                           */
/*                                                                           */
/* Specific compile notes for IBM SP2 :                                      */
/*    1) Delete the line "extern char *malloc();".                           */
/*    2) Use at least xlc-options : -DNO_IEEE_INFINITY -langlvl=classic      */
/*    3) Calls from Fortran programs should directly use integer function    */
/*       time_2d_(), NOT time_2d() (xlf does NOT append a '_' to             */
/*       C-function names, although it may be instructed to with some        */
/*       option I don't know...)                                             */
/*                                                                           */
/* ARGUMENTS (C description; all FORTRAN arguments are pointers)             */
/*                                                                           */
/*       (int)     NX,NY           : dimensions of the time field.           */
/*                                   (number of grid points)                 */
/*                                                                           */
/*       (float *) HS,T            : 1D arrays of NX*NY elements.            */
/*                                   T will contain the computed time field  */
/*                                   arranged as a succession of columns     */
/*                                   (x=ct,y=0,NY-1). HS contains model      */
/*                                   data (slowness*mesh spacing) organized  */
/*                                   in the same manner. Meshes located at   */
/*                                   coordinates x=NX-1 or y=NY-1 are        */
/*                                   dummy meshes (out of the real model).   */
/*                                   The corresponding values will be        */
/*                                   treated as IEEE (float)infinity().      */
/*                                                                           */
/*       (float)   XS,YS         :   Point source coordinates.               */
/*                                   Licit range [0,NX-1][0,NY-1] from C,    */
/*                                   but [1,NX][1,NY] from FORTRAN (cf.      */
/*                                   indexation standards).                  */
/*                                   If source point is found to be out      */
/*                                   of the licit range, the timefield is    */
/*                                   treated as already initialized in the   */
/*                                   calling program (multiple source).      */
/*                                                                           */
/*       (float)   EPS_INIT      :   tolerance on relative inhomogeneity     */
/*                                   used (only) during initialization.      */
/*                                   Licit values in [0,1].                  */
/*                                   Relevant values should be <0.01.        */
/*                                                                           */
/*       (int)     MESSAGES        : A few messages are printed if non-zero. */
/*                                                                           */
/* VALUE :                                                                   */
/*                                                                           */
/*       time_2d() returns a nonzero value if an error was detected.         */
/*       An explicit message is printed on 'stderr'.                         */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define    M_SQRT2     1.41421356237309504880
#define    M_PI        3.14159265358979323846
#ifdef     NO_IEEE_INFINITY
#define    INFINITY    1.e+06
/* Note: INFINITY must be lower than square root of the highest   */
/* acceptable float value (machine dependent) to prevent overflow */
#define    ISINF(x)    ((x)>1.0e+06)
#define    NINT(x)     (int)floor((x)+0.5)
/* NINT should strictly be "nearest integer" */
#else
#define    INFINITY    (float)infinity()
#define    ISINF(x)    isinf(x)
#define    NINT(x)     nint(x)
#endif
#define    min(x,y)    ( ((x)<(y))? (x):(y) )
#define    max(x,y)    ( ((x)>(y))? (x):(y) )
/*extern char *malloc();*/

static int pre_init();

/*------------------------------------------------Static Functions----------*/

static void
    error(),
    init_nearest(),
    propagate_point(),
    send_x_headwave(),send_y_headwave();

static int
    init_point(),
    recursive_init(),
    x_side(),y_side();

/*------------------------------------------------Static Declarations-------*/

/* M0DEL DESCRIPTION */

static    int       nmesh_x,nmesh_y;        /* Model dimensions            */
static    float    *hs_buf,**hs;            /* 1D and 2D slowness arrays   */
static    float    *hs_keep;                /* to save boundary values     */
static    float     eps_init;

/* TIMEFIELD DESCRIPTION */

static    int       nx,ny;                  /* Timefield dimensions        */
static    float    *t_buf,**t;              /* 1D and 2D time arrays       */

/* SOURCE */

static    float     fxs,fys;                /* Point source coordinates    */
static    int       xs,ys;                  /* Nearest node                */

/* STATIC PARAMETERS */

#ifndef INIT_MIN
#define INIT_MIN        10              /* Initialize by re-discretized FD  */
                                        /* if radius of uniform region near */
                                        /* the source < INIT_MIN meshes.    */ 
                                        /* (ADJUSTABLE).                    */
#endif /* !INIT_MIN */
#define N_INIT_X    (4*INIT_MIN+3)
#define N_INIT      N_INIT_X*N_INIT_X

#ifndef INIT_RECURS_LIMIT
#define INIT_RECURS_LIMIT    1          /* Maximum level of recursivity     */
                                        /* during initialization.           */
                                        /* (ADJUSTABLE).                    */
#endif /* !INIT_RECURS_LIMIT */

static int  init_stage=0,               /* recursivity of initialization    */
            source_at_node=0,
            reverse_order=0,            /* recursivity of side computation  */
            current_side_limit,         /* current boundary of computations */
            X0,X1,Y0,Y1,                /* incl. boundaries of computed box */
            messages;                   /* message flag: 0: silent
                                                        !0: verbose.        */

/*----------------------------------------------- ERROR FLAGS ------------*/

#define NO_ERROR      0
#define ERR_INIT    (-1)
#define ERR_MULT    (-2)
#define ERR_MALLOC  (-3)
#define ERR_RECURS  (-4)
#define ERR_EPS     (-5)

static char    *err_msg[]=
            {
                "\ntime_2d: Computations terminated normally.\n",
                "\ntime_2d: Source point is in a zero velocity zone.\n",
                "\ntime_2d: Multiple source but no source at finite time.\n",
                "\ntime_2d: Memory Allocation failed.\n",
                "\ntime_2d: Fatal error during recursive init.\n",
                "\ntime_2d: [Init] Illegal tolerance on inhomogeneity.\n",
            };

/*---------------------------------------------------Error()---------------*/

static void
error(flag)
int    flag;
{    
    if(messages || flag) fprintf(stderr,"%s",err_msg[-flag]);
}

#define VERBOSE messages>1

/*---------------------------------------------------TIME_2D()-------------*/

int
time_2d(HS,T,NX,NY,XS,YS,EPS_INIT,MESSAGES)

float     *HS,*T,XS,YS,EPS_INIT;
int        NX,NY,MESSAGES;

{    
    int    i,j,signal;
    float  *pf;

#ifdef DEBUG_ARGS
    fprintf(stderr,"Entering time_2d() in C-style.\n");
    fprintf(stderr,"Args: NX=%d NY=%d XS=%g YS=%g EPS_INIT=%g MESSAGES=%d\n",
                    NX,NY,XS,YS,EPS_INIT,MESSAGES);
    fflush(stderr);
    fprintf(stderr,"First elements of input arrays: *HS=%g, *T=%g\n",*HS,*T);
#endif

    hs_buf=HS;
    t_buf=T;
    nx=NX;
    ny=NY;
    fxs=XS;
    fys=YS;
    eps_init=EPS_INIT;
    messages=MESSAGES;
    if(eps_init<0.0 || eps_init>1.0){
        error(ERR_EPS);
        return(ERR_EPS);
    }

    if((signal=pre_init())!=NO_ERROR){
        error(signal);
        return(signal);
    }
    if((signal=init_point())==NO_ERROR) propagate_point();

    if(!init_stage){
        for(i=0,pf=hs_keep;i<nx;i++,pf++) hs[i][nmesh_y]= *pf;
        for(j=0;j<nmesh_y;j++,pf++) hs[nmesh_x][j]= *pf;
        free((char *)hs_keep);
    }
    free((char *)hs);
    free((char *)t);
    error(signal);
    return(signal);
}

/*---------------------------------------------------TIME_2D_()-------------*/
/*---------------------------------------------FORTRAN INTERFACE------------*/
 
int
time_2d_(HS,T,NX,NY,XS,YS,EPS_INIT,MESSAGES)
 
float   *HS,*T,*XS,*YS,*EPS_INIT;
int     *NX,*NY,*MESSAGES;
 
{   
    int    i,j,signal;
    float  *pf;

#ifdef DEBUG_ARGS
    fprintf(stderr,"Entering time_2d() in FORTRAN-style.\n");
    fprintf(stderr,"Args: NFAST=%d NSLOW=%d S_FAST=%g S_SLOW=%g EPS_INIT=%g MESSAGES=%d\n",
                    *NX,*NY,*XS,*YS,*EPS_INIT,*MESSAGES);
    fflush(stderr);
    fprintf(stderr,"First elements of input arrays: *HS=%g, *T=%g\n",*HS,*T);
#endif
 
    hs_buf=HS;
    t_buf=T;
    nx= *NY;
    ny= *NX;
    fys= *XS-1; /* Fortran->C standard indexation */ 
    fxs= *YS-1;
 /*   fxs = 0.;
    fys = 0.; */
    eps_init= *EPS_INIT;
    messages= *MESSAGES;
    if(eps_init<0.0 || eps_init>1.0){
        error(ERR_EPS);
        return(ERR_EPS);
    }
 
    if((signal=pre_init())!=NO_ERROR){
        error(signal);
        return(signal);
    }
    if((signal=init_point())==NO_ERROR) propagate_point();

    if(!init_stage){
        for(i=0,pf=hs_keep;i<nx;i++,pf++) hs[i][nmesh_y]= *pf;
        for(j=0;j<ny-1;j++,pf++) hs[nmesh_x][j]= *pf;
        free((char *)hs_keep);
    }
    free((char *)hs);
    free((char *)t);
    error(signal);
    return(signal);
}

/*---------------------------------------------------Pre_init()------------*/

static int
pre_init()

{    
    int      i,j;
    float    *pf;

    nmesh_x=nx-1;
    nmesh_y=ny-1;

/* allocate 2D arrays */
    if(!(hs=(float **)malloc((unsigned)nx*sizeof(float *))))
        return(ERR_MALLOC);
    if(!(t=(float **)malloc((unsigned)nx*sizeof(float *)))){
        free((char *)hs);
        return(ERR_MALLOC);
    }
    for(i=j=0;i<nx;i++,j+=ny) {
        hs[i]=hs_buf+j;
        t[i]=t_buf+j;
    }

    if(init_stage) return(NO_ERROR);

/* assign zero velocity borders to the dummy meshes */
/* and keep masked values in hs_keep[].             */
    if(!(hs_keep=(float *)malloc((unsigned)(nx+nmesh_y)*sizeof(float )))){
        free((char *)hs);
        free((char *)t);
        return(ERR_MALLOC);
    }
    for(i=0,pf=hs_keep;i<nx;i++,pf++){
        *pf=hs[i][nmesh_y];
        hs[i][nmesh_y]=INFINITY;
    }
    for(j=0;j<nmesh_y;j++,pf++){
        *pf=hs[nmesh_x][j];
        hs[nmesh_x][j]=INFINITY;
    }

    return(NO_ERROR);
}

/*---------------------------------------------------Init_point()----------*/

static int
init_point()

{
    int     x,y,signal,ntime,
            test,test_X0,test_X1,test_Y0,test_Y1;
    float   min_t,hs0,sq_dist,allowed_delta_hs,*pf;
 
/* test relevance of source position or locate minimum time source point */
    if(fxs>=0.0 && fxs<=nmesh_x && fys>=0.0 && fys<=nmesh_y){
/* valid single source */
/* if first pass, assign infinity to all times */
        ntime=nx*ny;
        if(!init_stage)
            for(x=0,pf=t_buf;x<ntime;x++) *pf++=INFINITY;
 
        xs=NINT(fxs);
        ys=NINT(fys);
        if(xs==fxs && ys==fys) source_at_node=1;
        if(messages)
            printf("\nInitLevel%d: Point source [%g,%g]. Nearest node [%d,%d].",                        init_stage,fxs,fys,xs,ys);
    }
    else {
/* multiple source */
        for(x=0,min_t=INFINITY;x<nx;x++)
            for(y=0;y<ny;y++)
                if(t[x][y]<min_t){
                    min_t=t[x][y];
                    xs=x;
                    ys=y;
                }
        if(ISINF(min_t)) return(ERR_MULT);
        source_at_node=1;
        X0=X1=xs;
        Y0=Y1=ys;
        if(messages)
            printf("\nInitLevel%d: Multiple source starting at node [%d,%d],\
 at time %g.",init_stage,xs,ys,min_t);
        return(NO_ERROR);
    }

/* test relevance of model properties around the source */
    if(source_at_node){
        hs0=hs[xs][ys];
        if(ISINF(hs0) && xs) hs0=hs[xs-1][ys];
        if(ISINF(hs0) && xs && ys) hs0=hs[xs-1][ys-1];
        if(ISINF(hs0) && ys) hs0=hs[xs][ys-1];
    }
    else{
        x=(fxs<xs) ? xs-1:xs;
        y=(fys<ys) ? ys-1:ys;
        hs0=hs[x][y];
        if(ISINF(hs0) && fxs==xs && xs){
            hs0=hs[x-1][y];
            if(ISINF(hs0) && fys==ys && ys) hs0=hs[x-1][y-1];
        }
        if(ISINF(hs0) && fys==ys && ys) hs0=hs[x][y-1];
    }
    if(ISINF(hs0)) return(ERR_INIT);
    /* Source point is located in a zero velocity zone ! */

/* search for the largest square box with constant slowness */
    X0=max(xs-1,0);
    X1=min(xs+1,nmesh_x-1);
    Y0=max(ys-1,0);
    Y1=min(ys+1,nmesh_y-1);
    allowed_delta_hs=hs0*eps_init;
    test_X0=test_X1=test_Y0=test_Y1=0;
    do{
        test=0;
        if(X0 && !test_X0){
            test++;
            x=X0;
            for(y=Y0;y<=Y1 && y<nmesh_y && !test_X0;y++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_X0=1;
            if(!test_X0) X0--;
        }
        if(Y0 && !test_Y0){
            test++;
            y=Y0;
            for(x=X0;x<=X1 && x<nmesh_x && !test_Y0;x++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_Y0=1;
            if(!test_Y0) Y0--;
        }
        if(X1<nmesh_x && !test_X1){
            test++;
            X1++;
            x=X1;
            for(y=Y0;y<=Y1 && y<nmesh_y && !test_X1;y++)
                if(fabs(hs[x][y]-hs0)>allowed_delta_hs) test_X1=1;
        }
        if(Y1<nmesh_y && !test_Y1){
            test++;
            Y1++;
            y=Y1;
            for(x=X0;x<=X1 && x<nmesh_x && !test_Y1;x++)
                if(fabs(hs[x][Y1]-hs0)>allowed_delta_hs) test_Y1=1;
        }
    } while(test);

/* decrement boundaries of homogeneous region (if not at model boundaries) */
/* so that heterogeneous interfaces are dealt with by the FD scheme.       */
    if(X0) X0++;
    if(Y0) Y0++;
    if(X1<nmesh_x) X1--;
    if(Y1<nmesh_y) Y1--;

/* initialize the time-field according to situation... */
    if(init_stage>=INIT_RECURS_LIMIT ||
        (   (X0==0 || (xs-X0)>=INIT_MIN) &&
            (Y0==0 || (ys-Y0)>=INIT_MIN) &&
            (X1==nmesh_x || (X1-xs)>=INIT_MIN) &&
            (Y1==nmesh_y || (Y1-ys)>=INIT_MIN)
        ) ) {
        if((X1-X0+1)*(Y1-Y0+1)==1) init_nearest();
        else {
            for(x=X0;x<=X1;x++)
                for(y=Y0;y<=Y1;y++) {
                    sq_dist=(x-fxs)*(x-fxs)+(y-fys)*(y-fys);
                    t[x][y]=hs0*sqrt(sq_dist);
                }
            if(!init_stage && messages) 
                printf("\nHomogeneous region: x(%d->%d),y(%d->%d)",X0,X1,Y0,Y1);
        }
        signal=NO_ERROR;
    }
    else{
        if((signal=recursive_init())!=NO_ERROR) return(signal);
        X0=max(xs-INIT_MIN,0);
        Y0=max(ys-INIT_MIN,0);
        X1=min(xs+INIT_MIN,nmesh_x);
        Y1=min(ys+INIT_MIN,nmesh_y);
    }

    return(signal);

}

/*----------------------------------------------------init_nearest()---------*/
/* last changes in this function 060593; could fail for sources located at */
/* model boundaries */

static void
init_nearest()

/* initialize the 1|4|6 nearest nodes when model is immediately inhomogeneous */

{
    int x,y;
    float dx,dy,hs0,hs1,d;

    if(messages) printf("\nInitializing closest node(s).");
    if(source_at_node){
        t[xs][ys]=0.0;
        return;
    }
    x=(fxs<xs) ? xs-1:xs;
    y=(fys<ys) ? ys-1:ys;
    hs0=hs[x][y];
    dx=fabs(fxs-xs);
    dy=fabs(fys-ys);
    if(xs==fxs){
        if(x) hs1=hs[x-1][y];
        else hs1=INFINITY;
        t[xs][ys]=min(hs0,hs1)*dy;
        t[xs][y]=min(hs0,hs1)*(1.0-dy);
        d=sqrt(1.0+dy*dy);
        t[x][ys]=hs0*d;
        if(x) t[x-1][ys]=hs1*d;
        d=sqrt(1.0+(1.0-dy)*(1.0-dy));
        t[x][y]=hs0*d;
        if(x) t[x-1][y]=hs1*d;
    }/* source is located on a x=Ct interface */
    else if(ys==fys){
        if(y) hs1=hs[x][y-1];
        else hs1=INFINITY;
        t[xs][ys]=min(hs0,hs1)*dx;
        t[x][ys]=min(hs0,hs1)*(1.0-dx);
        d=sqrt(1.0+dx*dx);
        t[xs][y]=hs0*d;
        if(y) t[xs][y-1]=hs1*d;
        d=sqrt(1.0+(1.0-dx)*(1.0-dx));
        t[x][y]=hs0*d;
        if(y) t[x][y-1]=hs1*d;
    }/* source is located on a y=Ct interface */
    else {
        t[xs][ys]=hs0*sqrt(dx*dx+dy*dy);
        t[x][ys]=hs0*sqrt((1.0-dx)*(1.0-dx)+dy*dy);
        t[xs][y]=hs0*sqrt((1.0-dy)*(1.0-dy)+dx*dx);
        t[x][y]=hs0*sqrt((1.0-dx)*(1.0-dx)+(1.0-dy)*(1.0-dy));
    }/* source is located within a cell */
}

/*----------------------------------------------------recursive_init()-------*/

static int
recursive_init()

/* time_2d() is used recursively to initialize the timefield on a */
/* small region around the source with a smaller mesh spacing.    */

{
    int     signal,
            nx_,ny_,
            xs_,ys_,
            n,d,
            i,ii,ihs,i0,
            j,jj,jhs,j0;
    float   
            fxs_,fys_,
            *hs_buf_,*t_buf_,
            HS[N_INIT],T[N_INIT];

/* save static parameters at this stage */
    nx_=nx;
    ny_=ny;
    hs_buf_=hs_buf;
    t_buf_=t_buf;
    xs_=xs;
    ys_=ys;
    fxs_=fxs;
    fys_=fys;

/* increment count of recursivity level */
    init_stage++;
    if(messages)
        printf("\nRecursive initialization: level %d",init_stage);

/* free (float **) pointers */
    free((char *)hs);
    free((char *)t);

/* build the re-discretized local model */
    for(i=0;i<N_INIT;i++) HS[i]=T[i]=INFINITY;
    nx=ny=N_INIT_X;
    xs=ys=2*INIT_MIN+1;
    i0=j0=1;
    ihs=xs_-INIT_MIN-1;
    if((d=INIT_MIN-xs_)>=0){
        ihs+=d+1;
        d=1+2*d;
        nx-=d;
        xs-=d;
        i0=0;
    }
    if((d=xs_+INIT_MIN-nx_+1)>=0) nx-=1+2*d;
    jhs=ys_-INIT_MIN-1;
    if((d=INIT_MIN-ys_)>=0){
        jhs+=d+1;
        d=1+2*d;
        ny-=d;
        ys-=d;
        j0=0;
    }
    if((d=ys_+INIT_MIN-ny_+1)>=0) ny-=1+2*d;
    for(i=ihs,n=ii=0;ii<nx;ii++){
        for(j=jhs,jj=0;jj<ny;jj++,n++){
            HS[n]=0.5*hs_buf_[i*ny_+j];
            if(jj%2!=j0) j++;
        }
        if(ii%2!=i0) i++;
    }/* No smoothing is associated with this re-discretization */

/* compute new source coordinates */
    fxs=xs+2.0*(fxs_-xs_);
    fys=ys+2.0*(fys_-ys_);

/* recursively compute times on this model */
    signal=time_2d(HS,T,nx,ny,fxs,fys,eps_init,messages);

/* assign relevant times to final timefield */
    if(signal==NO_ERROR){
        for(i=ihs+i0,ii=i0;ii<nx;ii+=2,i++)
            for(j=jhs+j0,jj=j0;jj<ny;jj+=2,j++){
                t_buf_[i*ny_+j]=T[ii*ny+jj];
            }
    }
    else error(signal);

/* retrieve old static parameters */
    nx=nx_;
    ny=ny_;
    hs_buf=hs_buf_;
    t_buf=t_buf_;
    xs=xs_;
    ys=ys_;
    fxs=fxs_;
    fys=fys_;

/* reallocate (float **) pointers but do not re-initialize ! */
    i=pre_init();
    signal=(signal==NO_ERROR)?i:ERR_RECURS;

/* decrement count of recursivity level */
    init_stage--;

    return signal;
}

/*---------------------------------------------------Propagate_point()-----*/

static void
propagate_point()

{
    int msg,
        test;

/* make recursive init silent */
    msg=messages;
    if(init_stage) messages=0;

    do {
        test=0;
        if(X0>0){
            X0--;
            if(VERBOSE) printf("\nside x=%d.",X0);
            y_side(X0,-1,Y0,Y1);
            test++;
        }
        if(Y0>0){
            Y0--;
            if(VERBOSE) printf("\nside y=%d.",Y0);
            x_side(Y0,-1,X0,X1);
            test++;
        }
        if(X1<nmesh_x){
            X1++;
            if(VERBOSE) printf("\nside x=%d.",X1);
            y_side(X1,1,Y0,Y1);
            test++;
        }
        if(Y1<nmesh_y){
            Y1++;
            if(VERBOSE) printf("\nside y=%d.",Y1);
            x_side(Y1,1,X0,X1);
            test++;
        }

    } while(test);

    messages=msg;

}

/*---------------------------------------------------y_side()---------------*/

static int
y_side(x,future,y_begin,y_end)

int    x,future,y_begin,y_end;

/* propagates computations from row x-future to row x */
/* between inclusive y_begin and y_end coordinates.   */
/* y_side() returns the number of stencils adopted.   */

{
    int
        x0,          /* past side coordinate          */
        x_s0,        /* current mesh coordinate       */
        y,           /* current point coordinate      */
        y_old,       /* last local minimum coordinate */
        past,        /* ...opposite to future !       */
        updated,     /* counts accepted stencils      */
        longhead,    /* counts longitudinal headwaves */
        alert;       /* longhead is alert ? (0/1)     */
    float
        hs0,hs1,hs2, /* current slownesses            */
        t_est,       /* current time estimate         */
        dt;          /* time difference               */

    if(reverse_order==0) current_side_limit=x+future;

    updated=longhead=0;
    x0=x-future;
    if(future==1) x_s0=x0;
    else x_s0=x;

    for(y=y_begin;y<=y_end;){

/* Search the next local time minimum on row x-future, save its position */
        while(y<y_end && t[x0][y+1]<t[x0][y]) y++;
        y_old=y;

/* Compute time in front of this local minimum */
        hs1=hs[x_s0][y];
        if(y==0) hs0=INFINITY;
        else hs0=hs[x_s0][y-1];
        if((t_est=t[x0][y]+min(hs0,hs1))<t[x][y]) {
            t[x][y]=t_est;
            updated++;
        }/* timed by 1D transmission */

/* Proceed backwards until the last local time maximum (or y_begin) */
        y--;
        alert=0;
        while(y>=y_begin && (dt=t[x0][y]-t[x0][y+1])>=0.0) {
            hs0=hs[x_s0][y];
            if(dt<hs0/M_SQRT2 
                && (t_est=t[x0][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            dt=t[x][y+1]-t[x0][y+1];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x][y+1]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            if(y!=0) {
                hs1=hs[x_s0][y-1];
                if((t_est=t[x0][y]+hs1)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                }/* 1D transmission towards future */
            }
            if((t_est=t[x0][y+1]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(x_s0+future>=0) {
                hs2=hs[x_s0+future][y];
                if((t_est=t[x][y+1]+hs2)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                    if(!alert){
                        longhead++;
                        send_y_headwave(y,y_begin,
                            x,x_s0,x_s0+future);
                        alert=1;
                    }
                }/* 1D transmission along the current side. */
                /* if adopted, implies reverse propagation  */
                else alert=0;
            }
            y--;
        }

/* Proceed forwards until the next local time maximum (or y_end) */
        if((y=y_old)==y_end) break;
        y++;
        alert=0;
        while(y<=y_end && (dt=t[x0][y]-t[x0][y-1])>=0.0) {
            hs0=hs[x_s0][y-1];
            if(dt<hs0/M_SQRT2
                && (t_est=t[x0][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            dt=t[x][y-1]-t[x0][y-1];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x][y-1]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            hs1=hs[x_s0][y];
            if((t_est=t[x0][y]+hs1)<t[x][y]){
                t[x][y]=t_est;
                updated++;
            }/* 1D transmission towards future */
            if((t_est=t[x0][y-1]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(x_s0+future>=0) {
                hs2=hs[x_s0+future][y-1];
                if((t_est=t[x][y-1]+hs2)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                    if(!alert){
                        longhead++;
                        send_y_headwave(y,y_end,
                            x,x_s0,x_s0+future);
                        alert=1;
                    }
                }/* 1D transmission along the current side. */
                /* if adopted, implies reverse propagation  */
                else alert=0;
            }
            y++;
        }
    }

/* times are now computed at every point on the current side.  */
/* Reverse propagation may be necessary if a headwave has been */
/* generated along this side...                                */

    if(longhead) {

        reverse_order++;
        if(messages) 
            printf("\nReverse-propagation(order %d) from row %d:",
                reverse_order,x);
        past= -future;
        for(x=x0;x!=current_side_limit;x+=past) {
            if(x<0 || x>nmesh_x) break; /* The Horrific Case ! (see Note) */
            if(VERBOSE) printf("\nrow #%d: ",x);
            if(y_side(x,past,y_begin,y_end)==0) break;
            if(messages) printf(" (%d):updated.",x);
        }
        reverse_order--;

    }

    return(updated);

}
/* Note(May 1993); New simple tests were added to prevent The Horrific Case. */
/* The Horrific Case may arise in very awkward (non-connex) models when      */
/* several sources are simulated in areas separated by zero-velocity regions */ 
/* and if one of these regions of propagation reaches model boundaries with  */
/* a particular geometry... Whaow !...                                       */
/* This case WAS encountered in a simulation of rupture along a fault-plane. */

/*---------------------------------------------------send_y_headwave()------*/

static void
send_y_headwave(from,to,x,x_s_now,x_s_future)

int    from,to,x,x_s_now,x_s_future;

/* enforces propagation of a headwave to the end of the current row.   */
/* (this is only needed in very severe models where such headwave may  */
/* provide first arrival at points located FAR from the region where   */
/* critical conditions were encountered...e.g., the slow disk example).*/
/* These headwaves may "jump" over local maxima of the past side.      */
/* In such cases, this information would be lost.                      */

{
    int      y;
    float    hsnow,hsmin,t_est;

    if(from<to) for(y=from;y<to;y++){
        hsnow=hs[x_s_now][y];
        if(x_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x_s_future][y]);
        if((t_est=t[x][y]+hsmin)<t[x][y+1]) t[x][y+1]=t_est;
    }
    else for(y=from;y>to;y--){
        hsnow=hs[x_s_now][y-1];
        if(x_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x_s_future][y-1]);
        if((t_est=t[x][y]+hsmin)<t[x][y-1]) t[x][y-1]=t_est;
    }
}

/*---------------------------------------------------x_side()---------------*/

static int
x_side(y,future,x_begin,x_end)

int    y,future,x_begin,x_end;

/* propagates computations from row y-future to row y */
/* between inclusive x_begin and x_end coordinates.   */
/* x_side() returns the number of stencils adopted.   */

{
    int
        y0,          /* past side coordinate          */
        y_s0,        /* current mesh coordinate       */
        x,           /* current point coordinate      */
        x_old,       /* last local minimum coordinate */
        past,        /* ...opposite to future !       */
        updated,     /* counts accepted stencils      */
        longhead,    /* counts longitudinal headwaves */
        alert;       /* longhead is alert ? (0/1)     */
    float
        hs0,hs1,hs2, /* current slownesses            */
        t_est,       /* current time estimate         */
        dt;          /* time difference               */

    if(reverse_order==0) current_side_limit=y+future;

    updated=longhead=0;
    y0=y-future;
    if(future==1) y_s0=y0;
    else y_s0=y;

    for(x=x_begin;x<=x_end;){

/* Search for the next local time minimum on row y-future, save its position */
        while(x<x_end && t[x+1][y0]<t[x][y0]) x++;
        x_old=x;

/* Compute time in front of this local minimum */
        hs1=hs[x][y_s0];
        if(x==0) hs0=hs1;
        else hs0=min(hs1,hs[x-1][y_s0]);
        if((t_est=t[x][y0]+hs0)<t[x][y]) {
            t[x][y]=t_est;
            updated++;
        }/* timed by 1D transmission */

/* Proceed backwards until the last local time maximum (or x_begin) */
        x--;
        alert=0;
        while(x>=x_begin && (dt=t[x][y0]-t[x+1][y0])>=0.0) {
            hs0=hs[x][y_s0];
            if(dt<hs0/M_SQRT2 
                && (t_est=t[x][y0]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            dt=t[x+1][y]-t[x+1][y0];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x+1][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            if(x!=0) {
                hs1=hs[x-1][y_s0];
                if((t_est=t[x][y0]+hs1)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                }/* 1D transmission towards future */
            }
            if((t_est=t[x+1][y0]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(y_s0+future>=0) {
                hs2=hs[x][y_s0+future];
                if((t_est=t[x+1][y]+hs2)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                    if(!alert){
                        longhead++;
                        send_x_headwave(x,x_begin,
                            y,y_s0,y_s0+future);
                        alert=1;
                    }
                }/* 1D transmission along the current side. */
                /* if adopted, implies reverse propagation  */
                else alert=0;
            }
            x--;
        }

/* Proceed forwards until the next local time maximum (or x_end) */
        if((x=x_old)==x_end) break;
        x++;
        alert=0;
        while(x<=x_end && (dt=t[x][y0]-t[x-1][y0])>=0.0) {
            hs0=hs[x-1][y_s0];
            if(dt<hs0/M_SQRT2
                && (t_est=t[x][y0]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through horizontal interface */
            dt=t[x-1][y]-t[x-1][y0];
            if(dt>=0.0 && dt<hs0/M_SQRT2
                && (t_est=t[x-1][y]+sqrt(hs0*hs0-dt*dt))<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* 2D transmission through vertical interface */
            hs1=hs[x][y_s0];
            if((t_est=t[x][y0]+hs1)<t[x][y]){
                t[x][y]=t_est;
                updated++;
            }/* 1D transmission towards future */
            if((t_est=t[x-1][y0]+hs0*M_SQRT2)<t[x][y]) {
                t[x][y]=t_est;
                updated++;
            }/* diffraction */
            if(y_s0+future>=0) {
                hs2=hs[x-1][y_s0+future];
                if((t_est=t[x-1][y]+hs2)<t[x][y]){
                    t[x][y]=t_est;
                    updated++;
                    if(!alert){
                        longhead++;
                        send_x_headwave(x,x_end,
                            y,y_s0,y_s0+future);
                        alert=1;
                    }
                }/* 1D transmission along the current side. */
                /* if adopted, implies reverse propagation  */
                else alert=0;
            }
            x++;
        }
    }

/* times are now computed at every point on the current side.  */
/* Reverse propagation may be necessary if a headwave has been */
/* generated along this side...                                */

    if(longhead) {

        reverse_order++;
        if(messages) 
            printf("\nReverse-propagation(order %d) from row %d:",
                reverse_order,y);
        past= -future;
        for(y=y0;y!=current_side_limit;y+=past) {
            if(y<0 || y>nmesh_y) break; /* The Horrific Case ! (see Note) */
            if(VERBOSE) printf("\nrow #%d: ",y);
            if(x_side(y,past,x_begin,x_end)==0) break;
            if(messages) printf(" (%d):updated.",y);
        }
        reverse_order--;

    }

    return(updated);

}

/*---------------------------------------------------send_x_headwave()------*/

static void
send_x_headwave(from,to,y,y_s_now,y_s_future)

int    from,to,y,y_s_now,y_s_future;

/* enforces propagation of a headwave to the end of the current row.   */
/* (this is only needed in very severe models where such headwave may  */
/* provide first arrival at points located FAR from the region where   */
/* critical conditions were encountered...e.g., the slow disk example).*/
/* These headwaves may "jump" over local maxima of the past side.      */
/* In such cases, this information would be lost.                      */

{
    int    x;
    float  hsnow,hsmin,t_est;

    if(from<to) for(x=from;x<to;x++){
        hsnow=hs[x][y_s_now];
        if(y_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x][y_s_future]);
        if((t_est=t[x][y]+hsmin)<t[x+1][y]) t[x+1][y]=t_est;
    }
    else for(x=from;x>to;x--){
        hsnow=hs[x-1][y_s_now];
        if(y_s_future<0) hsmin=hsnow;
        else hsmin=min(hsnow,hs[x-1][y_s_future]);
        if((t_est=t[x][y]+hsmin)<t[x-1][y]) t[x-1][y]=t_est;
    }
}
/*----------------------------------------------END: TIME_2D----------------*/
