// 8 bins in play are 1,3,7,9,10,11,12,14
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <fcntl.h>
#include <sys/soundcard.h>
#include <alsa/asoundlib.h>
#include <speex/speex_preprocess.h>


#define MAXCB 10000 // maximum words times average word size

#define WINDOW 320 
#define STEP 80 // for 50fps at 8K
#define LSP_MAX_ORDER 10
#define ORDER 10
#define FFT 32
#define FF2 17
int bin[FF2]={0,16,8,24,4,20,11,19,2,18,10,21,5,25,9,17,1};

#define TILT 20.0

#define ETOPBIT 1024 // 2^10, i.e. energy buffer size is 10 samples
#define TOL_ON  450000 // senn=1000000
#define TOL_OFF 152500 // senn=550000
#define DELTA -20000
#define EMPH 0.0
#define FLOOR 100.0
#define TRACE 0

// globals...
char report[8];
float realtwiddle[16],imtwiddle[16];
float cb[MAXCB][FF2];
int start[MAXCB],size[MAXCB],cbsize=0;
float word[1000][FF2];
int wsize=0;
int training=0,count=0,on=0;

char target[37]="thequickbrownfoxoxjumpsoverthelazydog";

float lpc[ORDER+1],last[FF2],last2[FF2];

//----------------------------------------------------------------------
void setup_twiddles(float *realtwiddle, float *imtwiddle, int N)
{
    int n;

    for(n=1; n<N/2; n++)
    {
       realtwiddle[n] = cos(2.0*M_PI*n/N);
       imtwiddle[n] = -sin(2.0*M_PI*n/N);
    }
}

//-----------------------------------------------------------------------
void simple_fft(float *real, float *im, float *realtwiddle, float *imtwiddle, int N)
{
    unsigned int even, odd, span, log=0, rootindex;    // indexes
    float temp;
    int i;

    for(span=N>>1; span; span>>=1, log++)   
    {
        for(odd=span; odd<N; odd++)         // iterate over the dual nodes
        {
	    
            odd |= span;                    // iterate over odd blocks only
            even = odd ^ span;              // even part of the dual node pair
                      
            temp = real[even] + real[odd];       
            real[odd] = real[even] - real[odd];
            real[even] = temp;
           
            temp = im[even] + im[odd];           
            im[odd] = im[even] - im[odd];
            im[even] = temp;
           
            rootindex = (even<<log) & (N-1); // find root of unity index
            if(rootindex)                    // skip rootindex[0] (has an identity)
            {
                temp=realtwiddle[rootindex]*real[odd]-imtwiddle[rootindex]*im[odd];
                im[odd]=realtwiddle[rootindex]*im[odd]+imtwiddle[rootindex]*real[odd];
                real[odd] = temp;
            }
   
        } // end of loop over n
     } // end of loop over FFT stages
} //end of function


//----------------------------------------------------------------------------

void autocorrelate(
  float Sn[],	/* frame of Nsam windowed speech samples */
  float Rn[],	/* array of P+1 autocorrelation coefficients */
  int Nsam,	/* number of windowed samples to use */
  int order	/* order of LPC analysis */
)
{
  int i,j;	/* loop variables */
  float x;

  for(j=0; j<order+1; j++) {
    x = 0.0;
    for(i=0; i<Nsam-j; i++) x+=Sn[i]*Sn[i+j];
    Rn[j]=x;
  }
  if (Rn[0]<1.0) Rn[0]=1.0; // prevent div-by-0 later
}

//------------------------------------------------------
// Levinson-Durbin function based on Speex 1.0.5 lpc.c....
void wld(
    float       * lpc, /*      [0...p-1] LPC coefficients      */
    const float * ac,  /*  in: [0...p] autocorrelation values  */
    int p
    )
{
   int i, j;  float r, error = ac[0];

   for (i = 0; i < p; i++) {

      /* Sum up this iteration's reflection coefficient.
       */
      r = -ac[i + 1];
      for (j = 0; j < i; j++) r -= lpc[j] * ac[i - j];
      r /= error;

      /*  Update LPC coefficients and total error.
       */
      lpc[i] = r;
      for (j = 0; j < i/2; j++) {
         float tmp  = lpc[j];
         lpc[j]     += r * lpc[i-1-j];
         lpc[i-1-j] += r * tmp;
      }
      if (i % 2) lpc[j] += lpc[j] * r;

      error *= 1.0 - r * r;
   }
}

//---------------------------------------------------------------
// add word into the codebook, write+exit if all words trained.
int train()
{
  FILE *fp;
  int i,n,k;
  
  bcopy(word,&cb[start[cbsize]][0],wsize*FF2*4);
  size[cbsize]=wsize; start[cbsize+1]=start[cbsize]+wsize;
  fprintf(stderr,"*** TRAINED %i or '%c'\n",cbsize+1,cbsize+'a'); cbsize++;
  if (cbsize==training) {  // write codebook, exit...
    fp=fopen("cb.txt","w");
    for (i=0; i<training; i++) {
      fprintf(fp,"%i\n",size[i]);
      for (n=0; n<size[i]; n++) {
        for (k=1; k<FF2-1; k++) fprintf(fp,"%i,",(int)cb[start[i]+n][k]);
	fprintf(fp,"\n");
      }
    }
    fclose(fp);
    exit(0);
  }
}

//---------------------------------------------------------------
// search codebook, see what word[] most closely matches...
int closest()
{
   double dist,bestdist,d;
   int best=0,c,i,n,x;
   FILE *fp_hid;
   
   bestdist=9999999999999.0;
   for (c=0; c<cbsize; c++) {
     dist=0;
     for (i=0; i<wsize; i++) {
       //x=i*size[c]/wsize; if (x>=size[c]) x=size[c]-1;
       x=i; if (x>=size[c]) x=size[c]-1;
       for (n=1; n<FF2-1; n++)
        {d=word[i][n]-cb[start[c]+x][n]; dist+=d*d;}
     }
     if (dist<bestdist) {bestdist=dist; best=c;}
   }
   fprintf(stderr,"BEST=%i	or %c, dist=%.f\n",best+1, best+'a',bestdist);
   if (best==26) best=40; // space
   if (best==27) best=51; // fullstop
   if (best==28) best=36; // newline
   if (best==29) best=38; // no / backsapce
   if (best==30) {on=1; system("i2cset -f -y 0 0x34 0x93 0x1");  best=36;} // switch on.
   if (best==31) {on=0; system("i2cset -f -y 0 0x34 0x93 0x0"); best=36;} // switch off.
   if (on) {
     fp_hid=fopen("/dev/hidg0","w");
     report[2]=best+4;
     fwrite(report,8,1,fp_hid); report[2]=0; fwrite(report,8,1,fp_hid);
     fclose(fp_hid);
  }
}



//---------------------------------------------------------------

main(int argc, char *argv[])
{ 

  FILE *fp;
  int fd,arg;
  snd_pcm_t *handle;
  snd_pcm_hw_params_t *hw_params;
  int rate=8000;
  float f[WINDOW],hann[WINDOW],w[WINDOW],w2[WINDOW],s0,s1=0,tot;
  float ac[ORDER+1],lcp[ORDER+1],lsp[ORDER],l[ORDER],weight[ORDER],delta,d;
  short sample,s[160],buf[2000];
  int i,j,n,b,toggle=1;
  float e,laste=0;
  int ebit=ETOPBIT, ebuff=0;
  int sound=0; // boolean start/stop
  float f2[FFT],min;
  float real[FFT],imag[FFT];
  int dummy[100000];
  float amp[WINDOW],pha[WINDOW];
  int frame=0;
  SpeexPreprocessState *st;
 
  for (i=0; i<8; i++) report[i]=0; 
  st = speex_preprocess_state_init(160, 8000);
  i=1;
  speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DENOISE, &i);
//  i=1;
//  speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB, &i);
//  e=.0;
//  speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_DECAY, &e);
//  e=.0;
//  speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, &e);
  
  
  setup_twiddles(realtwiddle,imtwiddle,FFT);
  
  for (i=0; i<WINDOW; i++) f[i]=0;
  for (i=0; i<ORDER; i++) {last[i]=0; last2[i]=0;}
  
  for(i=0; i<WINDOW; i++) hann[i]=0.5-0.5*cos(2.0*M_PI*i/(WINDOW-1));
  
  if (argc==2) training=atoi(argv[1]);
  fprintf(stderr,"training=%i\n",training); //  exit(0);
  
  cbsize=0; start[0]=0; size[0]=0;
  if (training==0) if (fp=fopen("cb.txt","r")) {
    while(!feof(fp)) {
      fscanf(fp,"%i\n",&size[cbsize]);
      for (i=start[cbsize]; i<start[cbsize]+size[cbsize]; i++) {
        for (n=1; n<FF2-1; n++) fscanf(fp,"%f,",&cb[i][n]); fscanf(fp,"\n");
      }
      start[cbsize+1]=start[cbsize]+size[cbsize];  cbsize++; 
    }
    fclose(fp);
  }
  //for (i=0; i<cbsize; i++) printf("%i,\n",size[i]); exit(0);
  
  //---------------------------------
  
  
  fp=fopen("/tmp/b.raw","w");
  snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0);	
  snd_pcm_hw_params_malloc(&hw_params);			 
  snd_pcm_hw_params_any(handle, hw_params);
  snd_pcm_hw_params_set_access(handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);	
  snd_pcm_hw_params_set_format(handle, hw_params, SND_PCM_FORMAT_S16_LE);	
  snd_pcm_hw_params_set_rate_near(handle, hw_params, &rate, 0);
  snd_pcm_hw_params_set_channels(handle, hw_params, 2);
  snd_pcm_hw_params(handle, hw_params);
  snd_pcm_hw_params_free(hw_params);
  snd_pcm_prepare(handle);
  
  //printf("sleep 1...\n"); sleep(1); printf("OK, go....\n");
  while(1) {
    for (i=0; i<WINDOW-STEP; i++) f[i]=f[i+STEP]; // shift samples down
    if (toggle) {
      //read(fd,s,160*2);
      snd_pcm_readi(handle, buf, 160);
      for (i=0; i<160; i++) s[i]=buf[i*2];
      speex_preprocess_run(st, s);
    }
    else bcopy(&s[80],s,80*2);
    toggle=!toggle;
    for (i=WINDOW-STEP,j=0; i<WINDOW; i++,j++) {
      sample=s[j]; s0=(float)sample;
      f[i]=s0-s1*EMPH; s1=s0; // 1.0 pre-emphasis
      fwrite(&sample,2,1,fp);
  
    }
    for (i=0; i<WINDOW; i++) w[i]=f[i];
    // remove any DC level....
    tot=0; for (i=0; i<WINDOW; i++) tot+=w[i]; tot/=WINDOW;
    for (i=0; i<WINDOW; i++) w[i]-=tot;
    
    for (i=0; i<WINDOW; i++) w[i]*=hann[i]; // window data
 
    autocorrelate(w,ac,WINDOW,ORDER);
    wld(&lpc[1],ac,ORDER); lpc[0]=1.0;
        // e=ac[0];
    e=0;for(i=0; i<=ORDER; i++) e+=ac[i]*lpc[i];   if (e<0) e=0;

    if (e>TOL_OFF) ebuff|=ebit; else ebuff&=~ebit; // update energy bit-buffer
    ebit>>=1; if (ebit==0) ebit=ETOPBIT; // circular shift
 
    for (i=0; i<FFT; i++) {real[i]=0; imag[i]=0;}
    for (i=0; i<=ORDER; i++) real[i]=lpc[i];
    simple_fft(real,imag,realtwiddle,imtwiddle,FFT);
    for (i=0; i<FF2; i++) {
      b=bin[i];
      f2[i]=powf(real[b]*real[b]+imag[b]*imag[b],-0.5);
      //f2[i]=powf(f2[i],0.333333);
      //f2[i]=powf(f2[i],1.2);
      f2[i]=logf(f2[i]);
    }

   
    // spectral tilt compensation...
    for (i=1; i<FF2; i++) f2[i]=f2[i]*(float)(i+TILT)/TILT;

 
    // fold down to 9 bins...
/*
    if (f2[FF2-2]>f2[FF2-3]) f2[FF2-3]=f2[FF2-2]; f2[FF2-2]=0;
    if (f2[FF2-4]>f2[FF2-5]) f2[FF2-5]=f2[FF2-4]; f2[FF2-4]=0;
    if (f2[FF2-9]>f2[FF2-10]) f2[FF2-10]=f2[FF2-9]; f2[FF2-9]=0;
    if (f2[FF2-11]>f2[FF2-12]) f2[FF2-12]=f2[FF2-11]; f2[FF2-11]=0;
    if (f2[FF2-13]>f2[FF2-14]) f2[FF2-14]=f2[FF2-13]; f2[FF2-13]=0;
    if (f2[FF2-15]>f2[FF2-16]) f2[FF2-16]=f2[FF2-15]; f2[FF2-15]=0;
*/
    for (i=0; i<FF2; i++) {
      if (f2[i]>6.0) f2[i]=6.0;
      f2[i]*=100.0;
    }
    
if (TRACE) { fprintf(stderr,"%.f,",e);  for (i=1; i<FF2-1; i++) fprintf(stderr,"%.f,",f2[i]); fprintf(stderr,"\n");}

   
    // calculate frame delta....
    delta=0; for (i=1; i<FF2-1; i++) {d=f2[i]-last[i]; delta+=d*d;}
    //printf("delta=%f\n",delta);
 
    if (sound==0 && e>TOL_ON && frame>200)  { // start recording...
      bcopy(last2,&word[wsize],FF2*4); wsize++;
      bcopy(last,&word[wsize],FF2*4); wsize++;
      sound=1; wsize=0;
      bcopy(f2,&word[wsize],FF2*4); wsize++;
      bcopy(last,last2,FF2*4);
      bcopy(f2,last,FF2*4);
    }
    else if (sound==1 && e>TOL_OFF) { // continue reading word...
      bcopy(f2,&word[wsize],FF2*4); wsize++; if (wsize>200) wsize=200;
      bcopy(last,last2,FF2*4);
      bcopy(f2,last,FF2*4);
    }
    else if (sound==1 && ebuff==0) { // finised reading word
       // wsize-=8; // remove training silence (2 frame buffer)
       if (wsize>4 && wsize<50) {
         if (training>0) train();
         else   closest(); 
       }     
       sound=0; wsize=0; bcopy(last,last2,FF2*4); bcopy(f2,last,FF2*4);
    }
    
   //for (i=1; i<FF2-1; i++) printf("%.0f,",f2[i]); printf("  e=%f\n",e);
   laste=e;
   frame++; if (frame==37800) exit(0);
  }
}
