/*
 * This file is part of the Kalman1D library
 * Usage: A template library for the implementation
 *        of one dimensional Kalman filters for arduino/teensy
 *        with asynchronous sensor fusion
 * Version 1.0.1
 * Developed by Evan https://github.com/halsw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 * Classes:
 *   K1Filter the one dimesional Kalman filter with asynchronous data fusion
 *            the process noise here may be defined each iteration and does not assume a fixed period
 *   K1FilterLinear a cluster of one dimesional Kalman filters with asynchronous data fusion for dynamic models in the form X(t+Δt) = X(t) + X'(t)Δt+X"(t)Δt^2/2! + ...
 *            the update period defined is used to normalize the process noise to tha actual state update interval
 *   K1FilterMobile a cluster of one dimesional Kalman filters with asynchronous data fusion for holonomic mobile robots with a state vector containing
 *            0:distance travelled, 1:linear velocity, 2:linear acceleration, 3:angle travelled, 4:angular velocity, 5: angular acceleration,
 *            6:x coordinate, 7: y coordinate
 *            it also overrides the digest() method to convert a linear acceleration measurement to angular acceleration, this is not the best practice to measure
 *            angular acceleration and is implemented here for test and demonstration purposes
 */
#ifndef KALMAN1D_H
#define KALMAN1D_H

#define K1_STATE 0
#define K1_STATE_PREVIOUS 1
#define K1_STATE_ESTIMATED 2
#define K1_VARIANCE 3
#define K1_VARIANCE_PREVIOUS 4
#define K1_VARIANCE_ESTIMATED 5
#define K1_PROCESS_NOISE 6
#define K1_CONTROL 7
#define K1_PERIOD 8
#define K1_SENSORS 9
#define K1_SENSORS_PREVIOUS 10
#define K1_SENSOR_VARIANCE 11
#define K1_SENSOR_VARIANCE_PREVIOUS 12
#define K1_SENSORS_NUM 13
#define K1_SENSORS_NUM_PREVIOUS 14

#define K1STATE_DISPLACEMENT 0
#define K1STATE_VELOCITY 1
#define K1STATE_ACCELERATION 2
#define K1STATE_ANGLE 3
#define K1STATE_HEADING 3
#define K1STATE_ANGVEL 4
#define K1STATE_ANGACC 5
#define K1STATE_XCOORD 6
#define K1STATE_YCOORD 7


template <class K1> class K1FilterLinear;
template <class K1> class K1FilterMobile;

//The one dimensional Kalman Filter
template <class K1>
class K1Filter {
  friend class K1FilterLinear<K1>;
  friend class K1FilterMobile<K1>;
  protected:
    K1Filter *parent, *child;  
    unsigned long tp;
    K1 xPrv, xEst, xCur;
    K1 vPrv, vEst, vCur;
    K1 gain, xVar, xCon;
    K1 input0, input1;
    K1 invar, invar0, invar1;
    K1 dt;
    unsigned long intim0, intim1;  
    unsigned char sensors0, sensors1;
    char nam[8];
    static int registry;

    //Evaluate processnoise
    virtual K1 processnoise( K1 v ) {
      if (v != 0.0) this->xVar = v;
      return this->xVar;
    }

    //Estimate the current state
    virtual K1 estimate() {
      return this->xPrv+this->xCon;
    }
  public:
    K1Filter (K1 processNoise, K1 sensorNoise=0.0, K1Filter* parent=NULL, K1 initState = 0.0, K1 initNoise = 0.0, char* nam=NULL) {
      this->dt = 0.0;
      this->child = NULL;  
      this->parent = parent;  
      this->xCur = initState;  
      this->vCur = initNoise;
      this->invar = sensorNoise;
      this->xCon = 0.0;
      this->input1 = initState;
      this->invar1 = initNoise;
      this->input0 = 0.0;
      this->invar0 = 0.0;
      this->sensors0 = 0;
      this->sensors1 = 1;
      this->gain = 0.0;
      this->intim0 = 0;
      this->intim1 = 0;
      this->dt = 0.001;
      if (nam)
        strncpy(this->nam, nam, 8);
      else
        snprintf(this->nam,8, "Kf%02i", registry++);
    }

    //Sensor measurement
    K1 measure( K1 value, unsigned char state = 0, unsigned long t=0, K1 variance = 0.0) {
      if (state>0) return(owner(state)->measure(value, 0, t, variance));
      this->intim0 += ((t?t:millis()) - this->intim0) / ++this->sensors0;
      variance = this->getSensVar(0, variance);
      if ( this->sensors0 == 1 ) {
        this->input0 = value;
        this->invar0 = variance;
      } else {  
        this->input0 += (value - this->input0)*this->invar0/(this->invar0 + variance);
        this->invar0 *= variance/(variance + this->invar0);
      }
      return(this->input0);
    }    

    //Sensor measurement relative to current state
    K1 difMeasure( K1 dvalue, unsigned char state = 0, unsigned long t=0, K1 variance = 0.0) {
      if (state>0) return(owner(state)->difMeasure(dvalue, 0, t, variance));
      variance = this->getSensVar(0, variance);
      dvalue += (this->xCur - dvalue)*variance/(variance + this->vCur);
      variance *= -variance/(variance + this->vCur);
      return(measure(dvalue, state, t, variance));
    }    

    //Process sensor measurement
    virtual K1 digest(K1 value, unsigned char state = 0, unsigned long t=0, K1 variance = 0.0) {
      return measure(value, state, t, variance>0.0?-variance:variance);
    }
    
    //Update the current state
    virtual void update( K1 procnoise = 0.0 ) {
      if (!this->sensors0) {
        this->dt = (millis() - this->intim1)/1000.0;
        procnoise=processnoise(procnoise);
        this->gain = 0.0;
        this->vCur += procnoise;
        this->vEst = this->vCur;
        this->xCur = this->xEst = estimate();
        return;
      }
      this->dt = (this->intim0 - this->intim1)/1000.0;
      procnoise=processnoise(procnoise);
      this->xPrv = this->xCur;
      this->vPrv = this->vCur;
      this->vEst = this->vPrv + procnoise;
      this->gain = this->vEst / (this->vEst + this->invar0);
      this->vCur = (1.0 - this->gain)*this->vEst;
      this->xEst = estimate();
      this->xCur = this->xEst + this->gain*(this->input0 - this->xEst);
      this->xCon = 0.0;
      this->input1 = this->input0;
      this->intim1 = this->intim0;
      this->sensors1 = this->sensors0;
      this->sensors0 = 0;
    }    

    //Get the linked filter that handles the state index given
    K1Filter* whois(unsigned char state=0) {
      K1Filter* f=this;
      while (f->parent) f=f->parent;
      while (state--) f=f->child;
      return(f);
    }    

    //Get the state index
    unsigned char whoami() {
      unsigned char i=0;
      K1Filter* f=this;
      while (f=f->parent) i++;
      return(i);
    }    

    //Get the number of states of linked filters
    unsigned char dimension() {
      unsigned char i=1;
      K1Filter* f=whois(0);
      while (f=f->child) i++;
      return(i);
    }    

    //Get the number of downwards linked filters 
    unsigned char scions() {
      unsigned char i=0;
      K1Filter* f=this;
      while (f=f->child) i++;
      return(i);
    }    

    //Get a downwards linked filter 
    inline K1Filter* owner(unsigned char state=0) {
      K1Filter* f=this;
      while (state--) f=f->child;
      return(f);
    }    

    //Get fused sensor reading 
    inline K1 getSensor(unsigned char state=0) {
      return(owner(state)->input0);
    }    

    //Get fused sensor variance 
    inline K1 getMeasurementVar(unsigned char state=0) {
      return(owner(state)->invar0);
    }    

    //Get fused sensor difference between consecutive updates 
    K1 getMeasurementDif(unsigned char state=0) {
      K1Filter* f=owner(state);
      if (f->intim0) return(f->input0 - f->input1);
      return(f->xCur - f->input1);
    }    

    //Get the variacne of fused sensor difference between consecutive updates 
    K1 getMeasurementDifVar(unsigned char state=0) {
      K1Filter* f=owner(state);
      if (f->intim0) return(f->invar0 + f->invar1);
      return(f->vCur + f->invar1);
    }    

    //Get fused sensor derivative between consecutive updates 
    K1 getMeasurementDrv(unsigned char state=0) {
      K1Filter* f=owner(state);
      return(1000.0*f->getMeasurementDif()/((f->intim0?f->intim0 : millis()) - f->intim1));
    }    

    //Get the default sensor variance (=0) or override it (<0) or set it (>0) for a given state
    K1 getSensVar(unsigned char state=0, K1 variance = 0.0) {
      K1Filter* f = owner(state);
      if (variance > 0 )
        f->invar= variance;
      else
        variance = variance == 0.0 ? f->invar: -variance;
      return(variance);
    }    

    //Set the default sensor variance for a given state
    inline void setSensVar(K1 variance, unsigned char state=0) {
      owner(state)->invar = variance;
    }    

    //Get the current sensor readings
    inline K1 getSensNum(unsigned char state=0) {
      return(owner(state)->sensors);
    }    

    //Get an array of linked filter member variables
    virtual void getArray(K1 *arr, unsigned char type=K1_STATE) {
      K1Filter* f=this;
      while (f) {
        switch (type) {
          case K1_STATE: {*arr = f->xCur; break;}
          case K1_STATE_PREVIOUS: {*arr = f->xPrv; break;}
          case K1_STATE_ESTIMATED: {*arr = f->xEst; break;}
          case K1_VARIANCE: {*arr = f->vCur; break;}
          case K1_VARIANCE_PREVIOUS: {*arr = f->vPrv; break;}
          case K1_VARIANCE_ESTIMATED: {*arr = f->vEst; break;}
          case K1_PROCESS_NOISE: {*arr = f->xVar; break;}
          case K1_CONTROL: {*arr = f->xCon; break;}
          case K1_PERIOD: {*arr = f->dt; break;}
          case K1_SENSORS: {*arr = f->input0; break;}
          case K1_SENSORS_PREVIOUS: {*arr = f->input1; break;}
          case K1_SENSOR_VARIANCE: {*arr = f->invar0; break;}
          case K1_SENSOR_VARIANCE_PREVIOUS: {*arr = f->invar1; break;}
          case K1_SENSORS_NUM: {*arr = f->sensors0; break;}
          case K1_SENSORS_NUM_PREVIOUS: {*arr = f->sensors1; break;}
        }
        arr++;
        f=f->child;
      }
    }    

    //Set linked filter member variables from an array
    virtual void setArray(K1 *arr, unsigned char type) {
      K1Filter* f=this;
      unsigned char d; 
      if (type == K1_SENSOR_VARIANCE) d=this->dimension();
      while (f) {
        switch (type) {
          case K1_STATE: {f->xCur = *arr; break;}
          case K1_VARIANCE: {f->vCur = *arr; break;}
          case K1_CONTROL: {f->xCon = *arr; break;}
          case K1_SENSORS: {f->measure(*arr); break;}
          case K1_SENSOR_VARIANCE: {if ( arr[d] > 0.0 ) f->measure(*arr, 0, 0, arr[d]); break;}
        }
        arr++;
        f=f->child;
      }
    }    

    //Process sensor readings from an array
    void digestArray(K1 *arr, unsigned char type) {
      unsigned char i,d=this->dimension(); 
      for (i=0; i<d; i++)
        if (type == K1_SENSOR_VARIANCE) {
          if (arr[i+d] == 0.0) continue;
          digest(arr[i], 0, 0, arr[i+d]);
        } else
          digest(arr[i]);
    }      

    //Get the current state
    inline K1 getState(unsigned char state=0) {
      return(owner(state)->xCur);
    }    

    //Set current state
    inline void setState(K1 value, unsigned char state=0) {
      owner(state)->xCur = value;
    }    

    //Get previous state
    inline K1 getStatePrev(unsigned char state=0) {
      return(owner(state)->xPrv);
    }    

    //Get current state estimate
    inline K1 getStateEst(unsigned char state=0) {
      return(owner(state)->xEst);
    }    

    //Get current state variance
    inline K1 getVariance(unsigned char state=0) {
      return(owner(state)->vCur);
    }    

    //Set current state variance
    inline void setVariance(K1 variance, unsigned char state=0) {
      owner(state)->vCur = variance;
    }    

    //Get current control input
    inline K1 getControl(unsigned char state=0) {
      return(owner(state)->xCon);
    }    

    //Set current control input
    inline void setControl(K1 controlInput, unsigned char state=0) {
      owner(state)->xCon = controlInput;
    }    

    //Get previous state variance
    inline K1 getVarPrev(unsigned char state=0) {
      return(owner(state)->vPrv);
    }    

    //Get estimated state variance
    inline K1 getVarEst(unsigned char state=0) {
      return(owner(state)->vEst);
    }    

    //Get process noise variance
    inline K1 getProcessNoise(unsigned char state=0) {
      return(owner(state)->xVar);
    }    

    //Set process noise variance
    inline void setProcessNoise(K1 noise, unsigned char state=0) {
      owner(state)->xVar = noise;
    }        

    //Set time elapsed between the two last updates
    inline K1 getPeriod(unsigned char state=0) {
      return(owner(state)->dt);
    }    

    //Link with a previous filter
    inline void setParent(K1Filter *parent, unsigned char state=0) {
      owner(state)->parent = parent;
    }        

    //Get upwards link to previous filter
    inline K1Filter* getParent(unsigned char state=0) {
      return(owner(state)->parent);
    }

    //Link with a subsequent filter
    inline void setChild(K1Filter *child, unsigned char state=0) {
      owner(state)->child = child;
    }        

    //Get downwards link to subsequent filter
    inline K1Filter* getChild(unsigned char state=0) {
      return(owner(state)->child);
    }

    //Set the name of the filter
    inline void setName(char* nam, unsigned char state=0) {
      strncpy(owner(state)->nam,8,nam);
    }

    //Set the name of the filter from a String type
    void setName(String nam, unsigned char state=0) {
      setName(const_cast<char*>(nam.c_str()), state);
    }

    //Get the name of the filter
    inline const char* getName(unsigned char state=0) {
      return(owner(state)->nam);
    }

  };

template <class K1> int K1Filter<K1>::registry = 0;

//A cluster of 1D Kalman filters for linear system equation x(t+dt) = x(t) + x'(t).dt + 1/2!.x"(t)dt^2 + ...
template <class K1>
class K1FilterLinear: public K1Filter<K1> {
  friend class K1FilterMobile<K1>;
  protected:
    unsigned char order;
    K1FilterLinear* master;
    K1 period;

    virtual K1 processnoise( K1 v ) {
      K1FilterLinear* f;
      if (v != 0.0 && this == master) {
        f = this->owner(this->order-1);
        f->xVar = v;
        while (f!=this) {
          f = f->parent;
          f->xVar = f->period / (f->order - 1.0);
          f->xVar *= f->xVar * f->child->xVar; 
        }
      }
      v = this->dt / this->period;
      return this->xVar*v*v;
    }

    virtual K1 estimate() {
      int i=0;
      K1FilterLinear *f=this;
      K1 est = this->xPrv;
      K1 coef = 1.0;
      while (++i < this->order) {
        f = f->child;
        coef *= this->period / i;
        est += coef * f->xCur;
      }
      if ( this->order > 1 ) {        
        if ( this->sensors0 )
          this->child->measure((this->input0 - this->input1)/this->dt, 0, this->intim0, (this->invar0 + this->invar1)/this->dt/this->dt);
        this->child->update();
      }
      est+=this->xCon;
      return est;
    }
  public:
    K1FilterLinear (unsigned char order, K1 periodMS, K1 processNoise, K1Filter<K1> *parent=NULL, char* nam=NULL):K1Filter<K1>(0.0, 0.0, parent, 0.0, 0.0, nam) {
      K1 ms=abs(periodMS);
      char buf[8], *c, i=1;
      this->order=order;
      this->parent=parent; 
      this->period=ms<1?ms:ms/1000.0;
      if (periodMS>0.0)
        this->master = this;
      else
        this->master = static_cast<K1FilterLinear*>(this->parent)->master;  
      if (order>1) {
        strcpy(buf,this->nam);
        c = strchr(buf,'.');
        if ( c )
          i = atoi(c+1) + 1;
        else
          c = strchr(buf,'\0');
        sprintf(c,".%i",i);
        this->child=new K1FilterLinear<K1>(order-1, -ms, processNoise, this, buf);
        this->xVar = this->period / (order - 1.0);
        this->xVar *= this->xVar * this->child->xVar;
      } else {
        this->xVar = processNoise;
        this->child = NULL;
      }
    }
  };

//A cluster of 1D Kalman filters for holonomic mobile robot
template <class K1>
class K1FilterMobile: public K1FilterLinear<K1> {
  protected:
    K1FilterLinear<K1> *ang;
    K1Filter<K1> *xc;
    K1Filter<K1> *yc;
    K1 d, uvar;
    virtual K1 processnoise( K1 v ) {
      K1 n;
      if (v != 0.0 ) {
        n = K1FilterLinear<K1>::processnoise( v );
        this->uvar = 2.0 * v / this->d / this->d;
      } else {
        n = this->dt / this->period;
        n *= n * this->xVar;
      }  
      xc->processnoise( n );
      yc->processnoise( n );
      return n;
    }

    virtual K1 estimate() {
      K1 est = K1FilterLinear<K1>::estimate();
      K1 ds, dx, dy, vs, vh, v, h, ch, sh;
      unsigned long t;
      if ( this->sensors0 && this->ang->sensors0 ) {
        ds=getMeasurementDif();
        h = ang->xCur + ang->getMeasurementDif() / 2.0;
        h = trunc(h/TWO_PI)*TWO_PI;
        sh = sin(h);
        ch = cos(h);
        dx = ds*ch;
        dy = ds*sh;
        vs = getMeasurementDifVar(); //V(ds)
        vh=ang->vCur + ang->getMeasurementDifVar()/2.0; //V(h+dh)
        v = vh*(1.0 - (1 - vh/4.0)*ch*ch); //Vcos(h+dh)
        t = (this->intim0+ang->intim0)>>1;          
        this->xc->difMeasure(dx,0, t, vs*v+vs*ds*ds/4.0+v*h*h);
        v = vh*(1.0 - (1 - vh/4.0)*sh*sh); //Vsin(h+dh)          
        this->yc->difMeasure(dy,0, t, vs*v+vs*ds*ds/4.0+v*h*h);
      }
      this->ang->update(uvar);
      this->uvar = 0;
      this->xc->update();
      this->yc->update();
      return est;
    }
  public:
    K1FilterMobile (K1 wheelbase, K1 periodMS, K1 processNoise, K1Filter<K1> *parent=NULL, char* nam=NULL):K1FilterLinear<K1>(3, periodMS, processNoise, parent, nam) {
      char buf[8];
      this->uvar = 0.0;
      this->d = wheelbase;
      wheelbase *= wheelbase;
      snprintf(buf,8,"%s.3",this->nam);
      this->ang = new K1FilterLinear<K1>(3, periodMS, 2.0*processNoise/wheelbase, owner(K1STATE_ACCELERATION), buf);
      owner(K1STATE_ACCELERATION)->child = this->ang;
      periodMS *= periodMS;
      periodMS *= periodMS;
      snprintf(buf,8,"%s.6",this->nam);
      this->xc  = new K1Filter<K1>(processNoise*periodMS, 0.0, owner(K1STATE_ANGACC), 0.0, 0.0, buf);      
      snprintf(buf,8,"%s.7",this->nam);
      this->yc  = new K1Filter<K1>(processNoise*periodMS, 0.0, owner(K1STATE_ANGACC), 0.0, 0.0, buf);
      owner(K1STATE_ANGACC)->child = this->xc;
      this->xc->child = this->yc;
    }

    virtual K1 digest(K1 value, unsigned char state = 0, unsigned long t=0, K1 variance = 0.0) {
      K1Filter<K1> *vf, *wf, *af;
      K1 e, v, r;
      if (state != K1STATE_ANGACC ) return measure(value, state, t, variance);
      vf = owner(K1STATE_VELOCITY);
      wf = owner(K1STATE_ANGVEL);
      af = wf->child;
      variance = af->getSensVar(0, variance);
      if (vf->xCur == 0.0) return(0.0);
      r = wf->xCur / vf->xCur; //1/r
      if (abs(r) > 100.0 || abs(r) < 0.01) return(0.0);
      v = 1 / vf->xCur;
      v *= v;
      e = (1.0 + vf->vCur*v) / vf->xCur; //E(1/w)
      v *= vf->vCur*(1+v*vf->vCur); // V(1/w)
      v =  wf->vCur*(v + wf->xCur*wf->xCur) + v*e*e; //V(w/v)
      e *= wf->xCur; //E(w/v);
      v = af->vCur*(variance + af->xCur*af->xCur) + v*e*e;//V(linear accel/turn radius)
      return af->measure(value*r, 0, t, -v);
    }
  };

#endif
