/*
 * This file is part of the Kalman1D library
 * Usage: Provide an example use of the library
 * 
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
 */
#include "Kalman1D.h"
#include <limits.h>

#define NOISEDEV_ACCEL 0.009
#define NOISEDEV_GYRO 0.52
#define NOISEDEV_MAGN 0.05
#define NOISEDEV_GPS 1.7
#define PERIOD_MS 50
#define STATE_DIM 8

typedef union {
  struct {
    double s; //distance (m)
    double v; //speed (m/s)
    double g; //acceleration (m/s^2)
    double h; //heading (rad)
    double w; //angular velocity (rad/s)
    double a; //angular acceleration (rad/s^2)
    double x; //displacement X axis (m)
    double y; //displacement Y axis (m)
  };
  double state[STATE_DIM];
} State;

K1FilterMobile<double> *filter;

unsigned int wait() {
  static unsigned int load = 0;
  static unsigned long loopMS=PERIOD_MS;
  unsigned long now=millis();
  if (now>loopMS) {
    load = 100<<8;
    loopMS = now + PERIOD_MS;
  } else {
    load += (25600 - (loopMS-now)*(25600/PERIOD_MS) - load)>>3;
    while (loopMS>millis());
    loopMS += PERIOD_MS;
  }
  return(load);
}

const double r = 0.5; //(m) radius of circular movement
const double w = 0.5; //(rad/s) angular speed
const uint16_t cpm = 20000; //clicks per meter to simulate quadrature encoder
const double d = 0.1; //(m) wheelbase
const double d2 = d/2.0; //(m) wheelbase

void setup() {
  Serial.begin(115200);
  randomSeed(analogRead(0)); //assuming A0 is not connected
  filter=new K1FilterMobile<double>(d, PERIOD_MS, sq(w*r*(PERIOD_MS/1000.0)));
}

void loop() {
  static double t=0;
  static State mean = {0.0};
  static State variance = {0.0};
  static unsigned long sample = 1;
  State actual, recover, error, accumulate;
  State measure[2] = {0.0};
  double encLft, encRgt, compass, err, test;
  int i, n, c, load;
  t += PERIOD_MS/1000.0; //update simulation time
  actual.v = w*t;
  actual.x = sin(actual.v);
  actual.y = w*r;
  actual.h = actual.v - actual.x; //(rads) heading (also the angle along the circle) h=ωt-sinωt
  actual.s = r*actual.h; //the actual distance travelled s=ωr(t-sinωt/ω)
  actual.g = actual.y*actual.x; //the actual acceleration along the circle γ=ωr.sinωt
  actual.a = w*actual.x; //(rads/s^2) the actual angular acceleration
  actual.v = actual.y*(1.0-cos(actual.v)); //the actual velocity along the circle 
  actual.w = actual.v/r; //(rads/s) the actual angular velocity
  actual.x = r*sin(actual.h); //the actual x position
  actual.y = r*(cos(actual.h)-1.0); //the actual y position
  load = wait()>>8; //Sample at exact intervals
  encLft = trunc(actual.s*(1.0+d2/r)*cpm)/cpm; //Left encoder (only quantization noise)
  encRgt = trunc(actual.s*(1.0-d2/r)*cpm)/cpm; //Right encoder (only quantization noise)
  measure[0].s = (encLft + encRgt) / 2.0; 
  measure[1].s = sq(1.0/cpm)/12.0;
  measure[0].g = actual.g + NOISEDEV_ACCEL*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX; 
  measure[1].g = sq(NOISEDEV_ACCEL);
  measure[0].h = (encLft - encRgt) / d; 
  measure[1].h = measure[1].s / sq(d2);
  measure[0].w = actual.w + NOISEDEV_GYRO*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX;
  measure[1].w = sq(NOISEDEV_GYRO);
#ifdef NOISEDEV_GPS
  measure[0].x = actual.x + NOISEDEV_GPS*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX; 
  measure[1].x = sq(NOISEDEV_GPS);
  measure[0].y = actual.y + NOISEDEV_GPS*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX; 
  measure[1].y = sq(NOISEDEV_GPS);
#endif  
  filter->setArray(measure[0].state,K1_SENSOR_VARIANCE); //update multiple sensors

  //simulate compass
  compass = actual.h + NOISEDEV_MAGN*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX;
  if (compass>0)
    compass -= TWO_PI*trunc(compass/TWO_PI);
  else
    compass += TWO_PI*(1+trunc(compass/TWO_PI));

  //recover rotation angle from compass
  n=trunc(recover.h/TWO_PI);
  if (recover.h<0) n--;
  test = 100.0;
  for (i=-1; i<2; i++) {
    err = abs((n-i)*TWO_PI+compass - recover.h);
    if (err > test) break;
    c = i;
    test = err;
  }
  filter->measure((n-c)*TWO_PI+compass, K1STATE_HEADING, 0,sq(NOISEDEV_MAGN));
    
  filter->digest(actual.a + NOISEDEV_ACCEL*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX, K1STATE_ANGACC, 0, sq(NOISEDEV_ACCEL)); //compute angular acceleration from linear  

  filter->update(); // apply filter
  filter->getArray(recover.state);//get current state vector(array)

  // calculate mean and variance of error for all states
  for (i=0; i<STATE_DIM; i++) {
    error.state[i] = recover.state[i] - actual.state[i] - mean.state[i]; 
    accumulate.state[i] = error.state[i] / sample;
    if (sample>1) variance.state[i] -= variance.state[i] / ( sample - 1);
    variance.state[i] += error.state[i] * accumulate.state[i];
    mean.state[i] += accumulate.state[i];
  }
  
  if (sample++ % 100 == 0) {
    Serial.print("Stdev s:");
    Serial.print(sqrt(variance.s),4);
    Serial.print(" v:");
    Serial.print(sqrt(variance.v),4);
    Serial.print(" g:");
    Serial.println(sqrt(variance.g),4);
    Serial.print("...   h:");
    Serial.print(sqrt(variance.h),4);
    Serial.print(" w:");
    Serial.print(sqrt(variance.w),4);
    Serial.print(" a:");
    Serial.println(sqrt(variance.a),4);
    Serial.print("...   x:");
    Serial.print(sqrt(variance.x),4);
    Serial.print(" y:");
    Serial.print(sqrt(variance.y),4);
    Serial.print(" cur_pos_err:");
    Serial.print(sqrt(sq(actual.x-recover.x)+sq(actual.y-recover.y)),4);
    Serial.print(" load:");
    Serial.print(load);
    Serial.println("%");
  }
}
