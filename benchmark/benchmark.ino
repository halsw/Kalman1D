/*
 * This file is part of the Kalman1D library
 * Usage: Benchmatk the iteration time of the filter
 * 
 * Version 1.0
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
#define BENCHMARK 1000
const double d = 0.1; //(m) wheelbase
const double dt = 0.003; //(s) update period
const double np = 100.0; //process noise
const double nm = 0.0001; //measurement noise
const double ds = sq(dt)*sqrt(2.0*np); //step displacement for test
int counter;
unsigned long bmtime;
K1FilterMobile <double> filter(d,dt,np);
double s=0;

void setup()  { 
  Serial.begin(115200);
  while (! Serial);
  Serial.println("1D Kalman Filter benchmark");
  filter.setSensVar(K1STATE_DISPLACEMENT, nm);
  filter.setSensVar(K1STATE_HEADING, nm);
  counter = 0;
  bmtime = millis();
}

void loop() {
  if (counter++ < BENCHMARK) {
    s += ds;
    filter.measure(s, K1STATE_DISPLACEMENT);
    filter.measure(s, K1STATE_HEADING);
    filter.update();
  } else {
    float t=(millis()-bmtime) / (BENCHMARK*1.0);
    Serial.print("Filter update in ");
    Serial.print(t);
    Serial.println(" ms");
    counter = 0;
    bmtime = millis();
  }
}
