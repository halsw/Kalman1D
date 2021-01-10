# [One Dimension Kalman Filter](https://github.com/halsw/Kalman1D) for Arduino and Teensy
A library with template classes for the 1D Kalman Filter with asynchronous sensor fusion

The 1D Kalman Filter is, as the name implies, for states described by a scalar, although this library caters for arrays of 1D Kalman Filters that can be used on linear and non-linear systems.

These do not perform as well as the Kalman Filter because they don't account for the cross covariance between state elements but the simplification also means a performance gain in terms of speed. This way a state containing numerous elements can be updated in few msecs.

## Process Noise
Process noise is the [variance](https://en.wikipedia.org/wiki/Variance) of the filter input between iterations which must be defined either in the constructor or in each iteration by calling the method **update()**

## Measurement Noise
Measurement noise is the [variance](https://en.wikipedia.org/wiki/Variance) of the error between the measured and actual values which must be defined in the filter either in the constructor or in each sensor measurement by calling one of the methods **measure()**, **difMeasure()** or **digest()**

## Example
The example .ino file simulates a holonomic 2-wheel robot that runs around in a circle with sinusoidal acceleration profile and has a 8D state containing its linear distance travelled, velocity and acceleration, its angular displacement travelled, velocity and acceleration, as well as its x and y coordinates
