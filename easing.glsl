/**
 * - Penner's Equations
 * - Cubic Bezier
 */

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef HALF_PI
#define HALF_PI 1.5707963267948966
#endif

//////////////////////////////////////////////////
// Easings

float backIn(float t) {
  return pow(t, 3.0) - t * sin(t * PI);
}

float backOut(float t) {
  float f = 1.0 - t;
  return 1.0 - (pow(f, 3.0) - f * sin(f * PI));
}

float backInOut(float t) {
  float f = t < 0.5
    ? 2.0 * t
    : 1.0 - (2.0 * t - 1.0);

  float g = pow(f, 3.0) - f * sin(f * PI);

  return t < 0.5
    ? 0.5 * g
    : 0.5 * (1.0 - g) + 0.5;
}

float bounceOut(float t) {
  const float a = 4.0 / 11.0;
  const float b = 8.0 / 11.0;
  const float c = 9.0 / 10.0;

  const float ca = 4356.0 / 361.0;
  const float cb = 35442.0 / 1805.0;
  const float cc = 16061.0 / 1805.0;

  float t2 = t * t;

  return t < a
    ? 7.5625 * t2
    : t < b
      ? 9.075 * t2 - 9.9 * t + 3.4
      : t < c
        ? ca * t2 - cb * t + cc
        : 10.8 * t * t - 20.52 * t + 10.72;
}

float bounceIn(float t) {
  return 1.0 - bounceOut(1.0 - t);
}

float bounceInOut(float t) {
  return t < 0.5
    ? 0.5 * (1.0 - bounceOut(1.0 - t * 2.0))
    : 0.5 * bounceOut(t * 2.0 - 1.0) + 0.5;
}

float circIn(float t) {
  return 1.0 - sqrt(1.0 - t * t);
}

float circOut(float t) {
  return sqrt((2.0 - t) * t);
}

float circInOut(float t) {
  return t < 0.5
    ? 0.5 * (1.0 - sqrt(1.0 - 4.0 * t * t))
    : 0.5 * (sqrt((3.0 - 2.0 * t) * (2.0 * t - 1.0)) + 1.0);
}

float cubicIn(float t) {
  return t * t * t;
}

float cubicOut(float t) {
  float f = t - 1.0;
  return f * f * f + 1.0;
}

float cubicInOut(float t) {
  return t < 0.5
    ? 4.0 * t * t * t
    : 0.5 * pow(2.0 * t - 2.0, 3.0) + 1.0;
}

float elasticIn(float t) {
  return sin(13.0 * t * HALF_PI) * pow(2.0, 10.0 * (t - 1.0));
}

float elasticOut(float t) {
  return sin(-13.0 * (t + 1.0) * HALF_PI) * pow(2.0, -10.0 * t) + 1.0;
}

float elasticInOut(float t) {
  return t < 0.5
    ? 0.5 * sin(+13.0 * HALF_PI * 2.0 * t) * pow(2.0, 10.0 * (2.0 * t - 1.0))
    : 0.5 * sin(-13.0 * HALF_PI * ((2.0 * t - 1.0) + 1.0)) * pow(2.0, -10.0 * (2.0 * t - 1.0)) + 1.0;
}

float expoIn(float t) {
  return t == 0.0 ? t : pow(2.0, 10.0 * (t - 1.0));
}

float expoOut(float t) {
  return t == 1.0 ? t : 1.0 - pow(2.0, -10.0 * t);
}

float expoInOut(float t) {
  return t == 0.0 || t == 1.0
    ? t
    : t < 0.5
      ? +0.5 * pow(2.0, (20.0 * t) - 10.0)
      : -0.5 * pow(2.0, 10.0 - (t * 20.0)) + 1.0;
}

float linear(float t) {
  return t;
}

float quadIn(float t) {
  return t * t;
}

float quadOut(float t) {
  return -t * (t - 2.0);
}

float quadInOut(float t) {
  float p = 2.0 * t * t;
  return t < 0.5 ? p : -p + (4.0 * t) - 1.0;
}

float quartIn(float t) {
  return pow(t, 4.0);
}

float quartOut(float t) {
  return pow(t - 1.0, 3.0) * (1.0 - t) + 1.0;
}

float quartInOut(float t) {
  return t < 0.5
    ? +8.0 * pow(t, 4.0)
    : -8.0 * pow(t - 1.0, 4.0) + 1.0;
}

float quinticIn(float t) {
  return pow(t, 5.0);
}

float quinticOut(float t) {
  return 1.0 - (pow(t - 1.0, 5.0));
}

float quinticInOut(float t) {
  return t < 0.5
    ? +16.0 * pow(t, 5.0)
    : -0.5 * pow(2.0 * t - 2.0, 5.0) + 1.0;
}

float sineIn(float t) {
  return sin((t - 1.0) * HALF_PI) + 1.0;
}

float sineOut(float t) {
  return sin(t * HALF_PI);
}

float sineInOut(float t) {
  return -0.5 * (cos(PI * t) - 1.0);
}

//////////////////////////////////////////////////
// Bezier

float slopeFromT(float t, float A, float B, float C) {
  float dtdx = 1.0/(3.0*A*t*t + 2.0*B*t + C);
  return dtdx;
}

float xFromT(float t, float A, float B, float C, float D) {
  float x = A*(t*t*t) + B*(t*t) + C*t + D;
  return x;
}

float yFromT(float t, float E, float F, float G, float H) {
  float y = E*(t*t*t) + F*(t*t) + G*t + H;
  return y;
}
float lineSegment(vec2 p, vec2 a, vec2 b) {
  vec2 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return smoothstep(0.0, 1.0 / u_resolution.x, length(pa - ba*h));
}

float cubicBezier(float x, vec2 a, vec2 b) {
  float y0a = 0.0; // initial y
  float x0a = 0.0; // initial x
  float y1a = a.y;    // 1st influence y
  float x1a = a.x;    // 1st influence x
  float y2a = b.y;    // 2nd influence y
  float x2a = b.x;    // 2nd influence x
  float y3a = 1.0; // final y
  float x3a = 1.0; // final x

  float A =   x3a - 3.0*x2a + 3.0*x1a - x0a;
  float B = 3.0*x2a - 6.0*x1a + 3.0*x0a;
  float C = 3.0*x1a - 3.0*x0a;
  float D =   x0a;

  float E =   y3a - 3.0*y2a + 3.0*y1a - y0a;
  float F = 3.0*y2a - 6.0*y1a + 3.0*y0a;
  float G = 3.0*y1a - 3.0*y0a;
  float H =   y0a;

  // Solve for t given x (using Newton-Raphelson), then solve for y given t.
  // Assume for the first guess that t = x.
  float currentt = x;
  for (int i=0; i < 5; i++){
    float currentx = xFromT (currentt, A,B,C,D);
    float currentslope = slopeFromT (currentt, A,B,C);
    currentt -= (currentx - x)*(currentslope);
    currentt = clamp(currentt,0.0,1.0);
  }

  float y = yFromT (currentt,  E,F,G,H);
  return y;
}

//////////////////////////////////////////////////
// Export

#pragma glslify: export(backIn)
#pragma glslify: export(backOut)
#pragma glslify: export(backInOut)
#pragma glslify: export(bounceIn)
#pragma glslify: export(bounceOut)
#pragma glslify: export(bounceInOut)
#pragma glslify: export(circIn)
#pragma glslify: export(circOut)
#pragma glslify: export(circInOut)
#pragma glslify: export(cubicIn)
#pragma glslify: export(cubicOut)
#pragma glslify: export(cubicInOut)
#pragma glslify: export(elasticIn)
#pragma glslify: export(elasticOut)
#pragma glslify: export(elasticInOut)
#pragma glslify: export(expoIn)
#pragma glslify: export(expoOut)
#pragma glslify: export(expoInOut)
#pragma glslify: export(linear)
#pragma glslify: export(quadIn)
#pragma glslify: export(quadOut)
#pragma glslify: export(quadInOut)
#pragma glslify: export(quartIn)
#pragma glslify: export(quartOut)
#pragma glslify: export(quartInOut)
#pragma glslify: export(quinticIn)
#pragma glslify: export(quinticOut)
#pragma glslify: export(quinticInOut)
#pragma glslify: export(sineIn)
#pragma glslify: export(sineOut)
#pragma glslify: export(sineInOut)
#pragma glslify: export(cubicBezier)
