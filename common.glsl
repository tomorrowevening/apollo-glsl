/**
 * - Random / Noise
 * - Rotate2D / Scale2D
 * - RGB to HSB / HSB to RGB
 * - 2D Shapes
 */

#ifndef TWO_PI
#define TWO_PI 6.28318530718
#endif

//////////////////////////////////////////////////
// Common

/**
 * Actual pixel value
 */

vec2 getPixel(vec2 uv, vec2 res) {
  return vec2(uv.x * res.x, res.y - (uv.y * res.y));
}

//////////////////////////////////////////////////
// Random / Noise

float random(vec2 st) {
  return fract(sin(dot(st.xy,
                       vec2(12.9898,78.233)))*
      43758.5453123);
}

// Based on thebookofshaders
// https://github.com/patriciogonzalezvivo/thebookofshaders/blob/master/18/grain.frag
vec3 random3(vec3 c) {
  float j = 4096.0*sin(dot(c,vec3(17.0, 59.4, 15.0)));
  vec3 r;
  r.z = fract(512.0*j);
  j *= .125;
  r.x = fract(512.0*j);
  j *= .125;
  r.y = fract(512.0*j);
  return r-0.5;
}

// 2D Noise based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise(in vec2 st) {
  vec2 i = floor(st);
  vec2 f = fract(st);

  // Four corners in 2D of a tile
  float a = random(i);
  float b = random(i + vec2(1.0, 0.0));
  float c = random(i + vec2(0.0, 1.0));
  float d = random(i + vec2(1.0, 1.0));

  // Smooth Interpolation

  // Cubic Hermine Curve.  Same as SmoothStep()
  vec2 u = f*f*(3.0-2.0*f);
  // u = smoothstep(0.,1.,f);

  // Mix 4 coorners porcentages
  return mix(a, b, u.x) +
          (c - a)* u.y * (1.0 - u.x) +
          (d - b) * u.x * u.y;
}

const float F3 =  0.3333333;
const float G3 =  0.1666667;
float snoise(vec3 p) {
  vec3 s = floor(p + dot(p, vec3(F3)));
  vec3 x = p - s + dot(s, vec3(G3));
  vec3 e = step(vec3(0.0), x - x.yzx);
  vec3 i1 = e*(1.0 - e.zxy);
  vec3 i2 = 1.0 - e.zxy*(1.0 - e);
  vec3 x1 = x - i1 + G3;
  vec3 x2 = x - i2 + 2.0*G3;
  vec3 x3 = x - 1.0 + 3.0*G3;
  vec4 w, d;
  w.x = dot(x, x);
  w.y = dot(x1, x1);
  w.z = dot(x2, x2);
  w.w = dot(x3, x3);
  w = max(0.6 - w, 0.0);
  d.x = dot(random3(s), x);
  d.y = dot(random3(s + i1), x1);
  d.z = dot(random3(s + i2), x2);
  d.w = dot(random3(s + 1.0), x3);
  w *= w;
  w *= w;
  d *= w;
  return dot(d, vec4(52.0));
}

//  Author: Matt DesLauriers
//  https://github.com/mattdesl/glsl-film-grain/blob/master/index.glsl
//
float grain(vec2 texCoord, vec2 resolution, float frame, float multiplier) {
  vec2 mult = texCoord * resolution;
  float offset = snoise(vec3(mult / multiplier, frame));
  float n1 = snoise(vec3(mult, offset));
  return n1 / 2.0 + 0.5;
}

float grain(vec2 texCoord, vec2 resolution, float frame) {
  return grain(texCoord, resolution, frame, 2.5);
}

float grain(vec2 texCoord, vec2 resolution) {
  return grain(texCoord, resolution, 0.0);
}

//////////////////////////////////////////////////
// Matrices

mat2 rotate2D(float _angle){
  return mat2(cos(_angle),-sin(_angle),
              sin(_angle),cos(_angle));
}

mat2 scale2D(vec2 _scale){
  return mat2(_scale.x,0.0,
              0.0,_scale.y);
}

//////////////////////////////////////////////////
// Colors

vec3 rgb2hsb( in vec3 c ) {
  vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
  vec4 p = mix(vec4(c.bg, K.wz),
               vec4(c.gb, K.xy),
               step(c.b, c.g));
  vec4 q = mix(vec4(p.xyw, c.r),
               vec4(c.r, p.yzx),
               step(p.x, c.r));
  float d = q.x - min(q.w, q.y);
  float e = 1.0e-10;
  return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)),
              d / (q.x + e),
              q.x);
}

//  Function from Iñigo Quiles
//  https://www.shadertoy.com/view/MsS3Wc
vec3 hsb2rgb( in vec3 c ) {
  vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                           6.0)-3.0)-1.0,
                   0.0,
                   1.0 );
  rgb = rgb*rgb*(3.0-2.0*rgb);
  return c.z * mix(vec3(1.0), rgb, c.y);
}

//////////////////////////////////////////////////
// Shapes

/**
 * Smooths your float, a cheap anti-aliasing, you'll need this for most of these shape functions
 */
float smoothedge(float value, float smoothAmount, vec2 res) {
  return smoothstep(0.0, smoothAmount / res.x, value);
}

float line(vec2 _st, vec2 _p1, vec2 _p2, float _width, float _spread) {
  _width = 1.0 / _width;
  vec2 p2p1 = _p1 - _p2;
  vec2 p1p2 = -(p2p1);
  vec2 p2p = _st - _p2;
  vec2 p1p = _st - _p1;
  vec2 pp1 = -(p1p);
  vec2 pd = normalize(vec2(p2p1.y, -p2p1.x));
  float proj = dot(pd, pp1);
  float pr1 = dot(p2p1, p2p);
  float pr2 = dot(p1p2, p1p);

  if(pr1 > 0.0 && pr2 > 0.0) {
      return pow(1.0 / abs(proj * _width), _spread);
  } else {
      return 0.0;
  }
}

float circle(vec2 p, float radius) {
  return length(p) - radius;
}

float rect(vec2 p, vec2 size) {
  vec2 d = abs(p) - size;
  return min(max(d.x, d.y), 0.0) + length(max(d,0.0));
}

float roundRect(vec2 p, vec2 size, float radius) {
  vec2 d = abs(p) - size;
  return min(max(d.x, d.y), 0.0) + length(max(d,0.0))- radius;
}

float ring(vec2 p, float radius, float width) {
  return abs(length(p) - radius * 0.5) - width;
}

float hexagon(vec2 p, float radius) {
  vec2 q = abs(p);
  return max(abs(q.y), q.x * 0.866025 + q.y * 0.5) - radius;
}

float triangle(vec2 p, float size) {
  vec2 q = abs(p);
  return max(q.x * 0.866025 + p.y * 0.5, -p.y * 0.5) - size * 0.5;
}

float ellipse(vec2 p, vec2 r, float s) {
  return (length(p / r) - s);
}

float capsule(vec2 p, vec2 a, vec2 b, float r) {
  vec2 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

//http://thndl.com/square-shaped-shaders.html
float polygon(vec2 p, int vertices, float size) {
  float a = atan(p.x, p.y);
  float b = 6.28319 / float(vertices);
  return cos(floor(0.5 + a / b) * b - a) * length(p) - size;
}

//////////////////////////////////////////////////
// Export

#pragma glslify: export(getPixel)

#pragma glslify: export(random)
#pragma glslify: export(random3)
#pragma glslify: export(noise)
#pragma glslify: export(snoise)
#pragma glslify: export(grain)

#pragma glslify: export(rotate2D)
#pragma glslify: export(scale2D)

#pragma glslify: export(rgb2hsb)
#pragma glslify: export(hsb2rgb)

#pragma glslify: export(line)
#pragma glslify: export(circle)
#pragma glslify: export(rect)
#pragma glslify: export(roundRect)
#pragma glslify: export(ring)
#pragma glslify: export(hexagon)
#pragma glslify: export(triangle)
#pragma glslify: export(ellipse)
#pragma glslify: export(capsule)
#pragma glslify: export(polygon)
