/**
 * - Palette mixing
 * - Raymarching Distance functions
 * - Additional easing functions
 * - Sphere Raymarching
 */

//////////////////////////////////////////////////
// http://www.iquilezles.org/www/articles/palettes/palettes.htm

/**
 * cosine based palette, 4 vec3 params
 */
vec3 palette( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
  return a + b*cos( 6.28318*(c*t+d) );
}

//////////////////////////////////////////////////
// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

/**
 * Sphere - signed - exact
 */
float sdSphere( vec3 p, float s ) {
  return length(p)-s;
}

/**
 * Box - unsigned - exact
 */
float udBox( vec3 p, vec3 b ) {
  return length(max(abs(p)-b,0.0));
}

/**
 * Round Box - unsigned - exact
 */
float udRoundBox( vec3 p, vec3 b, float r ) {
  return length(max(abs(p)-b,0.0))-r;
}

/**
 * Box - signed - exact
 */
float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

/**
 * Torus - signed - exact
 */
float sdTorus( vec3 p, vec2 t ) {
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

/**
 * Cylinder - signed - exact
 */
float sdCylinder( vec3 p, vec3 c ) {
  return length(p.xz-c.xy)-c.z;
}

/**
 * Cone - signed - exact
 */
float sdCone( vec3 p, vec2 c ) {
  // c must be normalized
  float q = length(p.xy);
  return dot(c,vec2(q,p.z));
}

/**
 * Plane - signed - exact
 */
float sdPlane( vec3 p, vec4 n ) {
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

/**
 * Hexagonal Prism - signed - exact
 */
float sdHexPrism( vec3 p, vec2 h ) {
  vec3 q = abs(p);
  return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

/**
 * Triangular Prism - signed - exact
 */
float sdTriPrism( vec3 p, vec2 h ) {
  vec3 q = abs(p);
  return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}

/**
 * Capsule / Line - signed - exact
 */
float sdCapsule( vec3 p, vec3 a, vec3 b, float r ) {
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

/**
 * Capped cylinder - signed - exact
 */
float sdCappedCylinder( vec3 p, vec2 h ) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

/**
 * Capped Cone - signed - bound
 */
float sdCappedCone( in vec3 p, in vec3 c ) {
  vec2 q = vec2( length(p.xz), p.y );
  vec2 v = vec2( c.z*c.y/c.x, -c.z );
  vec2 w = v - q;
  vec2 vv = vec2( dot(v,v), v.x*v.x );
  vec2 qv = vec2( dot(v,w), v.x*w.x );
  vec2 d = max(qv,0.0)*qv/vv;
  return sqrt( dot(w,w) - max(d.x,d.y) ) * sign(max(q.y*v.x-q.x*v.y,w.y));
}

/**
 * Ellipsoid - signed - bound
 */
float sdEllipsoid( in vec3 p, in vec3 r ) {
  return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

float dot2( in vec3 v ) { return dot(v,v); }

/**
 * Triangle - unsigned - exact
 */
float udTriangle( vec3 p, vec3 a, vec3 b, vec3 c ) {
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 ac = a - c; vec3 pc = p - c;
  vec3 nor = cross( ba, ac );

  return sqrt(
  (sign(dot(cross(ba,nor),pa)) +
   sign(dot(cross(cb,nor),pb)) +
   sign(dot(cross(ac,nor),pc))<2.0)
   ?
   min( min(
   dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
   dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
   dot2(ac*clamp(dot(ac,pc)/dot2(ac),0.0,1.0)-pc) )
   :
   dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

/**
 * Quad - unsigned - exact
 */
float udQuad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d ) {
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 dc = d - c; vec3 pc = p - c;
  vec3 ad = a - d; vec3 pd = p - d;
  vec3 nor = cross( ba, ad );

  return sqrt(
  (sign(dot(cross(ba,nor),pa)) +
   sign(dot(cross(cb,nor),pb)) +
   sign(dot(cross(dc,nor),pc)) +
   sign(dot(cross(ad,nor),pd))<3.0)
   ?
   min( min( min(
   dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
   dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
   dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
   dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
   :
   dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

/**
 * Torus82 - signed
 */
float sdTorus82( vec3 p, vec2 t ) {
  vec2 q = vec2(length2(p.xz)-t.x,p.y);
  return length8(q)-t.y;
}

/**
 * Torus88 - signed
 */
float sdTorus88( vec3 p, vec2 t ) {
  vec2 q = vec2(length8(p.xz)-t.x,p.y);
  return length8(q)-t.y;
}

/**
 * Union
 */
float opU( float d1, float d2 ) {
  return min(d1,d2);
}

/**
 * Substraction
 */
float opS( float d1, float d2 ) {
  return max(-d1,d2);
}

/**
 * Intersection
 */
float opI( float d1, float d2 ) {
  return max(d1,d2);
}

/**
 * Repetition
 */
float opRep( vec3 p, vec3 c ) {
  vec3 q = mod(p,c)-0.5*c;
  return primitve( q );
}

/**
 * Rotation/Translation
 */
vec3 opTx( vec3 p, mat4 m ) {
  vec3 q = invert(m)*p;
  return primitive(q);
}

/**
 * Scale
 */
float opScale( vec3 p, float s ) {
  return primitive(p/s)*s;
}

/**
 * Displacement
 */
float opDisplace( vec3 p ) {
  float d1 = primitive(p);
  float d2 = displacement(p);
  return d1+d2;
}

/**
 * Blend
 */
float opBlend( vec3 p ) {
  float d1 = primitiveA(p);
  float d2 = primitiveB(p);
  return smin( d1, d2 );
}

/**
 * Twist
 */
float opTwist( vec3 p ) {
  float c = cos(20.0*p.y);
  float s = sin(20.0*p.y);
  mat2  m = mat2(c,-s,s,c);
  vec3  q = vec3(m*p.xz,p.y);
  return primitive(q);
}

/**
 * Cheap Bend
 */
float opCheapBend( vec3 p ) {
  float c = cos(20.0*p.y);
  float s = sin(20.0*p.y);
  mat2  m = mat2(c,-s,s,c);
  vec3  q = vec3(m*p.xy,p.z);
  return primitive(q);
}

//////////////////////////////////////////////////
// http://www.iquilezles.org/www/articles/functions/functions.htm

/**
 * Almost Identity
 * Say you don't want to change a value unless it's too small and
 * screws some of your computations up. Then, rather than doing a
 * sharp conditional branch, you can blend your value with your
 * threshold, and do it smoothly (say, with a cubic polynomial).
 * Set m to be your threshold (anything above m stays unchanged),
 * and n the value things will take when your value is zero
 * Then:
 * p(0) = n
 * p(m) = m
 * p'(0) = 0
 * p'(m) = 1
 *
 * therefore,
 * if p(x) is a cubic,
 * then p(x) = (2n-m)(x/m)^3 + (2m-3n)(x/m)^2 + n
 */
float almostIdentity( float x, float m, float n ) {
  if( x>m ) return x;

  const float a = 2.0*n - m
  const float b = 2.0*m - 3.0*n;
  const float t = x/m;

  return (a*t + b)*t*t + n;
}

/**
 * Impulse
 * Great for triggering behaviours or making envelopes for music
 * or animation, and for anything that grows fast and then slowly
 * decays. Use k to control the stretching o the function. Btw,
 * it's maximum, which is 1.0, happens at exactly x = 1/k.
 */
float impulse( float k, float x ) {
  const float h = k*x;
  return h*exp(1.0-h);
}

/**
 * Cubic Pulse
 * Of course you found yourself doing
 * smoothstep(c-w,c,x)-smoothstep(c,c+w,x)
 * very often, probably cause you were trying to isolate some features.
 * Then this cubicPulse() is your friend. Also, why not, you can use it
 * as a cheap replacement for a gaussian.
 */
float cubicPulse( float c, float w, float x ) {
  x = fabs(x - c);
  if( x>w ) return 0.0;
  x /= w;
  return 1.0 - x*x*(3.0-2.0*x);
}

/**
 * Exponential Step
 * A natural attenuation is an exponential of a linearly decaying
 * quantity: yellow curve, exp(-x). A gaussian, is an exponential of a
 * quadratically decaying quantity: light green curve, exp(-x²). You can
 * go on increasing powers, and get a sharper and sharper smoothstep(),
 * until you get a step() in the limit.
 */
float expStep( float x, float k, float n ) {
  return exp( -k*pow(x,n) );
}

/**
 * Gain
 * Remapping the unit interval into the unit interval by expanding the
 * sides and compressing the center, and keeping 1/2 mapped to 1/2,
 * that can be done with the gain() function. This was a common function
 * in RSL tutorials (the Renderman Shading Language). k=1 is the identity
 * curve, k<1 produces the classic gain() shape, and k>1 produces "s"
 * shaped curces. The curves are symmetric (and inverse) for k=a and k=1/a.
 */
float gain(float x, float k) {
  float a = 0.5*pow(2.0*((x<0.5)?x:1.0-x), k);
  return (x<0.5)?a:1.0-a;
}

/**
 * Parabola
 * A nice choice to remap the 0..1 interval into 0..1, such that the
 * corners are remapped to 0 and the center to 1. In other words,
 * parabola(0) = parabola(1) = 0, and parabola(1/2) = 1.
 */
float parabola( float x, float k ) {
  return pow( 4.0*x*(1.0-x), k );
}

/**
 * Power curve
 * A nice choice to remap the 0..1 interval into 0..1, such that the
 * corners are remapped to 0. Very useful to skew the shape one side
 * or the other in order to make leaves, eyes, and many other
 * interesting shapes
 */
float pcurve( float x, float a, float b ) {
  float k = pow(a+b,a+b) / (pow(a,a)*pow(b,b));
  return k * pow( x, a ) * pow( 1.0-x, b );
}

/**
 * Sinc curve
 * A phase shifter sinc curve can be useful if it starts at zero and
 * ends at zero, for some bouncing behaviors (suggested by Hubert-Jan).
 * Give k different integer values to tweak the amount of bounces. It
 * peaks at 1.0, but that take negative values, which can make it
 * unusable in some applications.
 */
float sinc( float x, float k ) {
  const float a = PI * ((float(k)*x-1.0);
  return sin(a)/a;
}

//////////////////////////////////////////////////
// http://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm

/**
 * Sphere Intersection
 */
float sphIntersect( vec3 ro, vec3 rd, vec4 sph ) {
  vec3 oc = ro - sph.xyz;
  float b = dot( oc, rd );
  float c = dot( oc, oc ) - sph.w*sph.w;
  float h = b*b - c;
  if( h<0.0 ) return -1.0;
  h = sqrt( h );
  return -b - h;
}

/**
 * Sphere Density
 * Maths / Article: http://www.iquilezles.org/www/articles/spheredensity/spheredensity.htm
 * Example / Code: https://www.shadertoy.com/view/XljGDy
 */
float sphDensity( vec3 ro, vec3 rd, vec4 sph, float dbuffer ) {
  float ndbuffer = dbuffer/sph.w;
  vec3  rc = (ro - sph.xyz)/sph.w;

  float b = dot(rd,rc);
  float c = dot(rc,rc) - 1.0;
  float h = b*b - c;
  if( h<0.0 ) return 0.0;
  h = sqrt( h );
  float t1 = -b - h;
  float t2 = -b + h;

  if( t2<0.0 || t1>ndbuffer ) return 0.0;
  t1 = max( t1, 0.0 );
  t2 = min( t2, ndbuffer );

  float i1 = -(c*t1 + b*t1*t1 + t1*t1*t1/3.0);
  float i2 = -(c*t2 + b*t2*t2 + t2*t2*t2/3.0);
  return (i2-i1)*(3.0/4.0);
}

/**
 * Sphere Projection
 * Maths / Article: http://www.iquilezles.org/www/articles/sphereproj/sphereproj.htm
 * Example / Code: https://www.shadertoy.com/view/XdBGzd
 */
float projectSphere( vec4 sph, in mat4 cam, float fle ) {
  vec3  o = (cam*vec4(sph.xyz,1.0)).xyz;

  float r2 = sph.w*sph.w;
  float z2 = o.z*o.z;
  float l2 = dot(o,o);
  float k1 = l2-r2;
  float k2 = r2-z2;

  return -pi*fle*fle*r2*sqrt(abs(k1/k2))/k2;
}

/**
 * Sphere Visibility
 * Maths / Article: http://www.iquilezles.org/www/articles/sphereocc/sphereocc.htm
 * Example / Code: https://www.shadertoy.com/view/XdS3Rt
 */
int sphereVisibility( vec4 sA, vec4 sB, vec3 c ) {
  vec3 ac = sA.xyz - c;
  vec3 bc = sB.xyz - c;
  float ia = 1.0/length(ac);
  float ib = 1.0/length(bc);
  float k0 = dot(ac,bc)*ia*ib;
  float k1 = sA.w*ia;
  float k2 = sB.w*ib;
  float m1 = k0*k0 + k1*k1 + k2*k2;
  float m2 = 2.0*k0*k1*k2;

       if( m1 + m2 - 1.0 < 0.0 ) return 1;
  else if( m1 - m2 - 1.0 < 0.0 ) return 2;
  return 3;
}

/**
 * Sphere Ambient Occlusion
 * Maths / Article: http://www.iquilezles.org/www/articles/sphereao/sphereao.htm
 * Example / Code: https://www.shadertoy.com/view/4djSDy
 */
float sphOcclusion( vec3 pos, vec3 nor, vec4 sph ) {
  vec3  di = sph.xyz - pos;
  float l  = length(di);
  float nl = dot(nor,di/l);
  float h  = l/sph.w;
  float h2 = h*h;
  float k2 = 1.0 - h2*nl*nl;

  float res = max(0.0,nl)/h2;
  if( k2 > 0.0 ) { // approx. for penetration
    res = clamp(0.5*(nl*h+1.0)/h2,0.0,1.0);
    res = sqrt( res*res*res );
  }
  return res;
}

/**
 * Sphere Soft Shadow
 * Example / Code: https://www.shadertoy.com/view/4d2XWV
 */
float sphSoftShadow( vec3 ro, vec3 rd, vec4 sph, float k ) {
  vec3 oc = ro - sph.xyz;
  float r = sph.w*sph.w;
  float b = dot( oc, rd );
  float c = dot( oc, oc ) - r;
  float h = b*b - c;
  float d = -sph.w + sqrt( max(0.0,r-h));
  float t = -b     - sqrt( max(0.0,h) );
  return (t<0.0)?1.0:smoothstep(0.0,1.0,k*d/t);
}

//////////////////////////////////////////////////
// Export

// Palette
#pragma glslify: export(palette)

// Raymarching Distance
#pragma glslify: export(sdSphere)
#pragma glslify: export(udBox)
#pragma glslify: export(udRoundBox)
#pragma glslify: export(sdBox)
#pragma glslify: export(sdTorus)
#pragma glslify: export(sdCylinder)
#pragma glslify: export(sdCone)
#pragma glslify: export(sdPlane)
#pragma glslify: export(sdHexPrism)
#pragma glslify: export(sdTriPrism)
#pragma glslify: export(sdCapsule)
#pragma glslify: export(sdCappedCylinder)
#pragma glslify: export(sdCappedCone)
#pragma glslify: export(udTriangle)
#pragma glslify: export(udQuad)
#pragma glslify: export(sdTorus82)
#pragma glslify: export(sdTorus88)
#pragma glslify: export(opU)
#pragma glslify: export(opS)
#pragma glslify: export(opI)
#pragma glslify: export(opRep)
#pragma glslify: export(opTx)
#pragma glslify: export(opScale)
#pragma glslify: export(opDisplace)
#pragma glslify: export(opBlend)
#pragma glslify: export(opTwist)
#pragma glslify: export(opCheapBend)

// Useful functions
#pragma glslify: export(almostIdentity)
#pragma glslify: export(impulse)
#pragma glslify: export(cubicPulse)
#pragma glslify: export(expStep)
#pragma glslify: export(gain)
#pragma glslify: export(parabola)
#pragma glslify: export(pcurve)
#pragma glslify: export(sinc)
// Sphere
#pragma glslify: export(sphIntersect)
#pragma glslify: export(sphDensity)
#pragma glslify: export(projectSphere)
#pragma glslify: export(sphereVisibility)
#pragma glslify: export(sphOcclusion)
#pragma glslify: export(sphSoftShadow)
