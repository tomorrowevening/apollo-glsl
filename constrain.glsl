/**
 * - Circle / Rect
 */

bool inCircle(vec2 pos, float radius, vec2 pixel) {
  float halfRad = radius * 0.5;
  return distance(pixel, pos + halfRad) < halfRad;
}

bool inRect(vec2 pos, vec2 size, vec2 pixel) {
  return pixel.x > pos.x && pixel.y > pos.y && pixel.x < (size.x + pos.x) && pixel.y < (size.y + pos.y);
}

#pragma glslify: export(inCircle)
#pragma glslify: export(inRect)
