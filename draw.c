#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "math.h"
#include "gmath.h"
#include "symtab.h"
#include "hash.h"

void draw_phong(struct matrix *polygons, screen s, zbuffer zb,
                double *view, double light[2][3], color ambient,
                struct constants *reflect) {
  if (polygons->lastcol < 3) {
    printf("Need at least 3 points to draw a polygon!\n");
    return;
  }

  int point;
  double *n;
  for (point = 0; point < polygons->lastcol - 2; point += 3) {
    n = calculate_normal(polygons, point);
    insert(polygons->m[0][point],
	   polygons->m[1][point],
	   polygons->m[2][point],
	   n[0], n[1], n[2]);
    insert(polygons->m[0][point + 1],
	   polygons->m[1][point + 1],
	   polygons->m[2][point + 1],
	   n[0], n[1], n[2]);
    insert(polygons->m[0][point + 2],
	   polygons->m[1][point + 2],
	   polygons->m[2][point + 2],
	   n[0], n[1], n[2]);
  }

  int v;
  for (point = 0; point < polygons->lastcol - 2; point += 3) {
    n = calculate_normal(polygons, point);
    v = dot_product(n, view);
    if (v > 0) {
      interpolate_phong(polygons, point, s, zb,
		  view, light, ambient,
		  reflect);
	
    }
  }
}
void interpolate_phong(struct matrix *points, int point, screen s, zbuffer zb,
                 double *view, double light[2][3], color ambient,
                 struct constants *reflect) {
  int bot, mid, top, y;
  int t_to_b, m_to_b, t_to_m;
  double x0, x1, dx0, dx1, y0, y1, y2, z0, z1, dz0, dz1;
  double dn0[3], dn1[3];
  double *b_n, *m_n, *t_n; 
  int flip = 0;
  y0 = points->m[1][point];
  y1 = points->m[1][point + 1];
  y2 = points->m[1][point + 2];

  if (y0<=y1 && y0<=y2) {
    bot = point;
    if (y1 <= y2) {
      mid = point + 1;
      top = point + 2;
    }
    else {
      mid = point + 2;
      top = point + 1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = point+1;
    if (y0 <= y2) {
      mid = point;
      top = point + 2;
    }
    else {
      mid = point + 2;
      top = point;
    }
  }//end y1 bottom
  else {
    bot = point + 2;
    if (y0 <= y1) {
      mid = point;
      top = point + 1;
    }
    else {
      mid = point + 1;
      top = point;
    }
  }

  b_n = search(points->m[0][bot], points->m[1][bot], points->m[2][bot]);
  m_n = search(points->m[0][mid], points->m[1][mid], points->m[2][mid]);
  t_n = search(points->m[0][top], points->m[1][top], points->m[2][top]);
  double n0[3] = {b_n[0], b_n[1], b_n[2]};
  double n1[3] = {b_n[0], b_n[1], b_n[2]};
  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  y = (int)(points->m[1][bot]);

  int t_x, t_y, t_z, m_x, m_y, m_z, b_x, b_y, b_z;
  t_x = (int)(points->m[0][top]);
  m_x = (int)(points->m[0][mid]);
  b_x = (int)(points->m[0][bot]);
  t_y = (int)(points->m[1][top]);
  m_y = (int)(points->m[1][mid]);
  b_y = (int)(points->m[1][bot]);
  t_z = (int)(points->m[2][top]);
  m_z = (int)(points->m[2][mid]);
  b_z = (int)(points->m[2][bot]);


  t_to_b = t_y - y; 
  m_to_b = m_y - y; 
  t_to_m = t_y- m_y;

  dx0 = t_to_b > 0 ? (t_x - b_x) / t_to_b : 0;
  dx1 = m_to_b > 0 ? (m_x - b_x) / m_to_b : 0;
  dz0 = t_to_b > 0 ? (t_z - b_z) / t_to_b : 0;
  dz1 = m_to_b > 0 ? (m_z - b_z) / m_to_b : 0;
  dn0[0] = t_to_b > 0 ? (double)(t_n[0] - b_n[0]) / t_to_b : 0; 
  dn0[1] = t_to_b > 0 ? (double)(t_n[1] - b_n[1]) / t_to_b : 0; 
  dn0[2] = t_to_b > 0 ? (double)(t_n[2] - b_n[2]) / t_to_b : 0; 
  dn1[0] = m_to_b > 0 ? (double)(m_n[0] - b_n[0]) / m_to_b : 0; 
  dn1[1] = m_to_b > 0 ? (double)(m_n[1] - b_n[1]) / m_to_b : 0; 
  dn1[2] = m_to_b > 0 ? (double)(m_n[2] - b_n[2]) / m_to_b : 0; 

  while (y <= t_y) {
    double n0_temp[3], n1_temp[3];
    n0_temp[0] = n0[0];
    n0_temp[1] = n0[1];
    n0_temp[2] = n0[2];
    n1_temp[0] = n1[0];
    n1_temp[1] = n1[1];
    n1_temp[2] = n1[2];
    draw_line_phong(x0, y, z0, x1, y, z1, s, zb, n0_temp, n1_temp, view, light, ambient, reflect);

    x0 += dx0;
    x1 += dx1;
    z0 += dz0;
    z1 += dz1;

    n0[0] += dn0[0];
    n0[1] += dn0[1];
    n0[2] += dn0[2];
    n1[0] += dn1[0];
    n1[1] += dn1[1];
    n1[2] += dn1[2];

    y ++;
    if (!flip && y >= m_y) {
      //printf("Flipping\n");
      flip = 1;
      dx1 = t_to_m > 0 ? (t_x - m_x) / t_to_m : 0;
      dz1 = t_to_m > 0 ? (t_z - m_z) / t_to_m : 0;
      x1 = m_x;
      z1 = m_z;
      n1[0] = m_n[0];
      n1[1] = m_n[1];
      n1[2] = m_n[2];
      dn1[0] = t_to_m > 0 ? (t_n[0] - m_n[0]) * 1.0 / t_to_m : 0;
      dn1[1] = t_to_m > 0 ? (t_n[1] - m_n[1]) * 1.0 / t_to_m : 0;
      dn1[2] = t_to_m > 0 ? (t_n[2] - m_n[2]) * 1.0 / t_to_m : 0;
    }
  }
}

void draw_gouraud(struct matrix *polygons, screen s, zbuffer zb,
                  double *view, double light[2][3], color ambient,
                  struct constants *reflect) {
  if (polygons->lastcol < 3) {
    printf("Need at least 3 points to draw a polygon!\n");
    return;
  }

  int point;
  double *normal;
  //	create_hash(polygons->lastcol);

  for (point = 0; point < polygons->lastcol - 2; point += 3) {
    normal = calculate_normal(polygons, point);
    insert(polygons->m[0][point], polygons->m[1][point], polygons->m[2][point],
	   normal[0], normal[1], normal[2]);
    insert(polygons->m[0][point + 1], polygons->m[1][point + 1], polygons->m[2][point + 1],
	   normal[0], normal[1], normal[2]);
    insert(polygons->m[0][point + 2], polygons->m[1][point + 2], polygons->m[2][point + 2],
	   normal[0], normal[1], normal[2]);
  }
  for (point = 0; point < polygons->lastcol - 2; point += 3) {
    normal = calculate_normal(polygons, point);
    if (dot_product(normal, view) > 0) {

      interpolate_gouraud(polygons, point, s, zb, view, light, ambient, reflect);

      double *n0 = search(polygons->m[0][point], polygons->m[1][point], polygons->m[2][point]);
      double *n1 = search(polygons->m[0][point + 1], polygons->m[1][point + 1], polygons->m[2][point + 1]);
      double *n2 = search(polygons->m[0][point + 2], polygons->m[1][point + 2], polygons->m[2][point + 2]);
      color c0 = get_lighting(n0, view, ambient, light, reflect);
      color c1 = get_lighting(n1, view, ambient, light, reflect);
      color c2 = get_lighting(n2, view, ambient, light, reflect);
		
      draw_line_gouraud(polygons->m[0][point],
			polygons->m[1][point],
			polygons->m[2][point],
			polygons->m[0][point+1],
			polygons->m[1][point+1],
			polygons->m[2][point+1],
			s, zb, c0, c1);
      draw_line_gouraud(polygons->m[0][point+2],
			polygons->m[1][point+2],
			polygons->m[2][point+2],
			polygons->m[0][point+1],
			polygons->m[1][point+1],
			polygons->m[2][point+1],
			s, zb, c2, c1);
      draw_line_gouraud(polygons->m[0][point],
			polygons->m[1][point],
			polygons->m[2][point],
			polygons->m[0][point+2],
			polygons->m[1][point+2],
			polygons->m[2][point+2],
			s, zb, c0, c2);			
    }
  }
}


void interpolate_gouraud(struct matrix *points, int i, screen s, zbuffer zb,
                   double *view, double light[2][3], color ambient,
                   struct constants *reflect) {
  int bot, mid, top, y;
  int distance0, distance1, distance2;
  double x0, x1, dx0, dx1, y0, y1, y2, z0, z1, dz0, dz1;
  double *b_n, *m_n, *t_n; // Vertex normals
  color b_c, m_c, t_c, c0, c1; // Associated colors
  double c0_exact[3], c1_exact[3], dc0[3], dc1[3]; // Remember to round to ints
  int flip = 0;
  y0 = points->m[1][i];
  y1 = points->m[1][i + 1];
  y2 = points->m[1][i + 2];
  if (y0 <= y1 && y0 <= y2) {
    bot = i;
    if (y1 <= y2) {
      mid = i + 1;
      top = i + 2;
    }
    else {
      mid = i + 2;
      top = i + 1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = i+1;
    if (y0 <= y2) {
      mid = i;
      top = i + 2;
    }
    else {
      mid = i + 2;
      top = i;
    }
  }//end y1 bottom
  else {
    bot = i + 2;
    if (y0 <= y1) {
      mid = i;
      top = i + 1;
    }
    else {
      mid = i + 1;
      top = i;
    }
  }

  b_n = search(points->m[0][bot], points->m[1][bot], points->m[2][bot]);
  m_n = search(points->m[0][mid], points->m[1][mid], points->m[2][mid]);
  t_n = search(points->m[0][top], points->m[1][top], points->m[2][top]);

  b_c = get_lighting(b_n, view, ambient, light, reflect);
  m_c = get_lighting(m_n, view, ambient, light, reflect);
  t_c = get_lighting(t_n, view, ambient, light, reflect);
  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  c0.red = c1.red = b_c.red;
  c0.green = c1.green = b_c.green;
  c0.blue = c1.blue = b_c.blue;

  c0_exact[0] = c1_exact[0] = b_c.red;
  c0_exact[1] = c1_exact[1] = b_c.green;
  c0_exact[2] = c1_exact[2] = b_c.blue;
  y = (int)(points->m[1][bot]);

  distance0 = (int)(points->m[1][top]) - y; 
  distance1 = (int)(points->m[1][mid]) - y; 
  distance2 = (int)(points->m[1][top]) - (int)(points->m[1][mid]); 

  dx0 = distance0 > 0 ? (points->m[0][top] - points->m[0][bot]) / distance0 : 0;
  dx1 = distance1 > 0 ? (points->m[0][mid] - points->m[0][bot]) / distance1 : 0;
  dz0 = distance0 > 0 ? (points->m[2][top] - points->m[2][bot]) / distance0 : 0;
  dz1 = distance1 > 0 ? (points->m[2][mid] - points->m[2][bot]) / distance1 : 0;
  dc0[0] = distance0 > 0 ? (t_c.red - b_c.red) * 1.0 / distance0 : 0;
  dc0[1] = distance0 > 0 ? (t_c.green - b_c.green) * 1.0 / distance0 : 0;
  dc0[2] = distance0 > 0 ? (t_c.blue - b_c.blue) * 1.0 / distance0 : 0;
  dc1[0] = distance1 > 0 ? (m_c.red - b_c.red) * 1.0 / distance1 : 0;
  dc1[1] = distance1 > 0 ? (m_c.green - b_c.green) * 1.0 / distance1 : 0;
  dc1[2] = distance1 > 0 ? (m_c.blue - b_c.blue) * 1.0 / distance1 : 0;

  while (y <= (int)points->m[1][top]) {
    draw_line_gouraud(x0, y, z0, x1, y, z1, s, zb, c0, c1);
    x0 += dx0;
    x1 += dx1;
    z0 += dz0;
    z1 += dz1;
    c0_exact[0] += dc0[0];
    c0_exact[1] += dc0[1];
    c0_exact[2] += dc0[2];
    c1_exact[0] += dc1[0];
    c1_exact[1] += dc1[1];
    c1_exact[2] += dc1[2];
    c0.red = (int)(c0_exact[0]);
    c0.green = (int)(c0_exact[1]);
    c0.blue = (int)(c0_exact[2]);
    c1.red = (int)(c1_exact[0]);
    c1.green = (int)(c1_exact[1]);
    c1.blue = (int)(c1_exact[2]);
    y ++;
    if (!flip && y >= (int)(points->m[1][mid])) {
      //printf("Flipping\n");
      flip = 1;
      dx1 = distance2 > 0 ? (points->m[0][top] - points->m[0][mid]) / distance2 : 0;
      dz1 = distance2 > 0 ? (points->m[2][top] - points->m[2][mid]) / distance2 : 0;
      x1 = points->m[0][mid];
      z1 = points->m[2][mid];
      c1_exact[0] = c1.red = m_c.red;
      c1_exact[1] = c1.green = m_c.green;
      c1_exact[2] = c1.blue = m_c.blue;
      dc1[0] = distance2 > 0 ? (t_c.red - m_c.red) * 1.0 / distance2 : 0;
      dc1[1] = distance2 > 0 ? (t_c.green - m_c.green) * 1.0 / distance2 : 0;
      dc1[2] = distance2 > 0 ? (t_c.blue - m_c.blue) * 1.0 / distance2 : 0;
    }
  }
}


/*======== void shade_flat() ==========
  Inputs: struct matrix *points
  int i
  screen s
  zbuffer zb
  Returns:

  Fills in polygon i by drawing consecutive horizontal (or vertical) lines.

  Color should be set differently for each polygon.
  ====================*/
void scanline_convert(struct matrix *points, int i, screen s, zbuffer zb, color c) {
  int top, mid, bot, y;
  int distance0, distance1, distance2;
  double x0, x1, y0, y1, y2, dx0, dx1, z0, z1, dz0, dz1;
  int flip = 0;

  z0 = z1 = dz0 = dz1 = 0;

  y0 = points->m[1][i];
  y1 = points->m[1][i+1];
  y2 = points->m[1][i+2];

  //find bot, mid, top
  if ( y0 <= y1 && y0 <= y2) {
    bot = i;
    if (y1 <= y2) {
      mid = i+1;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i+1;
    }
  }//end y0 bottom
  else if (y1 <= y0 && y1 <= y2) {
    bot = i+1;
    if (y0 <= y2) {
      mid = i;
      top = i+2;
    }
    else {
      mid = i+2;
      top = i;
    }
  }//end y1 bottom
  else {
    bot = i+2;
    if (y0 <= y1) {
      mid = i;
      top = i+1;
    }
    else {
      mid = i+1;
      top = i;
    }
  }
  x0 = points->m[0][bot];
  x1 = points->m[0][bot];
  z0 = points->m[2][bot];
  z1 = points->m[2][bot];
  y = (int)(points->m[1][bot]);

  distance0 = (int)(points->m[1][top]) - y;
  distance1 = (int)(points->m[1][mid]) - y;
  distance2 = (int)(points->m[1][top]) - (int)(points->m[1][mid]);

  dx0 = distance0 > 0 ? (points->m[0][top]-points->m[0][bot])/distance0 : 0;
  dx1 = distance1 > 0 ? (points->m[0][mid]-points->m[0][bot])/distance1 : 0;
  dz0 = distance0 > 0 ? (points->m[2][top]-points->m[2][bot])/distance0 : 0;
  dz1 = distance1 > 0 ? (points->m[2][mid]-points->m[2][bot])/distance1 : 0;

  while ( y <= (int)points->m[1][top] ) {
    draw_line(x0, y, z0, x1, y, z1, s, zb, c);

    x0+= dx0;
    x1+= dx1;
    z0+= dz0;
    z1+= dz1;
    y++;

    if ( !flip && y >= (int)(points->m[1][mid]) ) {
      flip = 1;
      dx1 = distance2 > 0 ? (points->m[0][top]-points->m[0][mid])/distance2 : 0;
      dz1 = distance2 > 0 ? (points->m[2][top]-points->m[2][mid])/distance2 : 0;
      x1 = points->m[0][mid];
      z1 = points->m[2][mid];
    }//end flip code
  }//end scanline loop
}

/*======== void add_polygon() ==========
  Inputs:   struct matrix *surfaces
  double x0
  double y0
  double z0
  double x1
  double y1
  double z1
  double x2
  double y2
  double z2
  Returns:
  Adds the vertices (x0, y0, z0), (x1, y1, z1)
  and (x2, y2, z2) to the polygon matrix. They
  define a single triangle surface.
  ====================*/
void add_polygon( struct matrix *polygons,
		  double x0, double y0, double z0,
		  double x1, double y1, double z1,
		  double x2, double y2, double z2 ) {

  add_point(polygons, x0, y0, z0);
  add_point(polygons, x1, y1, z1);
  add_point(polygons, x2, y2, z2);
}

/*======== void draw_flat() ==========
  Inputs:   struct matrix *polygons
  screen s
  color c
  Returns:
  Goes through polygons 3 points at a time, drawing
  lines connecting each points to create bounding
  triangles
  ====================*/
void draw_polygons(struct matrix *polygons, screen s, zbuffer zb,
		   double *view, double light[2][3], color ambient,
		   struct constants *reflect) {
  if (polygons->lastcol < 3) {
    printf("Need at least 3 points to draw a polygon!\n");
    return;
  }

  int point;
  double *normal;

  for (point = 0; point < polygons->lastcol - 2; point += 3) {
    normal = calculate_normal(polygons, point);
    if (dot_product(normal, view) > 0) {
      color c = get_lighting(normal, view, ambient, light, reflect);
      scanline_convert(polygons, point, s, zb, c);
      draw_line(polygons->m[0][point],
		polygons->m[1][point],
		polygons->m[2][point],
		polygons->m[0][point+1],
		polygons->m[1][point+1],
		polygons->m[2][point+1],
		s, zb, c);
      draw_line(polygons->m[0][point+2],
		polygons->m[1][point+2],
		polygons->m[2][point+2],
		polygons->m[0][point+1],
		polygons->m[1][point+1],
		polygons->m[2][point+1],
		s, zb, c);
      draw_line(polygons->m[0][point],
		polygons->m[1][point],
		polygons->m[2][point],
		polygons->m[0][point+2],
		polygons->m[1][point+2],
		polygons->m[2][point+2],
		s, zb, c);
    }
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * edges
  double x
  double y
  double z
  double width
  double height
  double depth
  Returns:

  add the points for a rectagular prism whose
  upper-left corner is (x, y, z) with width,
  height and depth dimensions.
  ====================*/
void add_box( struct matrix * polygons,
	      double x, double y, double z,
	      double width, double height, double depth ) {

  double x1, y1, z1;
  x1 = x+width;
  y1 = y-height;
  z1 = z-depth;

  //front
  add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
  add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);

  //back
  add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
  add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

  //right side
  add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
  add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
  //left side
  add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
  add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

  //top
  add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
  add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
  //bottom
  add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
  add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);
}//end add_box

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r
  double step
  Returns:

  adds all the points for a sphere with center
  (cx, cy, cz) and radius r.

  should call generate_sphere to create the
  necessary points
  ====================*/
void add_sphere( struct matrix * edges,
		 double cx, double cy, double cz,
		 double r, int step ) {

  struct matrix *points = generate_sphere(cx, cy, cz, r, step);

  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 0;
  longStop = step;

  step++;
  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {

      p0 = lat * (step) + longt;
      p1 = p0+1;
      p2 = (p1+step) % (step * (step-1));
      p3 = (p0+step) % (step * (step-1));

      //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
      if (longt < step - 2)
	add_polygon( edges, points->m[0][p0],
		     points->m[1][p0],
		     points->m[2][p0],
		     points->m[0][p1],
		     points->m[1][p1],
		     points->m[2][p1],
		     points->m[0][p2],
		     points->m[1][p2],
		     points->m[2][p2]);
      if (longt > 0 )
	add_polygon( edges, points->m[0][p0],
		     points->m[1][p0],
		     points->m[2][p0],
		     points->m[0][p2],
		     points->m[1][p2],
		     points->m[2][p2],
		     points->m[0][p3],
		     points->m[1][p3],
		     points->m[2][p3]);
    }
  }
  free_matrix(points);
}

/*======== void generate_sphere() ==========
  Inputs:   double cx
  double cy
  double cz
  double r
  int step
  Returns: Generates all the points along the surface
  of a sphere with center (cx, cy, cz) and
  radius r.
  Returns a matrix of those points
  ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
				double r, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle <= circ_stop; circle++){
      circ = (double)circle / step;

      x = r * cos(M_PI * circ) + cx;
      y = r * sin(M_PI * circ) *
	cos(2*M_PI * rot) + cy;
      z = r * sin(M_PI * circ) *
	sin(2*M_PI * rot) + cz;

      /* printf("rotation: %d\tcircle: %d\n", rotation, circle); */
      /* printf("rot: %lf\tcirc: %lf\n", rot, circ); */
      /* printf("sphere point: (%0.2f, %0.2f, %0.2f)\n\n", x, y, z); */
      add_point(points, x, y, z);
    }
  }

  return points;
}

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r1
  double r2
  double step
  Returns:

  adds all the points required to make a torus
  with center (cx, cy, cz) and radii r1 and r2.

  should call generate_torus to create the
  necessary points
  ====================*/
void add_torus( struct matrix * edges,
		double cx, double cy, double cz,
		double r1, double r2, int step ) {

  struct matrix *points = generate_torus(cx, cy, cz, r1, r2, step);

  int p0, p1, p2, p3, lat, longt;
  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = step;
  longStart = 0;
  longStop = step;

  //printf("points: %d\n", points->lastcol);

  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {
      p0 = lat * step + longt;
      if (longt == step - 1)
	p1 = p0 - longt;
      else
	p1 = p0 + 1;
      p2 = (p1 + step) % (step * step);
      p3 = (p0 + step) % (step * step);

      //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
      add_polygon( edges, points->m[0][p0],
		   points->m[1][p0],
		   points->m[2][p0],
		   points->m[0][p3],
		   points->m[1][p3],
		   points->m[2][p3],
		   points->m[0][p2],
		   points->m[1][p2],
		   points->m[2][p2]);
      add_polygon( edges, points->m[0][p0],
		   points->m[1][p0],
		   points->m[2][p0],
		   points->m[0][p2],
		   points->m[1][p2],
		   points->m[2][p2],
		   points->m[0][p1],
		   points->m[1][p1],
		   points->m[2][p1]);
    }
  }
  free_matrix(points);
}
/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double cz
  double r
  int step
  Returns: Generates all the points along the surface
  of a torus with center (cx, cy, cz) and
  radii r1 and r2.
  Returns a matrix of those points
  ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
				double r1, double r2, int step ) {

  struct matrix *points = new_matrix(4, step * step);
  int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
  double x, y, z, rot, circ;

  rot_start = 0;
  rot_stop = step;
  circ_start = 0;
  circ_stop = step;

  for (rotation = rot_start; rotation < rot_stop; rotation++) {
    rot = (double)rotation / step;

    for(circle = circ_start; circle < circ_stop; circle++){
      circ = (double)circle / step;

      x = cos(2*M_PI * rot) *
	(r1 * cos(2*M_PI * circ) + r2) + cx;
      y = r1 * sin(2*M_PI * circ) + cy;
      z = -1*sin(2*M_PI * rot) *
	(r1 * cos(2*M_PI * circ) + r2) + cz;

      //printf("rotation: %d\tcircle: %d\n", rotation, circle);
      //printf("torus point: (%0.2f, %0.2f, %0.2f)\n", x, y, z);
      add_point(points, x, y, z);
    }
  }
  return points;
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * points
  double cx
  double cy
  double r
  double step
  Returns:

  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle( struct matrix * edges,
		 double cx, double cy, double cz,
		 double r, int step ) {
  double x0, y0, x1, y1, t;
  int i;
  x0 = r + cx;
  y0 = cy;

  for (i=1; i<=step; i++) {
    t = (double)i/step;
    x1 = r * cos(2 * M_PI * t) + cx;
    y1 = r * sin(2 * M_PI * t) + cy;

    add_edge(edges, x0, y0, cz, x1, y1, cz);
    x0 = x1;
    y0 = y1;
  }
}


/*======== void add_curve() ==========
  Inputs:   struct matrix *points
  double x0
  double y0
  double x1
  double y1
  double x2
  double y2
  double x3
  double y3
  double step
  int type
  Returns:

  Adds the curve bounded by the 4 points passsed as parameters
  of type specified in type (see matrix.h for curve type constants)
  to the matrix points
  ====================*/
void add_curve( struct matrix *edges,
		double x0, double y0,
		double x1, double y1,
		double x2, double y2,
		double x3, double y3,
		int step, int type ) {

  double t, x, y;
  struct matrix *xcoefs;
  struct matrix *ycoefs;
  int i;

  xcoefs = generate_curve_coefs(x0, x1, x2, x3, type);
  ycoefs = generate_curve_coefs(y0, y1, y2, y3, type);

  for (i=1; i<=step; i++) {

    t = (double)i/step;
    x = xcoefs->m[0][0] *t*t*t + xcoefs->m[1][0] *t*t+
      xcoefs->m[2][0] *t + xcoefs->m[3][0];
    y = ycoefs->m[0][0] *t*t*t + ycoefs->m[1][0] *t*t+
      ycoefs->m[2][0] *t + ycoefs->m[3][0];

    add_edge(edges, x0, y0, 0, x, y, 0);
    x0 = x;
    y0 = y;
  }

  free_matrix(xcoefs);
  free_matrix(ycoefs);
}


/*======== void add_point() ==========
  Inputs:   struct matrix * points
  int x
  int y
  int z
  Returns:
  adds point (x, y, z) to points and increment points.lastcol
  if points is full, should call grow on points
  ====================*/
void add_point( struct matrix * points, double x, double y, double z) {

  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][ points->lastcol ] = x;
  points->m[1][ points->lastcol ] = y;
  points->m[2][ points->lastcol ] = z;
  points->m[3][ points->lastcol ] = 1;
  points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
  Inputs:   struct matrix * points
  int x0, int y0, int z0, int x1, int y1, int z1
  Returns:
  add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
  should use add_point
  ====================*/
void add_edge( struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
  Inputs:   struct matrix * points
  screen s
  color c
  Returns:
  Go through points 2 at a time and call draw_line to add that line
  to the screen
  ====================*/
void draw_lines(struct matrix * points, screen s, zbuffer zb, color c) {

  if ( points->lastcol < 2 ) {
    printf("Need at least 2 points to draw a line!\n");
    return;
  }
  int point;
  for (point=0; point < points->lastcol-1; point+=2)
    draw_line( points->m[0][point],
	       points->m[1][point],
	       points->m[2][point],
	       points->m[0][point+1],
	       points->m[1][point+1],
	       points->m[2][point+1],
	       s, zb, c);
}// end draw_lines




void draw_line(int x0, int y0, double z0,
               int x1, int y1, double z1,
               screen s, zbuffer zb, color c) {


  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  int loop_start, loop_end;
  double distance;
  double z, dz;

  //swap points if going right -> left
  int xt, yt;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    z = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    x1 = xt;
    y1 = yt;
    z1 = z;
  }

  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) { //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;

  while ( loop_start < loop_end ) {

    plot( s, zb, c, x, y, z );
    if ( (wide && ((A > 0 && d > 0) ||
		   (A < 0 && d < 0)))
	 ||
	 (tall && ((A > 0 && d < 0 ) ||
		   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    loop_start++;
  } //end drawing loop
  plot( s, zb, c, x1, y1, z );
} //end draw_line

void draw_line_gouraud(int x0, int y0, double z0,
		       int x1, int y1, double z1,
		       screen s, zbuffer zb, color c0, color c1) {

  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  int loop_start, loop_end;
  double distance;
  double z, dz;

  //swap points if going right -> left
  int xt, yt;
  color ct;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    z = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    x1 = xt;
    y1 = yt;
    z1 = z;
    ct.red = c0.red;
    ct.green = c0.green;
    ct.blue = c0.blue;
    c0.red = c1.red;
    c0.green = c1.green;
    c0.blue = c1.blue;
    c1.red = ct.red;
    c1.green = ct.green;
    c1.blue = ct.blue;
  }

  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ){ //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;

  color c;
  c.red = c0.red;
  c.green = c0.green;
  c.blue = c0.blue;

  double cred = c0.red, cgreen = c0.green, cblue = c0.blue;
  double dred, dgreen, dblue;
  // printf("870\n");

  if (loop_end != loop_start){
    dred = (c1.red - c0.red) * 1.0 / (loop_end - loop_start);
    dgreen = (c1.green - c0.green) * 1.0 / (loop_end - loop_start);
    dblue = (c1.blue - c0.blue) * 1.0 / (loop_end - loop_start);
  }

  while (loop_start < loop_end ) {
    plot( s, zb, c, x, y, z );
    cred += dred;
    cgreen += dgreen;
    cblue += dblue;
    c.red = (int)cred;
    c.green = (int)cgreen;
    c.blue = (int)cblue;

    if ( (wide && ((A > 0 && d > 0) ||
		   (A < 0 && d < 0)))
	 ||
	 (tall && ((A > 0 && d < 0 ) ||
		   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    loop_start++;
  } //end drawing loop
  plot( s, zb, c1, x1, y1, z );
} //end draw_line

void draw_line_phong(int x0, int y0, double z0,
		     int x1, int y1, double z1,
		     screen s, zbuffer zb, double normal0[3], double normal1[3],
		     double *view, double light[2][3], color ambient,
		     struct constants *reflect) {

  int x, y, d, A, B;
  int dy_east, dy_northeast, dx_east, dx_northeast, d_east, d_northeast;
  int loop_start, loop_end;
  double distance;
  double z, dz;

  //swap points if going right -> left
  int xt, yt;
  double normal_t[3];
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    z = z0;
    x0 = x1;
    y0 = y1;
    z0 = z1;
    x1 = xt;
    y1 = yt;
    z1 = z;
    normal_t[0] = normal0[0];
    normal_t[1] = normal0[1];
    normal_t[2] = normal0[2];
    normal0[0] = normal1[0];
    normal0[1] = normal1[1];
    normal0[2] = normal1[2];
    normal1[0] = normal_t[0];
    normal1[1] = normal_t[1];
    normal1[2] = normal_t[2];
  }
  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);
  int wide = 0;
  int tall = 0;
  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ){ //octant 1/8
    wide = 1;
    loop_start = x;
    loop_end = x1;
    dx_east = dx_northeast = 1;
    dy_east = 0;
    d_east = A;
    distance = x1 - x;
    if ( A > 0 ) { //octant 1
      d = A + B/2;
      dy_northeast = 1;
      d_northeast = A + B;
    }
    else { //octant 8
      d = A - B/2;
      dy_northeast = -1;
      d_northeast = A - B;
    }
  }//end octant 1/8
  else { //octant 2/7
    tall = 1;
    dx_east = 0;
    dx_northeast = 1;
    distance = abs(y1 - y);
    if ( A > 0 ) {     //octant 2
      d = A/2 + B;
      dy_east = dy_northeast = 1;
      d_northeast = A + B;
      d_east = B;
      loop_start = y;
      loop_end = y1;
    }
    else {     //octant 7
      d = A/2 - B;
      dy_east = dy_northeast = -1;
      d_northeast = A - B;
      d_east = -1 * B;
      loop_start = y1;
      loop_end = y;
    }
  }

  z = z0;
  dz = (z1 - z0) / distance;

  double normal[3];

  normal[0] = normal0[0];
  normal[1] = normal0[1];
  normal[2] = normal0[2];

  double dnx, dny, dnz;

  if (loop_end != loop_start){
    dnx = (normal1[0] - normal0[0]) * 1.0 / (loop_end - loop_start);
    dny = (normal1[1] - normal0[1]) * 1.0 / (loop_end - loop_start);
    dnz = (normal1[2] - normal0[2]) * 1.0 / (loop_end - loop_start);
  }
	
  color c;
  while (loop_start < loop_end ) {
    double normal_temp[3];
    normal_temp[0] = normal[0];
    normal_temp[1] = normal[1];
    normal_temp[2] = normal[2];
    c = get_lighting(normal_temp, view, ambient, light, reflect);
    plot(s, zb, c, x, y, z );
    /*
      printf("normal: ");
      print_normal(normal);
      printf(" | dnormal: %0.2f %0.2f %0.2f | loop size: %d", dnx, dny, dnz, loop_end - loop_start);
      printf("\n");
    */
    normal[0] += dnx;
    normal[1] += dny;
    normal[2] += dnz;

    if ( (wide && ((A > 0 && d > 0) ||
		   (A < 0 && d < 0)))
	 ||
	 (tall && ((A > 0 && d < 0 ) ||
		   (A < 0 && d > 0) ))) {
      y+= dy_northeast;
      d+= d_northeast;
      x+= dx_northeast;
    }
    else {
      x+= dx_east;
      y+= dy_east;
      d+= d_east;
    }
    z+= dz;
    loop_start++;
  } //end drawing loop
  c = get_lighting(normal1, view, ambient, light, reflect);
  plot( s, zb, c, x1, y1, z );
} //end draw_line
