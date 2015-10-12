// Blobs
// Eric Galin

const int Steps = 1000;				//nb pas max
const float Epsilon = 0.05; // Marching epsilon
const float T=0.5;

//à modifier si les objets ne s'affichent pas
const float rA=10.0; // Maximum ray marching or sphere tracing distance from origin
const float rB=40.0; // Minimum

//***************************************//
//*****			Transformations		*****//
//***************************************//

vec3 rotateX(vec3 p, float a)
{
  float sa = sin(a);
  float ca = cos(a);
  return vec3(p.x, ca*p.y - sa*p.z, sa*p.y + ca*p.z);
}

vec3 rotateY(vec3 p, float a)
{
  float sa = sin(a);
  float ca = cos(a);
  return vec3(ca*p.x + sa*p.z, p.y, -sa*p.x + ca*p.z);
}

vec3 rotateZ(vec3 p, float a)
{
  float sa = sin(a);
  float ca = cos(a);
  return vec3(ca*p.x + sa*p.y, -sa*p.x + ca*p.y, p.z);
}

// Smooth falloff function
// r : small radius
// R : Large radius
float falloff( float r, float R ) //fonction G du cours
{
  float x = clamp(r/R,0.0,1.0);
  float y = (1.0-x*x);
  return y*y*y;
}

//***************************************//
//******	Primitive Operations	*****//
//***************************************//

// Blending
// a : field function of left sub-tree
// b : field function of right sub-tree
float Blend(float a,float b)
{
    return a+b;
}

// Union
// a : field function of left sub-tree
// b : field function of right sub-tree
float Union(float a,float b)
{
    return max(a,b);
}

// Intersection
// a : field function of left sub-tree
// b : field function of right sub-tree
float Intersection(float a,float b)
{
    return min(a,b);
}

// Différence
// a : field function of left sub-tree
// b : field function of right sub-tree
float Difference(float a,float b)
{
    return a-b;
}

//***************************************//
//******	Primitive functions		*****//
//***************************************//

// Point skeleton
// p : point
// c : center of skeleton
// e : energy associated to skeleton
// R : large radius
float point(vec3 p, vec3 c, float e,float R)
{
  return e*falloff(length(p-c),R);
}

// Segment
// p : point
// a : extrémité 1 du segment
// b : extrémité 2 du segment
// e : energy associated to skeleton
// R : segment radius
float seg(in vec3 p, in vec3 a, in vec3 b, float e, float R){
  float d;
  vec3 ab = (b-a)/length(b-a);
  if(dot(p-a,ab)<0.0){d=length(p-a);}
  else if(dot(p-b,ab)>0.0){d=length(p-b);}
  else{
	  float t = dot(p-a,ab);
	  d = length(a+t*ab-p);
  }
  return e*falloff(d,R);
}

// Cercle
// p : point
// c : centre du cercle
// n : direction du cercle
//rayon : rayon interne du cercle
// e : energy associated to skeleton
// R : segment radius
float cercle(vec3 p, vec3 c, vec3 n,float rayon, float e, float R)
{
    vec3 normal = normalize(n);
	float ph = dot(normal,p-c);
	float ch = sqrt(dot(p-c,p-c)-(ph*ph));
	float qh = ch - rayon;
	float d = sqrt((qh*qh)+(ph*ph));
	return e*falloff(d,R);
}

// Cube
// p : point
// c : center of cube
// e : energy associated to skeleton
// R : distance au cube
// cote : largeur cube
float cube(vec3 p, vec3 c, float e, float R, vec3 cote)
{
    p = p - c;
    cote = cote/2.0;
    p = abs(p);
    float val = 0.0;
    
    if(p.x > cote.x)
        val += (cote.x - p.x)*(cote.x - p.x);
    
    if(p.y > cote.y)
        val += (cote.y - p.y)*(cote.y - p.y);
    
    if(p.z > cote.z)
        val += (cote.z - p.z)*(cote.z - p.z);
    
    
    return e * falloff(val, R);
}

// Disque
// p : point
// c : centre du disque
// n : normal du disque
// h : hauteur du disque
// e : energy associated to skeleton
// R : segment radius
float disque(in vec3 p, in vec3 c, in vec3 n, float rayon, float e, float R){
  
  n = normalize(n);
  
  float ph = dot(n, p-c);
  float ch = sqrt( dot(p-c, p-c) - ph*ph );
  float d;
  
  if(ch < rayon){
	  d = abs(ph);
  }
  else{
	  float qh = ch - rayon;
	  d = sqrt(qh*qh + ph*ph);
  }
  return e*falloff(d,R);
}

// Potential field of the object
// p : point
float Humain(vec3 p)
{
  p.z=-p.z;

  float v = Blend(point(p, vec3(0.0, 0.0, 0.0), 1.0, 4.5),
                  point(p, vec3(0.0, 3.7, 0.0), 6.5, 2.0));
    
    //bras
  v=Blend(v,point(p,vec3(-4.0, 1.0,0.0),0.5,3.5));
  v=Blend(v,point(p,vec3(4.0, 1.0,0.0),0.5,3.5));
    
    //jambe
  v=Blend(v,point(p,vec3(-1.6, -3.0,0.0),0.5,2.5));
  v=Blend(v,point(p,vec3(1.6, -3.0,0.0),0.5,2.5));
    
  v=Blend(v,point(p,vec3(-2.0, -5.0,0.0),0.5,2.5));
  v=Blend(v,point(p,vec3(2.0, -5.0,0.0),0.5,2.5));
    
  return v;
}




// Potential field of the object
// p : point
float object(vec3 p) //c'est ici qu'on créer notre objet en faisant des unions, intersection etc
{
  p.z=-p.z; //pour afficher l'objet à l'endroit
  
  //cacahuète du prof
  /*
  float v = Blend(point(p,vec3( 0.0, 1.0, 1.0),1.0,4.5),
                  point(p,vec3( 2.0, 0.0,-3.0),1.0,4.5));

  v=Blend(v,point(p,vec3(-3.0, 2.0,-3.0),1.0,4.5));
  v=Union(v,point(p,vec3(-1.0, -1.0, 0.0),1.0,4.5));
  */

  //float v = Humain(p);
  
  /*float v = seg(p, vec3(0,0,0), vec3(3,0,0), 1.0, 1.0);  
  float x = cos(iGlobalTime*0.1*3.14), y = sin(iGlobalTime*0.1*3.14);
  v = Union(v, seg(p, vec3(0,0,0), vec3(x*3.0,y*3.0, 0), 1.0, 1.0));*/

  /*float v = cube(p, vec3(0.0, 2.0,3.0),3.0,vec3(4.,4.,4.));
  v = Blend(v,cube(p, vec3(1.0, 1.0,0.0),3.0,vec3(4.,4.,4.)));*/
  
  //float v = disque(p, vec3(0,0,0), vec3(1,0,0), 2.0, 1.0, 4.0 );

  float v = cercle(p, vec3(0.0, 0.0,0.0),vec3(0.5, 0.2,0.0),3.0,0.6,1.0);
  
  return v-T;
}


// Calculate object normal
// p : point
vec3 ObjectNormal(in vec3 p )
{
  float eps = 0.0001;
  vec3 n;
  float v = object(p);
  n.x = object( vec3(p.x+eps, p.y, p.z) ) - v;
  n.y = object( vec3(p.x, p.y+eps, p.z) ) - v;
  n.z = object( vec3(p.x, p.y, p.z+eps) ) - v;
  return normalize(n);
}

//***************************************//
//******		Tracing				*****//
//***************************************//

// Trace ray using ray marching
// o : ray origin
// u : ray direction
// h : hit
// s : Number of steps
float Trace(vec3 o, vec3 u, out bool h,out int s)
{
  h = false;

    // Don't start at the origin, instead move a little bit forward
    float t=rA;

  for(int i=0; i<Steps; i++)
  {
    s=i;
    vec3 p = o+t*u;
    float v = object(p);
    // Hit object
      if (v > 0.0)
      {
          s=i;
          h = true;
          break;
      }
      // Move along ray
      t += Epsilon;
      // Escape marched far away
      if (t>rB)
      {
          break;
      }
  }
  return t;
}

// Trace ray using ray marching
// o : ray origin
// u : ray direction
// h : hit
// s : Number of steps
float SphereTrace(vec3 o, vec3 u, out bool h,out int s)
{
  h = false;

    // Don't start at the origin, instead move a little bit forward
    float t=rA;

  for(int i=0; i<Steps; i++)
  {
    s=i;
    vec3 p = o+t*u;
    float v = object(p);
    // Hit object
      if (v > 0.0)
      {
          s=i;
          h = true;
          break;
      }
      // Move along ray
      t += max(Epsilon,abs(v)/4.0); //le 4 à modifier si il y a des erreurs de shading
      // Escape marched far away
      if (t>rB)
      {
          break;
      }
  }
  return t;
}


// Background color
vec3 background(vec3 rd)
{
  return mix(vec3(0.4, 0.3, 0.0), vec3(0.7, 0.8, 1.0), rd.y*0.5+0.5);
}

// Shading and lighting
// p : point,
// n : normal at point
vec3 Shade(vec3 p, vec3 n)
{
  // point light
  const vec3 lightPos = vec3(5.0, 5.0, 5.0);
  const vec3 lightColor = vec3(0.5, 0.5, 0.5);

  vec3 c = 0.25*background(n);
  vec3 l = normalize(lightPos - p);

  // Not even Phong shading, use weighted cosine instead for smooth transitions
  float diff = 0.5*(1.0+dot(n, l));

  c += diff*lightColor;

  return c;
}

// Shading with number of steps
vec3 ShadeSteps(int n)
{
   float t=float(n)/(float(Steps-1));
   return vec3(t,0.25+0.75*t,0.5-0.5*t);
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
  vec2 pixel = (gl_FragCoord.xy / iResolution.xy)*2.0-1.0;

  // compute ray origin and direction
  float asp = iResolution.x / iResolution.y;
  vec3 rd = normalize(vec3(asp*pixel.x, pixel.y, -4.0));
  vec3 ro = vec3(0.0, 0.0, 20.0);

  // vec2 mouse = iMouse.xy / iResolution.xy;
  float a=iGlobalTime*0.25;
  ro = rotateY(ro, a);
  rd = rotateY(rd, a);

  // Trace ray
  bool hit;

  // Number of steps
  int s;

  float t = SphereTrace(ro, rd, hit,s);
  vec3 pos=ro+t*rd;
  // Shade background
  vec3 rgb = background(rd);

  if (hit)
  {
    // Compute normal
    vec3 n = ObjectNormal(pos);

    // Shade object with light
    rgb = Shade(pos, n);
  }

  // Uncomment this line to shade image with false colors representing the number of steps
  //rgb = ShadeSteps(s);

  fragColor=vec4(rgb, 1.0);
}

