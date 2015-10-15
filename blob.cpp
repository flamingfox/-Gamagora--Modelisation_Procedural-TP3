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
//**************  WARP  *****************//
//***************************************//

float hash( float n )
{
    return fract(sin(n)*43758.5453123);
}

// 3d noise function
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    float res = mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
    return res;
}

vec3 warp(vec3 p)
{

 	float val = noise(p*0.75);//*sin(iGlobalTime));
    val += noise(p*0.25);
    p += 0.3*sin(val+cos(iGlobalTime));//(sin(iGlobalTime+val)*sin(iGlobalTime))+(cos(iGlobalTime)*cos(iGlobalTime+val-0.1));//-cos(iGlobalTime*1.3);//vec3(val,val,val);
    p -= 0.3*cos(val+sin(iGlobalTime+46.0));
    return p;
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
// rayon : rayon du disque
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


// clou
// p : point
// a : extrémité 1 du segment
// b : extrémité 2 du segment
// e : energy associated to skeleton
// R : segment radius
float clou(in vec3 p, in vec3 a, in vec3 b, float r, float e, float R){
  float d;
  vec3 ab = (b-a)/length(b-a);
  if(dot(p-a,ab)<-r){d=length(p-a)-r;}
  else if(dot(p-b,ab)>r){d=length(p-b);}
  else{
	  float t = dot(p-a,ab);
	  d = length(a+t*ab-p);
  }
  return e*falloff(d,R);
}

// Cylindre
float cylindre(in vec3 p, in vec3 a, in vec3 b, float rayon, float e, float R){
  float d;
  vec3 ab = normalize(b-a);
  if(dot(p-a,ab)<0.0){
      float ph = dot(p-a,ab);
      float ch2 = dot(p-a, p-a) - ph*ph;
      if(ch2 < rayon*rayon){
          d = abs(ph);}
      else{
          float qh = sqrt(ch2) - rayon;
          d = sqrt(qh*qh + ph*ph);
      }
  }
  
  else if(dot(p-b,ab)>0.0){
      float ph = dot(p-b,ab);
      float ch2 = dot(p-b, p-b) - ph*ph;
      if(ch2 < rayon*rayon){
          d = abs(ph);}
      else{
          float qh = sqrt(ch2) - rayon;
          d = sqrt(qh*qh + ph*ph);
      }
  }
  else{
	  float t = dot(p-a,ab);
	  d = length(a+t*ab-p)-rayon;
  }
  return e*falloff(d,R);
}



// Colonne
float colonne(in vec3 p, in vec3 a, in vec3 b, float rayon, float e, float R){
	float d;
	float offset = 0.28;
	float offsetCube = 0.15;
	d = cylindre(p, a, b, rayon, e, R);
	
	for(float i=0.0;i<12.0; i++){
        vec3 pos = vec3(cos(i*3.14/6.0)*(rayon+rayon*offset), 0.0, sin(i*3.14/6.0)*(rayon+rayon*offset) );
        
        vec3 newA = pos + a;//vec3(a.x + cos(i*3.14/6.0)*(rayon+rayon*offset), a.y + 0.0, a.z + sin(i*3.14/6.0)*(rayon+rayon*offset) );
		vec3 newB = pos + b;//vec3(b.x + cos(i*3.14/6.0)*(rayon+rayon*offset), b.y - 0.0, b.z + sin(i*3.14/6.0)*(rayon+rayon*offset) );
		
		d = Difference(d, cylindre( p, newA, newB , rayon/8.0, e, R/8.0 ) );
		
		
	}
    float d2 = disque(p, a-vec3(0,R*0.7,0), normalize(b-a),rayon*1.5, e, R);
    d2 = Blend(d2,disque(p, b+vec3(0,R*0.7,0), normalize(a-b),rayon*1.5, e, R));
    d2 = Blend(d2, cube(p, vec3(a.x, a.y-R, a.z), e, 0.1, vec3( (rayon)*4., offsetCube, (rayon)*4. ) ) );
	d2 = Blend(d2, cube(p, vec3(b.x, b.y+R, b.z), e, 0.1, vec3( (rayon)*4., offsetCube, (rayon)*4. ) ) );
	
    d = Union(d,d2);
    
	return d;
}

float vague(in vec3 p, in vec3 c, float rayon, float periode, float amplitude, float e, float R)
{
    float d = 0.0;
    vec3 pos = vec3(c.x,0,c.z)-vec3(p.x,0,p.z);
    float dist = length(pos);
    if(dist > rayon)
        d += dist-rayon;
    
    float val = sin(dist/periode+iGlobalTime)*amplitude;
    //float val = cos(pos.x)+sin(pos.z);
    d += abs((p.y-c.y)-val);
    
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

float bassin(vec3 p)
{
 	float v = cylindre(p, vec3(0,-4,0), vec3(0,-3,0), 10.0,1.0,1.0);
    v = Difference(v, cylindre(p, vec3(0,-4,0), vec3(0,-3,0), 5.0,1.0,1.0));
    v = Blend(v, vague(warp(vec3(p)), vec3(0,-3.5,0), 5.0, 0.15, 0.15,1.0,1.0));
    

              
    v = Union(v, colonne(p, vec3(8,-1.8,0), vec3(8,3.0,0), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(0,-1.8,8), vec3(0,3.0,8), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(-8,-1.8,0), vec3(-8,3.0,0), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(0,-1.8,-8), vec3(0,3.0,-8), 0.5, 1.0, 0.6));
	
    v = Union(v, cylindre(p, vec3(0,4.2,0), vec3(0,5,0), 10.0,1.0,1.0));
    
    return v;
}

// Potential field of the object
// p : point
float object(vec3 p) //c'est ici qu'on créer notre objet en faisant des unions, intersection etc
{
  p.z=-p.z; //pour afficher l'objet à l'endroit
  //p = warp(vec3(p));
    
  //cacahuète du prof
  
  /*float v = Blend(point(p,vec3( 0.0, 1.0, 1.0),1.0,4.5),
                  point(p,vec3( 2.0, 0.0,-3.0),1.0,4.5));

  v=Blend(v,point(p,vec3(-3.0, 2.0,-3.0),1.0,4.5));
  v=Union(v,point(p,vec3(-1.0, -1.0, 0.0),1.0,4.5));*/
  

  //float v = Humain(p);
  
  /*float v = seg(p, vec3(0,0,0), vec3(3,0,0), 1.0, 1.0);  
  float x = cos(iGlobalTime*0.1*3.14), y = sin(iGlobalTime*0.1*3.14);
  v = Union(v, seg(p, vec3(0,0,0), vec3(x*3.0,y*3.0, 0), 1.0, 1.0));*/

  /*float v = cube(p, vec3(0.0, 2.0,3.0),3.0,vec3(4.,4.,4.));
  v = Blend(v,cube(p, vec3(1.0, 1.0,0.0),3.0,vec3(4.,4.,4.)));*/
  
  //float v = disque(p, vec3(0,0,0), vec3(1,0,0), 2.0, 1.0, 4.0 );

  //float v = cylindre(p, vec3(-3,0,0), vec3(3,0,0),2.0, 1.0, 1.0);  
  //float v = disque(p, vec3(0.0, 0.0,0.0),vec3(0.5, 0.2,0.0),3.0,0.6,1.0);
  
  //float v = colonne(p, vec3(0,-3.5,0), vec3(0,3.5,0), 1.0, 1.0, 1.0); 
  //  	v+= colonne(p, vec3(0,-3.5,5), vec3(0,3.5,5), 1.0, 1.0, 1.0); 
    
  	float v = bassin(p);
    return v-T;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


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




