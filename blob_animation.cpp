const int Steps = 1000;
const float Epsilon = 0.05; // Marching epsilon
const float T=0.5;

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

// Smooth falloff function au carré
// r : small radius au carré
// R : Large radius
float falloff2( float r2, float R )
{
    float x2 = clamp(r2/(R*R),0.0,1.0);
    float y = (1.0-x2);		//x*x = (r/R)² == r²/R² == x2
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
    return min(a,2.0*T-b);
}

float Metamorphe(float a, float b)
{
    float time = sqrt((cos(iGlobalTime)+1.0)*0.5);
    return a*time+b*(1.0-time);
}

//***************************************//
//**************  WARP  *****************//
//***************************************//

float hash( float n )
{
    return fract(sin(n)*43758.5453123);
}

// 3d noise function reprit de "electron"
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
    vec3 pc = p-c;
  	return e*falloff2(dot(pc,pc),R);
}

// Segment
// p : point
// a : extrémité 1 du segment
// b : extrémité 2 du segment
// e : energy associated to skeleton
// R : segment radius
float seg(in vec3 p, in vec3 a, in vec3 b, float e, float R){
  	float d2;
  	vec3 ab = normalize(b-a);
    vec3 pa = p-a;
    float ph = dot(pa,ab);
    if(ph<=0.0)
        d2=dot(pa,pa);
    else
    {
     	vec3 pb = p-b;
        if(dot(pb,ab)>=0.0)
            d2=dot(pb,pb);
        else{
            vec3 tmp = a+ph*ab-p;	//vecteur de distance entre le segment ab et le point p
            d2 = dot(tmp,tmp);
        }
    }
  	return e*falloff2(d2,R);
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
    //vec3 normal = normalize(n);
    vec3 pc = p-c;
	float ph = abs(dot(n,pc));
    if(ph >= R)
        return 0.0;
    float pc2 = dot(pc,pc);
    if(pc2 >= (rayon+R)*(rayon+R))
        return 0.0;
	float ch = sqrt(pc2-(ph*ph));	//pythagore a² = hypothénuse²-b²
	float qh = ch - rayon;
	float d2 = (qh*qh)+(ph*ph);
	return e*falloff2(d2,R);
}

// Cube
// p : point
// c : center of cube
// e : energy associated to skeleton
// R : distance au cube
// cote : largeur cube
float cube(vec3 p, vec3 c, vec3 cote, float e, float R)
{
    p = abs(p - c);
    cote = cote/2.0;
    
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
  
  	//n = normalize(n);
  
    vec3 pc = p-c;
  	float ph = dot(n, pc);
  	float ch = sqrt(dot(pc, pc) - ph*ph);
  	
    float d;
  	if(ch < rayon){	//le point est en dessus ou en dessous du disque
	  	d = abs(ph);
  	}
    else{
        float qh = ch - rayon;
        d = sqrt(qh*qh + ph*ph);
  	}
 	return e*falloff(d,R);
}

//plus lent, tentative d'amélioration de la vitesse ratée.
float disque2(in vec3 p, in vec3 c, in vec3 n, float rayon, float e, float R){
    
    vec3 pc = p-c;
  	float ph = abs(dot(n, pc));
    if(ph >= R)
    	return 0.0;
    
  	float ch2 = dot(pc, pc) - ph*ph;
  	
  	if(ch2 <= rayon*rayon)	//le point est en dessus ou en dessous du disque
	  	return e*falloff(ph,R);
    else if(ch2 >= (rayon+R)*(rayon+R))
        return 0.0;
  	else{
        float qh = sqrt(ch2) - rayon; 	//le point est en dehors du rayon du disque
        if(qh >= R)
            return 0.0;
        float d2 = qh*qh + ph*ph;
 		return e*falloff2(d2,R);      
  	}
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
  	vec3 ab = normalize(b-a);
    
    vec3 pa = p-a;
	float ph = dot(pa,-ab);
	//en dessous du disque de centre a
    if(ph > 0.0){
        //algo du disque
		if(ph >= R)
			return 0.0;
        float ch2 = dot(pa, pa) - ph*ph;
        if(ch2 < rayon*rayon)
            return e*falloff(ph,R);
		else if(ch2 >= (rayon+R)*(rayon+R))
			return 0.0;
        else{
            float qh = sqrt(ch2) - rayon;
            float d2 = qh*qh + ph*ph;
			return e*falloff2(d2,R);
        }
    }
	else{
		vec3 pb = p-b;
		float ph = dot(pb, ab);
		//en dessus du disque de centre b
		if(ph>0.0){
            if(ph >= R)
				return 0.0;
			float ch2 = dot(pb, pb) - ph*ph;
			if(ch2 < rayon*rayon)
				return e*falloff(ph,R);
			else if(ch2 >= (rayon+R)*(rayon+R))
				return 0.0;
			else{
				float qh = sqrt(ch2) - rayon;
				float d2 = qh*qh + ph*ph;
				return e*falloff2(d2,R);
			}
		}
		//entre les deux disque
		else{
			vec3 tmp = b+ph*ab-p;
			float dist2 = dot(tmp,tmp);
			if(dist2 <= rayon*rayon)
				return e;
			else if(dist2 >= (rayon+R)*(rayon+R))
				return 0.0;
			float d = sqrt(dist2)-rayon;
			return e*falloff(d,R);
		}
	}
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
		
		d = Difference(d, cylindre( p, newA, newB, rayon/8.0, e, R/8.0 ) );
		
		
	}
    float d2 = disque(p, a-vec3(0,R*0.7,0), normalize(b-a),rayon*1.5, e, R);
    d2 = Blend(d2,disque(p, b+vec3(0,R*0.7,0), normalize(a-b),rayon*1.5, e, R));
    d2 = Blend(d2, cube(p, vec3(a.x, a.y-R, a.z), vec3( (rayon)*4., offsetCube, (rayon)*4. ), e, 0.1 ));
	d2 = Blend(d2, cube(p, vec3(b.x, b.y+R, b.z), vec3( (rayon)*4., offsetCube, (rayon)*4. ), e, 0.1 ));
	
    d = Union(d,d2);
    
	return d;
}

float vague(in vec3 p, in vec3 c, float rayon, float periode, float amplitude, float e, float R)
{
    float d = 0.0;
    vec3 pos = vec3(c.x,0,c.z)-vec3(p.x,0,p.z);
    float dist2 = dot(pos,pos);
    if(dist2 > (rayon+R)*(rayon+R))
        return 0.0;
    else if(dist2 > rayon*rayon)
        d += sqrt(dist2)-rayon;
    
    float val = sin(sqrt(dist2)/periode+iGlobalTime)*amplitude;
    //float val = cos(pos.x)+sin(pos.z);
    d += abs((p.y-c.y)-val);
    
    return e*falloff(d,R);
}

float bassin(vec3 p)
{
	//sol avec un trou
 	float v = cylindre(p, vec3(0,-5,0), vec3(0,-3,0), 10.0,1.0,1.0);
    v = Difference(v, cylindre(p, vec3(0,-4,0), vec3(0,-3,0), 5.0,1.0,1.0));
    
	//vague dans le trou
	v = Blend(v, vague(warp(vec3(p)), vec3(0,-3.5,0), 5.5, 0.15, 0.15,1.0,1.0));
    
	//piliers
    v = Union(v, colonne(p, vec3(8,-1.8,0), vec3(8,3.0,0), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(0,-1.8,8), vec3(0,3.0,8), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(-8,-1.8,0), vec3(-8,3.0,0), 0.5, 1.0, 0.6)); 
    v = Union(v, colonne(p, vec3(0,-1.8,-8), vec3(0,3.0,-8), 0.5, 1.0, 0.6));
	
	//toit
    v = Union(v, cylindre(p, vec3(0,4.2,0), vec3(0,5,0), 10.0,1.0,1.0));
    
    return v;
}

//pris sur le code d' "antialias"
vec3 hash3( float n ) { return fract(sin(vec3(n,n+1.0,n+2.0))*43758.5453123); }

float bouleMouv(vec3 p)
{
    float v = 0.0;
	float time = iGlobalTime;
		
    for(float i = 1.0;	i < 9.0;	i++)
    {
		float ra = pow(i/7.0,3.0);
	    vec3  pos = 1.0*cos( 6.2831*hash3(i*14.0) + 0.5*(1.0-0.7*ra)*hash3(i*7.0)*time );	
                      
    	v += point(p, pos*4.0, 1.0,3.0);
    }
    return v;
}

float bouleSnake(vec3 p)
{
    float v = 0.0;
	float time = iGlobalTime;
		
    for(float i = 1.0;	i < 16.0;	i++)
    {
        vec3 vi = vec3(i,i,i);
        float n = noise(vi)*i+time;
        float x = cos(n);
        float y = sin(n);
        float z = tan(n);
                      
    	v += point(p, vec3(x,y,z)*4.0, 1.0,2.0);
    }
    return v;
}

/*****************************************************************************************/

// Potential field of the object
// p : point
float object(vec3 p) //c'est ici qu'on créer notre objet en faisant des unions, intersection etc
{
    p.z=-p.z; //pour afficher l'objet à l'endroit    
    
    //cacahuète du prof
    float v;

    v = bouleMouv(p);   
    //v = vague(p, vec3(0,0,0), 15.0,0.3,0.2,1.0,1.0);
    //v = bouleSnake(p);
    //v = Metamorphe(testCercle(p),testSeg(p));
    
    return v-T;
}

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
      t += max(Epsilon,abs(v)/16.0); //le 4 à modifier si il y a des erreurs de shading
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



